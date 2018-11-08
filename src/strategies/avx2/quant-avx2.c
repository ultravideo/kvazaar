/*****************************************************************************
* This file is part of Kvazaar HEVC encoder.
*
* Copyright (C) 2013-2015 Tampere University of Technology and others (see
* COPYING file).
*
* Kvazaar is free software: you can redistribute it and/or modify it under
* the terms of the GNU Lesser General Public License as published by the
* Free Software Foundation; either version 2.1 of the License, or (at your
* option) any later version.
*
* Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along
* with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

/*
* \file
*/

#include "strategies/avx2/quant-avx2.h"

#if COMPILE_INTEL_AVX2 && defined X86_64
#include <immintrin.h>
#include <stdlib.h>

#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "rdo.h"
#include "scalinglist.h"
#include "strategies/generic/quant-generic.h"
#include "strategies/strategies-quant.h"
#include "strategyselector.h"
#include "tables.h"
#include "transform.h"

static INLINE int32_t reduce_mm256i(__m256i src)
{
  __m128i a = _mm256_extracti128_si256(src, 0);
  __m128i b = _mm256_extracti128_si256(src, 1);

  a = _mm_add_epi32(a, b);
  b = _mm_shuffle_epi32(a, _MM_SHUFFLE(0, 1, 2, 3));

  a = _mm_add_epi32(a, b);
  b = _mm_shuffle_epi32(a, _MM_SHUFFLE(2, 3, 0, 1));

  a = _mm_add_epi32(a, b);
  return _mm_cvtsi128_si32(a);
}

static INLINE int32_t reduce_16x_i16(__m256i src)
{
  __m128i a = _mm256_extracti128_si256(src, 0);
  __m128i b = _mm256_extracti128_si256(src, 1);
  __m256i c = _mm256_cvtepi16_epi32(a);
  __m256i d = _mm256_cvtepi16_epi32(b);

  c = _mm256_add_epi32(c, d);
  return reduce_mm256i(c);
}

// If ints is completely zero, returns 16 in *first and -1 in *last
static INLINE void get_first_last_nz_int16(__m256i ints, int32_t *first, int32_t *last)
{
  // Note that nonzero_bytes will always have both bytes set for a set word
  // even if said word only had one of its bytes set, because we're doing 16
  // bit wide comparisons. No big deal, just shift results to the right by one
  // bit to have the results represent indexes of first set words, not bytes.
  // Another note, it has to use right shift instead of division to preserve
  // behavior on an all-zero vector (-1 / 2 == 0, but -1 >> 1 == -1)
  const __m256i zero = _mm256_setzero_si256();

  __m256i zeros = _mm256_cmpeq_epi16(ints, zero);
  uint32_t nonzero_bytes = ~((uint32_t)_mm256_movemask_epi8(zeros));
  *first = (    (int32_t)_tzcnt_u32(nonzero_bytes)) >> 1;
  *last = (31 - (int32_t)_lzcnt_u32(nonzero_bytes)) >> 1;
}

// Rearranges a 16x32b double vector into a format suitable for a stable SIMD
// max algorithm:
// (abcd|efgh) (ijkl|mnop) => (aceg|ikmo) (bdfh|jlnp)
static INLINE void rearrange_512(__m256i *hi, __m256i *lo)
{
  __m256i tmphi, tmplo;

  tmphi = _mm256_shuffle_epi32(*hi, _MM_SHUFFLE(3, 1, 2, 0));
  tmplo = _mm256_shuffle_epi32(*lo, _MM_SHUFFLE(3, 1, 2, 0));

  tmphi = _mm256_permute4x64_epi64(tmphi, _MM_SHUFFLE(3, 1, 2, 0));
  tmplo = _mm256_permute4x64_epi64(tmplo, _MM_SHUFFLE(3, 1, 2, 0));

  *hi = _mm256_permute2x128_si256(tmplo, tmphi, 0x31);
  *lo = _mm256_permute2x128_si256(tmplo, tmphi, 0x20);
}

static INLINE void get_cheapest_alternative(__m256i costs_hi, __m256i costs_lo,
    __m256i ns, __m256i changes,
    int16_t *final_change, int32_t *min_pos)
{
  __m128i nslo, nshi, chlo, chhi;
  __m256i pllo, plhi; // Payload
  __m256i tmp1, tmp2;

  nshi = _mm256_extracti128_si256(ns, 1);
  nslo = _mm256_extracti128_si256(ns, 0);
  chhi = _mm256_extracti128_si256(changes, 1);
  chlo = _mm256_extracti128_si256(changes, 0);

  // Interleave ns and lo into 32-bit variables and to two 256-bit wide vecs,
  // to have the same data layout as in costs. Zero extend to 32b width, shift
  // changes 16 bits to the left, and store them into the same vectors.
  tmp1 = _mm256_cvtepu16_epi32(nslo);
  tmp2 = _mm256_cvtepu16_epi32(chlo);
  tmp2 = _mm256_bslli_epi128(tmp2, 2);
  pllo = _mm256_or_si256(tmp1, tmp2);

  tmp1 = _mm256_cvtepu16_epi32(nshi);
  tmp2 = _mm256_cvtepu16_epi32(chhi);
  tmp2 = _mm256_bslli_epi128(tmp2, 2);
  plhi = _mm256_or_si256(tmp1, tmp2);

  // Reorder to afford result stability (if multiple atoms tie for cheapest,
  // rightmost ie. the highest is the wanted one)
  rearrange_512(&costs_hi, &costs_lo);
  rearrange_512(&plhi, &pllo);

  // 0: pick hi, 1: pick lo (equality evaluates as 0)
  __m256i cmpmask1 = _mm256_cmpgt_epi32(costs_hi, costs_lo);
  __m256i cost1    = _mm256_blendv_epi8(costs_hi, costs_lo, cmpmask1);
  __m256i pl1      = _mm256_blendv_epi8(plhi, pllo, cmpmask1);

  __m256i cost2    = _mm256_shuffle_epi32(cost1, _MM_SHUFFLE(2, 3, 0, 1));
  __m256i pl2      = _mm256_shuffle_epi32(pl1,   _MM_SHUFFLE(2, 3, 0, 1));

  __m256i cmpmask2 = _mm256_cmpgt_epi32(cost2, cost1);
  __m256i cost3    = _mm256_blendv_epi8(cost2, cost1, cmpmask2);
  __m256i pl3      = _mm256_blendv_epi8(pl2,   pl1,   cmpmask2);

  __m256i cost4    = _mm256_shuffle_epi32(cost3, _MM_SHUFFLE(1, 0, 3, 2));
  __m256i pl4      = _mm256_shuffle_epi32(pl3,   _MM_SHUFFLE(1, 0, 3, 2));

  __m256i cmpmask3 = _mm256_cmpgt_epi32(cost4, cost3);
  __m256i cost5    = _mm256_blendv_epi8(cost4, cost3, cmpmask3);
  __m256i pl5      = _mm256_blendv_epi8(pl4,   pl3,   cmpmask3);

  __m256i cost6    = _mm256_permute4x64_epi64(cost5, _MM_SHUFFLE(1, 0, 3, 2));
  __m256i pl6      = _mm256_permute4x64_epi64(pl5,   _MM_SHUFFLE(1, 0, 3, 2));

  __m256i cmpmask4 = _mm256_cmpgt_epi32(cost6, cost5);
  __m256i pl7      = _mm256_blendv_epi8(pl6,   pl5,   cmpmask4);

  __m128i res128 = _mm256_castsi256_si128(pl7);
  uint32_t tmp = (uint32_t)_mm_extract_epi32(res128, 0);
  uint16_t n = (uint16_t)(tmp & 0xffff);
  uint16_t chng = (uint16_t)(tmp >> 16);

  *final_change = (int16_t)chng;
  *min_pos = (int32_t)n;
}

#define VEC_WIDTH 16
#define SCAN_SET_SIZE 16
#define LOG2_SCAN_SET_SIZE 4

static INLINE int32_t hide_block_sign(__m256i coefs, const coeff_t * __restrict q_coef_reord, __m256i deltas_h, __m256i deltas_l, coeff_t * __restrict q_coef, const uint32_t * __restrict scan, int32_t subpos, int32_t last_cg)
{
  // Ensure that the block is 256 bit (32 byte) aligned
  assert(subpos % (32 / sizeof(coeff_t)) == 0);
  assert(((size_t)q_coef_reord) % 32 == 0);
  assert(SCAN_SET_SIZE == 16);

  __m256i q_coefs = _mm256_load_si256((__m256i *)(q_coef_reord + subpos));

  int32_t first_nz_pos_in_cg, last_nz_pos_in_cg;
  int32_t abssum = 0;

  // Find first and last nonzero coeffs
  get_first_last_nz_int16(q_coefs, &first_nz_pos_in_cg, &last_nz_pos_in_cg);

  // Sum all kvz_quant coeffs between first and last
  abssum = reduce_16x_i16(q_coefs);

  if (last_nz_pos_in_cg >= 0 && last_cg == -1) {
    last_cg = 1;
  }

  if (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {

    uint32_t q_coef_signbits = _mm256_movemask_epi8(q_coefs);
    int32_t signbit = (q_coef_signbits >> (2 * first_nz_pos_in_cg + 1)) & 0x1;

    if (signbit != (abssum & 0x1)) { // compare signbit with sum_parity
      int32_t min_pos = -1;
      int16_t final_change = 0;

      const int32_t mask_max = (last_cg == 1) ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1;

      const __m256i zero = _mm256_setzero_si256();
      const __m256i ones = _mm256_set1_epi16(1);
      const __m256i maxiters = _mm256_set1_epi16(mask_max);
      const __m256i ff = _mm256_set1_epi8(0xff);

      const __m256i fnpics = _mm256_set1_epi16((int16_t)first_nz_pos_in_cg);
      const __m256i ns = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

      __m256i block_signbit = _mm256_set1_epi16(((int16_t)signbit) * -1);
      __m256i coef_signbits = _mm256_cmpgt_epi16(zero, coefs);
      __m256i signbits_equal_block = _mm256_cmpeq_epi16(coef_signbits, block_signbit);

      __m256i q_coefs_zero = _mm256_cmpeq_epi16(q_coefs, zero);

      __m256i dus_packed = _mm256_packs_epi32(deltas_l, deltas_h);
      __m256i dus_ordered = _mm256_permute4x64_epi64(dus_packed, _MM_SHUFFLE(3, 1, 2, 0));
      __m256i dus_positive = _mm256_cmpgt_epi16(dus_ordered, zero);

      __m256i q_coef_abss = _mm256_abs_epi16(q_coefs);
      __m256i q_coefs_plusminus_one = _mm256_cmpeq_epi16(q_coef_abss, ones);

      __m256i eq_fnpics = _mm256_cmpeq_epi16(fnpics, ns);
      __m256i lt_fnpics = _mm256_cmpgt_epi16(fnpics, ns);

      __m256i maxcost_subcond1s = _mm256_and_si256(eq_fnpics, q_coefs_plusminus_one);
      __m256i maxcost_subcond2s = _mm256_andnot_si256(signbits_equal_block, lt_fnpics);
      __m256i elsecond1s_inv = _mm256_or_si256(dus_positive, maxcost_subcond1s);
      __m256i elsecond1s = _mm256_andnot_si256(elsecond1s_inv, ff);

      __m256i outside_maxiters = _mm256_cmpgt_epi16(ns, maxiters);

      __m256i negdelta_cond1s = _mm256_andnot_si256(q_coefs_zero, dus_positive);
      __m256i negdelta_cond2s = _mm256_andnot_si256(maxcost_subcond2s, q_coefs_zero);
      __m256i negdelta_mask16s_part1 = _mm256_or_si256(negdelta_cond1s, negdelta_cond2s);
      __m256i negdelta_mask16s = _mm256_andnot_si256(outside_maxiters, negdelta_mask16s_part1);

      __m256i posdelta_mask16s_part1 = _mm256_andnot_si256(q_coefs_zero, elsecond1s);
      __m256i posdelta_mask16s = _mm256_andnot_si256(outside_maxiters, posdelta_mask16s_part1);

      __m256i maxcost_cond1_parts = _mm256_andnot_si256(dus_positive, maxcost_subcond1s);
      __m256i maxcost_cond1s = _mm256_andnot_si256(q_coefs_zero, maxcost_cond1_parts);
      __m256i maxcost_cond2s = _mm256_and_si256(q_coefs_zero, maxcost_subcond2s);
      __m256i maxcost_mask16s_parts = _mm256_or_si256(maxcost_cond1s, maxcost_cond2s);
      __m256i maxcost_mask16s = _mm256_or_si256(maxcost_mask16s_parts, outside_maxiters);

      __m128i tmp_l, tmp_h;
      tmp_l = _mm256_extracti128_si256(negdelta_mask16s, 0);
      tmp_h = _mm256_extracti128_si256(negdelta_mask16s, 1);
      __m256i negdelta_mask32s_l = _mm256_cvtepi16_epi32(tmp_l);
      __m256i negdelta_mask32s_h = _mm256_cvtepi16_epi32(tmp_h);

      tmp_l = _mm256_extracti128_si256(posdelta_mask16s, 0);
      tmp_h = _mm256_extracti128_si256(posdelta_mask16s, 1);
      __m256i posdelta_mask32s_l = _mm256_cvtepi16_epi32(tmp_l);
      __m256i posdelta_mask32s_h = _mm256_cvtepi16_epi32(tmp_h);

      tmp_l = _mm256_extracti128_si256(maxcost_mask16s, 0);
      tmp_h = _mm256_extracti128_si256(maxcost_mask16s, 1);
      __m256i maxcost_mask32s_l = _mm256_cvtepi16_epi32(tmp_l);
      __m256i maxcost_mask32s_h = _mm256_cvtepi16_epi32(tmp_h);

      // Output value generation
      // cur_change_max: zero
      // cur_change_negdelta: ff
      // cur_change_posdelta: ones
      __m256i costs_negdelta_h = _mm256_sub_epi32(zero, deltas_h);
      __m256i costs_negdelta_l = _mm256_sub_epi32(zero, deltas_l);
      // costs_posdelta_l and _h: deltas_l and _h
      __m256i costs_max_lh = _mm256_set1_epi32(0x7fffffff);

      __m256i change_neg = _mm256_and_si256(negdelta_mask16s, ones);
      __m256i change_pos = _mm256_and_si256(posdelta_mask16s, ff);
      __m256i change_max = _mm256_and_si256(maxcost_mask16s, zero);

      __m256i cost_neg_l = _mm256_and_si256(negdelta_mask32s_l, costs_negdelta_l);
      __m256i cost_neg_h = _mm256_and_si256(negdelta_mask32s_h, costs_negdelta_h);
      __m256i cost_pos_l = _mm256_and_si256(posdelta_mask32s_l, deltas_l);
      __m256i cost_pos_h = _mm256_and_si256(posdelta_mask32s_h, deltas_h);
      __m256i cost_max_l = _mm256_and_si256(maxcost_mask32s_l, costs_max_lh);
      __m256i cost_max_h = _mm256_and_si256(maxcost_mask32s_h, costs_max_lh);

      __m256i changes = _mm256_or_si256(change_neg, _mm256_or_si256(change_pos, change_max));
      __m256i costs_l = _mm256_or_si256(cost_neg_l, _mm256_or_si256(cost_pos_l, cost_max_l));
      __m256i costs_h = _mm256_or_si256(cost_neg_h, _mm256_or_si256(cost_pos_h, cost_max_h));

      get_cheapest_alternative(costs_h, costs_l, ns, changes, &final_change, &min_pos);

      coeff_t cheapest_q = q_coef_reord[min_pos + subpos];
      if (cheapest_q == 32767 || cheapest_q == -32768)
        final_change = -1;

      uint32_t coef_signs = _mm256_movemask_epi8(coef_signbits);
      uint32_t cheapest_coef_sign_mask = (uint32_t)(1 << (2 * min_pos));

      if (!(coef_signs & cheapest_coef_sign_mask))
        q_coef[scan[min_pos + subpos]] += final_change;
      else
        q_coef[scan[min_pos + subpos]] -= final_change;
    } // Hide
  }
  if (last_cg == 1)
    last_cg = 0;

  return last_cg;
}

/**
 * \brief quantize transformed coefficents
 *
 */
void kvz_quant_flat_avx2(const encoder_state_t * const state, coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t * const scan = kvz_g_sig_last_scan[scan_idx][log2_block_size - 1];

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6);
  const uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  const int32_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  const int32_t q_bits = QUANT_SHIFT + qp_scaled / 6 + transform_shift;
  const int32_t add = ((state->frame->slicetype == KVZ_SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;

  assert(quant_coeff[0] <= (1 << 15) - 1 && quant_coeff[0] >= -(1 << 15)); //Assuming flat values to fit int16_t

  uint32_t ac_sum = 0;
  int32_t last_cg = -1;

  __m256i v_ac_sum = _mm256_setzero_si256();
  __m256i v_quant_coeff = _mm256_set1_epi16(quant_coeff[0]);

  for (int32_t n = 0; n < width * height; n += VEC_WIDTH) {

    __m256i v_level = _mm256_loadu_si256((__m256i*)&(coef[n]));
    __m256i v_sign = _mm256_cmpgt_epi16(_mm256_setzero_si256(), v_level);
    v_sign = _mm256_or_si256(v_sign, _mm256_set1_epi16(1));

    v_level = _mm256_abs_epi16(v_level);
    __m256i low_a = _mm256_unpacklo_epi16(v_level, _mm256_set1_epi16(0));
    __m256i high_a = _mm256_unpackhi_epi16(v_level, _mm256_set1_epi16(0));

    __m256i low_b = _mm256_unpacklo_epi16(v_quant_coeff, _mm256_set1_epi16(0));
    __m256i high_b = _mm256_unpackhi_epi16(v_quant_coeff, _mm256_set1_epi16(0));

    __m256i v_level32_a = _mm256_madd_epi16(low_a, low_b);
    __m256i v_level32_b = _mm256_madd_epi16(high_a, high_b);

    v_level32_a = _mm256_add_epi32(v_level32_a, _mm256_set1_epi32(add));
    v_level32_b = _mm256_add_epi32(v_level32_b, _mm256_set1_epi32(add));

    v_level32_a = _mm256_srai_epi32(v_level32_a, q_bits);
    v_level32_b = _mm256_srai_epi32(v_level32_b, q_bits);

    v_level = _mm256_packs_epi32(v_level32_a, v_level32_b);
    v_level = _mm256_sign_epi16(v_level, v_sign);

    _mm256_storeu_si256((__m256i*)&(q_coef[n]), v_level);

    v_ac_sum = _mm256_add_epi32(v_ac_sum, v_level32_a);
    v_ac_sum = _mm256_add_epi32(v_ac_sum, v_level32_b);
  }

  __m128i temp = _mm_add_epi32(_mm256_castsi256_si128(v_ac_sum), _mm256_extracti128_si256(v_ac_sum, 1));
  temp = _mm_add_epi32(temp, _mm_shuffle_epi32(temp, _MM_SHUFFLE(1, 0, 3, 2)));
  temp = _mm_add_epi32(temp, _mm_shuffle_epi32(temp, _MM_SHUFFLE(0, 1, 0, 1)));
  ac_sum += _mm_cvtsi128_si32(temp);

  if (!encoder->cfg.signhide_enable || ac_sum < 2)
    return;

  coeff_t coef_reord[LCU_WIDTH * LCU_WIDTH >> 2] ALIGNED(32);
  coeff_t q_coef_reord[LCU_WIDTH * LCU_WIDTH >> 2] ALIGNED(32);

  // Reorder coef and q_coef for sequential access
  for (int32_t n = 0; n < width * height; n++) {
    coef_reord[n] = coef[scan[n]];
    q_coef_reord[n] = q_coef[scan[n]];
  }

  assert(VEC_WIDTH == SCAN_SET_SIZE);
  for (int32_t subpos = (width * height - 1) & (~(VEC_WIDTH - 1)); subpos >= 0; subpos -= VEC_WIDTH) {

    __m256i v_coef = _mm256_load_si256((__m256i *)(coef_reord + subpos));
    __m256i v_level = _mm256_abs_epi16(v_coef);
    __m256i low_a = _mm256_unpacklo_epi16(v_level, _mm256_set1_epi16(0));
    __m256i high_a = _mm256_unpackhi_epi16(v_level, _mm256_set1_epi16(0));

    __m256i low_b = _mm256_unpacklo_epi16(v_quant_coeff, _mm256_set1_epi16(0));
    __m256i high_b = _mm256_unpackhi_epi16(v_quant_coeff, _mm256_set1_epi16(0));

    __m256i v_level32_a = _mm256_madd_epi16(low_a, low_b);
    __m256i v_level32_b = _mm256_madd_epi16(high_a, high_b);

    v_level32_a = _mm256_add_epi32(v_level32_a, _mm256_set1_epi32(add));
    v_level32_b = _mm256_add_epi32(v_level32_b, _mm256_set1_epi32(add));

    v_level32_a = _mm256_srai_epi32(v_level32_a, q_bits);
    v_level32_b = _mm256_srai_epi32(v_level32_b, q_bits);

    v_level = _mm256_packs_epi32(v_level32_a, v_level32_b);

    __m256i v_coef_a = _mm256_unpacklo_epi16(_mm256_abs_epi16(v_coef), _mm256_set1_epi16(0));
    __m256i v_coef_b = _mm256_unpackhi_epi16(_mm256_abs_epi16(v_coef), _mm256_set1_epi16(0));
    __m256i v_quant_coeff_a = _mm256_unpacklo_epi16(v_quant_coeff, _mm256_set1_epi16(0));
    __m256i v_quant_coeff_b = _mm256_unpackhi_epi16(v_quant_coeff, _mm256_set1_epi16(0));
    v_coef_a = _mm256_madd_epi16(v_coef_a, v_quant_coeff_a);
    v_coef_b = _mm256_madd_epi16(v_coef_b, v_quant_coeff_b);
    v_coef_a = _mm256_sub_epi32(v_coef_a, _mm256_slli_epi32(_mm256_unpacklo_epi16(v_level, _mm256_set1_epi16(0)), q_bits) );
    v_coef_b = _mm256_sub_epi32(v_coef_b, _mm256_slli_epi32(_mm256_unpackhi_epi16(v_level, _mm256_set1_epi16(0)), q_bits) );
    v_coef_a = _mm256_srai_epi32(v_coef_a, q_bits8);
    v_coef_b = _mm256_srai_epi32(v_coef_b, q_bits8);
    
    __m256i deltas_h = _mm256_permute2x128_si256(v_coef_a, v_coef_b, 0x31);
    __m256i deltas_l = _mm256_permute2x128_si256(v_coef_a, v_coef_b, 0x20);

    last_cg = hide_block_sign(v_coef, q_coef_reord, deltas_h, deltas_l, q_coef, scan, subpos, last_cg);
  }

#undef VEC_WIDTH
#undef SCAN_SET_SIZE
#undef LOG2_SCAN_SET_SIZE
}

static INLINE __m128i get_residual_4x1_avx2(const kvz_pixel *a_in, const kvz_pixel *b_in){
  __m128i a = _mm_cvtsi32_si128(*(int32_t*)a_in);
  __m128i b = _mm_cvtsi32_si128(*(int32_t*)b_in);
  __m128i diff = _mm_sub_epi16(_mm_cvtepu8_epi16(a), _mm_cvtepu8_epi16(b) );
  return diff;
}

static INLINE __m128i get_residual_8x1_avx2(const kvz_pixel *a_in, const kvz_pixel *b_in){
  __m128i a = _mm_cvtsi64_si128(*(int64_t*)a_in);
  __m128i b = _mm_cvtsi64_si128(*(int64_t*)b_in);
  __m128i diff = _mm_sub_epi16(_mm_cvtepu8_epi16(a), _mm_cvtepu8_epi16(b) );
  return diff;
}

static INLINE int32_t get_quantized_recon_4x1_avx2(int16_t *residual, const kvz_pixel *pred_in){
  __m128i res = _mm_loadl_epi64((__m128i*)residual);
  __m128i pred = _mm_cvtsi32_si128(*(int32_t*)pred_in);
  __m128i rec = _mm_add_epi16(res, _mm_cvtepu8_epi16(pred));
  return _mm_cvtsi128_si32(_mm_packus_epi16(rec, rec));
}

static INLINE int64_t get_quantized_recon_8x1_avx2(int16_t *residual, const kvz_pixel *pred_in){
  __m128i res = _mm_loadu_si128((__m128i*)residual);
  __m128i pred = _mm_cvtsi64_si128(*(int64_t*)pred_in);
  __m128i rec = _mm_add_epi16(res, _mm_cvtepu8_epi16(pred));
  return _mm_cvtsi128_si64(_mm_packus_epi16(rec, rec));
}

static void get_residual_avx2(const kvz_pixel *ref_in, const kvz_pixel *pred_in, int16_t *residual, int width, int in_stride){

  __m128i diff = _mm_setzero_si128();
  switch (width) {
    case 4:
      diff = get_residual_4x1_avx2(ref_in + 0 * in_stride, pred_in + 0 * in_stride);
      _mm_storel_epi64((__m128i*)&(residual[0]), diff);
      diff = get_residual_4x1_avx2(ref_in + 1 * in_stride, pred_in + 1 * in_stride);
      _mm_storel_epi64((__m128i*)&(residual[4]), diff);
      diff = get_residual_4x1_avx2(ref_in + 2 * in_stride, pred_in + 2 * in_stride);
      _mm_storel_epi64((__m128i*)&(residual[8]), diff);
      diff = get_residual_4x1_avx2(ref_in + 3 * in_stride, pred_in + 3 * in_stride);
      _mm_storel_epi64((__m128i*)&(residual[12]), diff);
    break;
    case 8:
      diff = get_residual_8x1_avx2(&ref_in[0 * in_stride], &pred_in[0 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[0]), diff);
      diff = get_residual_8x1_avx2(&ref_in[1 * in_stride], &pred_in[1 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[8]), diff);
      diff = get_residual_8x1_avx2(&ref_in[2 * in_stride], &pred_in[2 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[16]), diff);
      diff = get_residual_8x1_avx2(&ref_in[3 * in_stride], &pred_in[3 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[24]), diff);
      diff = get_residual_8x1_avx2(&ref_in[4 * in_stride], &pred_in[4 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[32]), diff);
      diff = get_residual_8x1_avx2(&ref_in[5 * in_stride], &pred_in[5 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[40]), diff);
      diff = get_residual_8x1_avx2(&ref_in[6 * in_stride], &pred_in[6 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[48]), diff);
      diff = get_residual_8x1_avx2(&ref_in[7 * in_stride], &pred_in[7 * in_stride]);
      _mm_storeu_si128((__m128i*)&(residual[56]), diff);
    break;
    default:
      for (int y = 0; y < width; ++y) {
        for (int x = 0; x < width; x+=16) {
          diff = get_residual_8x1_avx2(&ref_in[x + y * in_stride], &pred_in[x + y * in_stride]);
          _mm_storeu_si128((__m128i*)&residual[x + y * width], diff);
          diff = get_residual_8x1_avx2(&ref_in[(x+8) + y * in_stride], &pred_in[(x+8) + y * in_stride]);
          _mm_storeu_si128((__m128i*)&residual[(x+8) + y * width], diff);
        }
      }
    break;
  }
}

static void get_quantized_recon_avx2(int16_t *residual, const kvz_pixel *pred_in, int in_stride, kvz_pixel *rec_out, int out_stride, int width){

  switch (width) {
    case 4:
      *(int32_t*)&(rec_out[0 * out_stride]) = get_quantized_recon_4x1_avx2(residual + 0 * width, pred_in + 0 * in_stride);
      *(int32_t*)&(rec_out[1 * out_stride]) = get_quantized_recon_4x1_avx2(residual + 1 * width, pred_in + 1 * in_stride);
      *(int32_t*)&(rec_out[2 * out_stride]) = get_quantized_recon_4x1_avx2(residual + 2 * width, pred_in + 2 * in_stride);
      *(int32_t*)&(rec_out[3 * out_stride]) = get_quantized_recon_4x1_avx2(residual + 3 * width, pred_in + 3 * in_stride);
      break;
    case 8:
      *(int64_t*)&(rec_out[0 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 0 * width, pred_in + 0 * in_stride);
      *(int64_t*)&(rec_out[1 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 1 * width, pred_in + 1 * in_stride);
      *(int64_t*)&(rec_out[2 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 2 * width, pred_in + 2 * in_stride);
      *(int64_t*)&(rec_out[3 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 3 * width, pred_in + 3 * in_stride);
      *(int64_t*)&(rec_out[4 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 4 * width, pred_in + 4 * in_stride);
      *(int64_t*)&(rec_out[5 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 5 * width, pred_in + 5 * in_stride);
      *(int64_t*)&(rec_out[6 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 6 * width, pred_in + 6 * in_stride);
      *(int64_t*)&(rec_out[7 * out_stride]) = get_quantized_recon_8x1_avx2(residual + 7 * width, pred_in + 7 * in_stride);
      break;
    default:
      for (int y = 0; y < width; ++y) {
        for (int x = 0; x < width; x += 16) {
          *(int64_t*)&(rec_out[x + y * out_stride]) = get_quantized_recon_8x1_avx2(residual + x + y * width, pred_in + x + y  * in_stride);
          *(int64_t*)&(rec_out[(x + 8) + y * out_stride]) = get_quantized_recon_8x1_avx2(residual + (x + 8) + y * width, pred_in + (x + 8) + y  * in_stride);
        }
      }
      break;
  }
}

/**
* \brief Quantize residual and get both the reconstruction and coeffs.
*
* \param width  Transform width.
* \param color  Color.
* \param scan_order  Coefficient scan order.
* \param use_trskip  Whether transform skip is used.
* \param stride  Stride for ref_in, pred_in and rec_out.
* \param ref_in  Reference pixels.
* \param pred_in  Predicted pixels.
* \param rec_out  Reconstructed pixels.
* \param coeff_out  Coefficients used for reconstruction of rec_out.
*
* \returns  Whether coeff_out contains any non-zero coefficients.
*/
int kvz_quantize_residual_avx2(encoder_state_t *const state,
  const cu_info_t *const cur_cu, const int width, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const kvz_pixel *const ref_in, const kvz_pixel *const pred_in,
  kvz_pixel *rec_out, coeff_t *coeff_out)
{
  // Temporary arrays to pass data to and from kvz_quant and transform functions.
  int16_t residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  coeff_t coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

  int has_coeffs = 0;

  assert(width <= TR_MAX_WIDTH);
  assert(width >= TR_MIN_WIDTH);

  // Get residual. (ref_in - pred_in -> residual)
  get_residual_avx2(ref_in, pred_in, residual, width, in_stride);

  // Transform residual. (residual -> coeff)
  if (use_trskip) {
    kvz_transformskip(state->encoder_control, residual, coeff, width);
  }
  else {
    kvz_transform2d(state->encoder_control, residual, coeff, width, color, cur_cu->type);
  }

  // Quantize coeffs. (coeff -> coeff_out)
  if (state->encoder_control->cfg.rdoq_enable &&
      (width > 4 || !state->encoder_control->cfg.rdoq_skip))
  {
    int8_t tr_depth = cur_cu->tr_depth - cur_cu->depth;
    tr_depth += (cur_cu->part_size == SIZE_NxN ? 1 : 0);
    kvz_rdoq(state, coeff, coeff_out, width, width, (color == COLOR_Y ? 0 : 2),
      scan_order, cur_cu->type, tr_depth);
  } else {
    kvz_quant(state, coeff, coeff_out, width, width, (color == COLOR_Y ? 0 : 2),
      scan_order, cur_cu->type);
  }

  // Check if there are any non-zero coefficients.
  for (int i = 0; i < width * width; i += 8) {
    __m128i v_quant_coeff = _mm_loadu_si128((__m128i*)&(coeff_out[i]));
    has_coeffs = !_mm_testz_si128(_mm_set1_epi8(0xFF), v_quant_coeff);
    if(has_coeffs) break;
  }

  // Do the inverse quantization and transformation and the reconstruction to
  // rec_out.
  if (has_coeffs) {

    // Get quantized residual. (coeff_out -> coeff -> residual)
    kvz_dequant(state, coeff_out, coeff, width, width, (color == COLOR_Y ? 0 : (color == COLOR_U ? 2 : 3)), cur_cu->type);
    if (use_trskip) {
      kvz_itransformskip(state->encoder_control, residual, coeff, width);
    }
    else {
      kvz_itransform2d(state->encoder_control, residual, coeff, width, color, cur_cu->type);
    }

    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    get_quantized_recon_avx2(residual, pred_in, in_stride, rec_out, out_stride, width);
  }
  else if (rec_out != pred_in) {
    // With no coeffs and rec_out == pred_int we skip copying the coefficients
    // because the reconstruction is just the prediction.
    int y, x;

    for (y = 0; y < width; ++y) {
      for (x = 0; x < width; ++x) {
        rec_out[x + y * out_stride] = pred_in[x + y * in_stride];
      }
    }
  }

  return has_coeffs;
}

void kvz_quant_avx2(const encoder_state_t * const state, coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type)
{
  if (state->encoder_control->scaling_list.enable){
    kvz_quant_generic(state, coef, q_coef, width, height, type, scan_idx, block_type);
  }
  else {
    kvz_quant_flat_avx2(state, coef, q_coef, width, height, type, scan_idx, block_type);
  }
}

/**
 * \brief inverse quantize transformed and quantized coefficents
 *
 */
void kvz_dequant_avx2(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  int32_t shift,add,coeff_q;
  int32_t n;
  int32_t transform_shift = 15 - encoder->bitdepth - (kvz_g_convert_to_bit[ width ] + 2);

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth-8)*6);

  shift = 20 - QUANT_SHIFT - transform_shift;

  if (encoder->scaling_list.enable)
  {
    uint32_t log2_tr_size = kvz_g_convert_to_bit[ width ] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    const int32_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
    shift += 4;

    if (shift >qp_scaled / 6) {
      add = 1 << (shift - qp_scaled/6 - 1);

      for (n = 0; n < width * height; n++) {
        coeff_q = ((q_coef[n] * dequant_coef[n]) + add ) >> (shift -  qp_scaled/6);
        coef[n] = (coeff_t)CLIP(-32768,32767,coeff_q);
      }
    } else {
      for (n = 0; n < width * height; n++) {
        // Clip to avoid possible overflow in following shift left operation
        coeff_q   = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q << (qp_scaled/6 - shift));
      }
    }
  } else {
    int32_t scale = kvz_g_inv_quant_scales[qp_scaled%6] << (qp_scaled/6);
    add = 1 << (shift-1);

    __m256i v_scale = _mm256_set1_epi32(scale);
    __m256i v_add = _mm256_set1_epi32(add);

    for (n = 0; n < width*height; n+=16) {
      __m128i temp0 = _mm_loadu_si128((__m128i*)&(q_coef[n]));
      __m128i temp1 = _mm_loadu_si128((__m128i*)&(q_coef[n + 8]));
      __m256i v_coeff_q_lo = _mm256_cvtepi16_epi32(_mm_unpacklo_epi64(temp0, temp1));
      __m256i v_coeff_q_hi = _mm256_cvtepi16_epi32(_mm_unpackhi_epi64(temp0, temp1));
      v_coeff_q_lo = _mm256_mullo_epi32(v_coeff_q_lo, v_scale);
      v_coeff_q_hi = _mm256_mullo_epi32(v_coeff_q_hi, v_scale);
      v_coeff_q_lo = _mm256_add_epi32(v_coeff_q_lo, v_add);
      v_coeff_q_hi = _mm256_add_epi32(v_coeff_q_hi, v_add);
      v_coeff_q_lo = _mm256_srai_epi32(v_coeff_q_lo, shift);
      v_coeff_q_hi = _mm256_srai_epi32(v_coeff_q_hi, shift);
      v_coeff_q_lo = _mm256_packs_epi32(v_coeff_q_lo, v_coeff_q_hi);
      _mm_storeu_si128((__m128i*)&(coef[n]), _mm256_castsi256_si128(v_coeff_q_lo) );
      _mm_storeu_si128((__m128i*)&(coef[n + 8]), _mm256_extracti128_si256(v_coeff_q_lo, 1) );
    }
  }
}

static uint32_t coeff_abs_sum_avx2(const coeff_t *coeffs, const size_t length)
{
  assert(length % 8 == 0);

  __m256i total = _mm256_abs_epi32(_mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i*) coeffs)));

  for (int i = 8; i < length; i += 8) {
    __m256i temp = _mm256_abs_epi32(_mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i*) &coeffs[i])));
    total = _mm256_add_epi32(total, temp);
  }

  __m128i result128 = _mm_add_epi32(
    _mm256_castsi256_si128(total),
    _mm256_extractf128_si256(total, 1)
  );

  uint32_t parts[4];
  _mm_storeu_si128((__m128i*) parts, result128);

  return parts[0] + parts[1] + parts[2] + parts[3];
}

#endif //COMPILE_INTEL_AVX2 && defined X86_64

int kvz_strategy_register_quant_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;

#if COMPILE_INTEL_AVX2 && defined X86_64
  success &= kvz_strategyselector_register(opaque, "quant", "avx2", 40, &kvz_quant_avx2);
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "quantize_residual", "avx2", 40, &kvz_quantize_residual_avx2);
    success &= kvz_strategyselector_register(opaque, "dequant", "avx2", 40, &kvz_dequant_avx2);
  }
  success &= kvz_strategyselector_register(opaque, "coeff_abs_sum", "avx2", 0, &coeff_abs_sum_avx2);
#endif //COMPILE_INTEL_AVX2 && defined X86_64

  return success;
}
