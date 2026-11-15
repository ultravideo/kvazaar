/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 *
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "strategies/neon/quant-neon.h"

#include "global.h" // IWYU pragma: keep

#if COMPILE_ARM_NEON

#include <stdlib.h>

#include "cabac.h"
#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "rdo.h"
#include "scalinglist.h"
#include "strategies/generic/quant-generic.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-encode.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-quant.h"
#include "strategyselector.h"
#include "tables.h"
#include "transform.h"

#define QUANT_SHIFT 14
#define SCAN_SET_SIZE 16
#define LOG2_SCAN_SET_SIZE 4

/*
 * NEON versions of the element-wise loops of the quantization and
 * dequantization. The sign-hiding search depends on the coefficients in a
 * data-dependent order, so it stays scalar and is the same code as the
 * generic version.
 *
 * The coefficient functions do not depend on the pixel type and are
 * registered for every bit depth. The value ranges do not grow with the bit
 * depth: the QP range is extended by QP_BD_OFFSET = 6 * (bitdepth - 8), but the
 * transform shift is smaller by bitdepth - 8, so q_bits = 14 + qp_scaled / 6 +
 * 15 - bitdepth - log2(size) is at most 27 at every bit depth.
 */

/* Quantization of four coefficients: returns the unsigned levels in 32-bit
 * lanes, and the 64-bit (abs * quant + add) >> q_bits values are kept for the
 * sign-hiding pass. */
static void kvz_quant_neon(const encoder_state_t * const state, const coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t * const scan = kvz_g_sig_last_scan[scan_idx][log2_block_size - 1];

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6, encoder->cfg.chroma_format == KVZ_CSP_420);
  const uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  const coeff_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  const int32_t q_bits = QUANT_SHIFT + qp_scaled / 6 + transform_shift;
  const int32_t add = ((state->frame->slicetype == KVZ_SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;
  const int32_t n_total = width * height;

  uint32_t ac_sum = 0;
  int32_t n = 0;

  /* 32-bit lanes are exact: |coef| <= 32768 and the scaling factor is an int16,
   * so |abs * q| <= 2^30, and adding add (< 2^26 as q_bits <= 27) still fits in int32.
   * vabsq_s16(-32768) wraps to -32768, which as u16 is the correct 32768. */
  const int32x4_t vadd = vdupq_n_s32(add);
  const int32x4_t vneg_qbits = vdupq_n_s32(-q_bits);
  int32x4_t vacc = vdupq_n_s32(0);
  for (; n + 8 <= n_total; n += 8) {
    const int16x8_t c = vld1q_s16(coef + n);
    const uint16x8_t a16 = vreinterpretq_u16_s16(vabsq_s16(c));
    const int16x8_t q = vld1q_s16(quant_coeff + n);
    int32x4_t ll = vmlaq_s32(vadd, vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(a16))), vmovl_s16(vget_low_s16(q)));
    int32x4_t lh = vmlaq_s32(vadd, vreinterpretq_s32_u32(vmovl_high_u16(a16)), vmovl_high_s16(q));
    ll = vshlq_s32(ll, vneg_qbits);
    lh = vshlq_s32(lh, vneg_qbits);
    vacc = vaddq_s32(vacc, vaddq_s32(ll, lh));

    // level * sign, with sign = -1 for negative coefficients (and +1 for zero).
    const int16x8_t neg = vshrq_n_s16(c, 15);
    const int32x4_t nl = vmovl_s16(vget_low_s16(neg));
    const int32x4_t nh = vmovl_high_s16(neg);
    ll = vsubq_s32(veorq_s32(ll, nl), nl);
    lh = vsubq_s32(veorq_s32(lh, nh), nh);
    vst1q_s16(q_coef + n, vqmovn_high_s32(vqmovn_s32(ll), lh));
  }
  ac_sum += (uint32_t)vaddvq_s32(vacc);
  for (; n < n_total; ++n) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t sign = (level < 0 ? -1 : 1);

    int32_t curr_quant_coeff = quant_coeff[n];
    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    ac_sum += level;

    level *= sign;
    q_coef[n] = (coeff_t)(CLIP(-32768, 32767, level));
  }

  if (!encoder->cfg.signhide_enable || ac_sum < 2) return;

  int32_t delta_u[LCU_WIDTH*LCU_WIDTH >> 2];

  n = 0;
  const int32x4_t vqbits = vdupq_n_s32(q_bits);
  const int32x4_t vneg_qbits8 = vdupq_n_s32(-q_bits8);
  for (; n + 8 <= n_total; n += 8) {
    const int16x8_t c = vld1q_s16(coef + n);
    const uint16x8_t a16 = vreinterpretq_u16_s16(vabsq_s16(c));
    const int16x8_t q = vld1q_s16(quant_coeff + n);
    const int32x4_t pl = vmulq_s32(vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(a16))), vmovl_s16(vget_low_s16(q)));
    const int32x4_t ph = vmulq_s32(vreinterpretq_s32_u32(vmovl_high_u16(a16)), vmovl_high_s16(q));
    const int32x4_t ll = vshlq_s32(vaddq_s32(pl, vadd), vneg_qbits);
    const int32x4_t lh = vshlq_s32(vaddq_s32(ph, vadd), vneg_qbits);
    // delta_u = (abs * q - (level << q_bits)) >> q_bits8; every term fits in int32.
    vst1q_s32(delta_u + n, vshlq_s32(vsubq_s32(pl, vshlq_s32(ll, vqbits)), vneg_qbits8));
    vst1q_s32(delta_u + n + 4, vshlq_s32(vsubq_s32(ph, vshlq_s32(lh, vqbits)), vneg_qbits8));
  }
  for (; n < n_total; ++n) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t curr_quant_coeff = quant_coeff[n];

    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    delta_u[n] = (int32_t)((abs_level * curr_quant_coeff - ((int64_t)level << q_bits)) >> q_bits8);
  }

  // Sign data hiding: sequential, same as the generic code.
  {
    const uint32_t * const scan_ = scan;
    int32_t n2, last_cg = -1, abssum = 0, subset, subpos;
    for (subset = (width*height - 1) >> LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg = -1;
      subpos = subset << LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n2 = SCAN_SET_SIZE - 1; n2 >= 0; n2--)  {
        if (q_coef[scan_[n2 + subpos]])  {
          last_nz_pos_in_cg = n2;
          break;
        }
      }

      // First coeff pos
      for (n2 = 0; n2 < SCAN_SET_SIZE; n2++) {
        if (q_coef[scan_[n2 + subpos]]) {
          first_nz_pos_in_cg = n2;
          break;
        }
      }

      // Sum all kvz_quant coeffs between first and last
      for (n2 = first_nz_pos_in_cg; n2 <= last_nz_pos_in_cg; n2++) {
        abssum += q_coef[scan_[n2 + subpos]];
      }

      if (last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan_[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1);
        if (signbit != (abssum & 0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff, min_pos = -1, cur_cost = 0x7fffffff;
          int16_t final_change = 0, cur_change = 0;
          for (n2 = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n2 >= 0; n2--) {
            uint32_t blkPos = scan_[n2 + subpos];
            if (q_coef[blkPos] != 0) {
              if (delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change = 1;
              }
              else if (n2 == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost = 0x7fffffff;
              }
              else {
                cur_cost = delta_u[blkPos];
                cur_change = -1;
              }
            }
            else if (n2 < first_nz_pos_in_cg && ((coef[blkPos] >= 0) ? 0 : 1) != signbit) {
              cur_cost = 0x7fffffff;
            }
            else {
              cur_cost = -delta_u[blkPos];
              cur_change = 1;
            }

            if (cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos = blkPos;
            }
          } // CG loop

          if (q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if (coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;
        } // Hide
      }
      if (last_cg == 1) last_cg = 0;
    }
  }
}

static void kvz_dequant_neon(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width,
  int32_t height, int8_t type, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  int32_t shift, add;
  const int32_t n_total = width * height;
  int32_t transform_shift = 15 - encoder->bitdepth - (kvz_g_convert_to_bit[width] + 2);

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6, encoder->cfg.chroma_format == KVZ_CSP_420);

  shift = 20 - QUANT_SHIFT - transform_shift;

  int32_t n = 0;
  if (encoder->scaling_list.enable) {
    uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    const coeff_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
    shift += 4;

    if (shift > qp_scaled / 6) {
      add = 1 << (shift - qp_scaled / 6 - 1);
      const int32x4_t vadd = vdupq_n_s32(add);
      const int32x4_t sh = vdupq_n_s32(-(shift - qp_scaled / 6));
      for (; n + 8 <= n_total; n += 8) {
        const int16x8_t q = vld1q_s16(q_coef + n);
        const int16x8_t d = vld1q_s16(dequant_coef + n);
        // (q * dq + add) >> (shift - qp/6); the int16 x int16 products are exact.
        const int32x4_t lo = vshlq_s32(vaddq_s32(vmull_s16(vget_low_s16(q), vget_low_s16(d)), vadd), sh);
        const int32x4_t hi = vshlq_s32(vaddq_s32(vmull_high_s16(q, d), vadd), sh);
        vst1q_s16(coef + n, vqmovn_high_s32(vqmovn_s32(lo), hi));
      }
      for (; n < n_total; ++n) {
        int32_t coeff_q = ((q_coef[n] * dequant_coef[n]) + add) >> (shift - qp_scaled / 6);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
      }
    } else {
      const int16x8_t ls = vdupq_n_s16((int16_t)(qp_scaled / 6 - shift));
      for (; n + 8 <= n_total; n += 8) {
        const int16x8_t q = vld1q_s16(q_coef + n);
        const int16x8_t d = vld1q_s16(dequant_coef + n);
        // Clip to 16 bits before the left shift, as the generic code does; the
        // saturating 16-bit shift is the clipped 32-bit shift.
        const int16x8_t p = vqmovn_high_s32(vqmovn_s32(vmull_s16(vget_low_s16(q), vget_low_s16(d))), vmull_high_s16(q, d));
        vst1q_s16(coef + n, vqshlq_s16(p, ls));
      }
      for (; n < n_total; ++n) {
        int32_t coeff_q = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q << (qp_scaled / 6 - shift));
      }
    }
  } else {
    int32_t scale = kvz_g_inv_quant_scales[qp_scaled % 6] << (qp_scaled / 6);
    add = 1 << (shift - 1);
    const int32x4_t vscale = vdupq_n_s32(scale);
    const int32x4_t sh = vdupq_n_s32(-shift);
    const int32x4_t vadd = vdupq_n_s32(add);
    for (; n + 8 <= n_total; n += 8) {
      // (q * scale + add) >> shift, with the products wrapping like the generic int32 code.
      const int16x8_t q = vld1q_s16(q_coef + n);
      const int32x4_t lo = vshlq_s32(vmlaq_s32(vadd, vmovl_s16(vget_low_s16(q)), vscale), sh);
      const int32x4_t hi = vshlq_s32(vmlaq_s32(vadd, vmovl_high_s16(q), vscale), sh);
      vst1q_s16(coef + n, vqmovn_high_s32(vqmovn_s32(lo), hi));
    }
    for (; n < n_total; ++n) {
      int32_t coeff_q = (q_coef[n] * scale + add) >> shift;
      coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}

static uint32_t coeff_abs_sum_neon(const coeff_t *coeffs, size_t length)
{
  uint32x4_t acc = vdupq_n_u32(0);
  size_t i = 0;
  for (; i + 8 <= length; i += 8) {
    const int16x8_t v = vld1q_s16(coeffs + i);
    acc = vpadalq_u16(acc, vreinterpretq_u16_s16(vabsq_s16(v)));
  }
  uint32_t sum = vaddvq_u32(acc);
  for (; i < length; ++i) {
    sum += abs(coeffs[i]);
  }
  return sum;
}

static double fast_coeff_cost_neon(const coeff_t *coeff, int32_t width, uint64_t weights)
{
  // Weights are four 16-bit values packed in a 64-bit word. They are looked up
  // with a byte table, once for the low byte and once for the high byte.
  uint16_t w[4];
  w[0] = (weights >>  0) & 0xffff;
  w[1] = (weights >> 16) & 0xffff;
  w[2] = (weights >> 32) & 0xffff;
  w[3] = (weights >> 48) & 0xffff;

  const uint8_t lo_bytes[8] = { (uint8_t)(w[0] & 0xff), (uint8_t)(w[1] & 0xff), (uint8_t)(w[2] & 0xff), (uint8_t)(w[3] & 0xff), 0, 0, 0, 0 };
  const uint8_t hi_bytes[8] = { (uint8_t)(w[0] >> 8), (uint8_t)(w[1] >> 8), (uint8_t)(w[2] >> 8), (uint8_t)(w[3] >> 8), 0, 0, 0, 0 };
  const uint8x8_t tbl_lo = vld1_u8(lo_bytes);
  const uint8x8_t tbl_hi = vld1_u8(hi_bytes);
  const uint8x8_t three = vdup_n_u8(3);

  const int32_t n_total = width * width;
  uint32x4_t acc = vdupq_n_u32(0);
  int32_t i = 0;
  for (; i + 8 <= n_total; i += 8) {
    const int16x8_t v = vld1q_s16(coeff + i);
    const uint8x8_t idx = vmin_u8(vqmovun_s16(vqabsq_s16(v)), three);  // vqabs saturates -32768 to 32767, as the generic abs does for the > 3 clamp
    const uint8x8_t lo = vtbl1_u8(tbl_lo, idx);
    const uint8x8_t hi = vtbl1_u8(tbl_hi, idx);
    const uint16x8_t wv = vorrq_u16(vshlq_n_u16(vmovl_u8(hi), 8), vmovl_u8(lo));
    acc = vpadalq_u16(acc, wv);
  }
  uint32_t sum = vaddvq_u32(acc);
  for (; i < n_total; ++i) {
    int16_t curr = coeff[i];
    uint32_t curr_abs = abs(curr);
    if (curr_abs > 3) {
      curr_abs = 3;
    }
    sum += w[curr_abs];
  }
  return (double)sum / 256.0;
}

/*
 * 4x4 scans of the coefficient groups (kvz_g_sig_last_scan[scan_mode][1]) as
 * byte tables for vqtbl1q_u8: scan index -> raster position in the group, and
 * raster position -> scan index. Filled when the strategies are registered.
 */
static uint8_t cg_scan[3][16];
static uint8_t cg_scan_inv[3][16];

static void init_cg_scan_tables(void)
{
  for (int mode = 0; mode < 3; ++mode) {
    for (int i = 0; i < 16; ++i) {
      const uint32_t raster = kvz_g_sig_last_scan[mode][1][i];
      cg_scan[mode][i] = (uint8_t)raster;
      cg_scan_inv[mode][raster] = (uint8_t)i;
    }
  }
}

/*
 * Finds the last coefficient in scan order that would quantize to a non-zero
 * level, as find_last_scanpos_generic. The coefficient groups are 4x4 blocks:
 * the first scan position of a group is its top-left corner and the group is
 * scanned in the 4x4 order of scan_mode, so the levels of a whole group are
 * computed from four row loads. In the group with the last level, the
 * non-zero mask is permuted to scan order to find the last level, and the
 * coefficients after it are cleared with one masked store per row.
 */
static void find_last_scanpos_neon(coeff_t* coef, coeff_t* dest_coeff, int8_t type, int32_t q_bits, const coeff_t* quant_coeff, struct kvz_sh_rates_t* sh_rates, const uint32_t cg_size,
  uint16_t* ctx_set, const uint32_t* scan, int32_t* cg_last_scanpos, int32_t* last_scanpos, uint32_t cg_num, int32_t* cg_scanpos, int32_t width, int8_t scan_mode)
{
  const int32x4_t vmax = vdupq_n_s32(MAX_INT - (1 << (q_bits - 1)));
  const int32x4_t vhalf = vdupq_n_s32(1 << (q_bits - 1));
  const int32x4_t vneg_qbits = vdupq_n_s32(-q_bits);
  assert(cg_size == 16);
  (void)cg_size;

  for (*cg_scanpos = (cg_num - 1); *cg_scanpos >= 0; (*cg_scanpos)--) {
    const int32_t scan_pos = *cg_scanpos * 16;
    const uint32_t block_pos = scan[scan_pos];
    coeff_t *dest[4];
    uint16x4_t nonzero[4];
    for (int r = 0; r < 4; ++r) {
      dest[r] = dest_coeff + block_pos + r * width;
      const int16x4_t c = vld1_s16(coef + block_pos + r * width);
      const int16x4_t q = vld1_s16(quant_coeff + block_pos + r * width);
      // abs(c) * q as c * q with the sign of c removed; |c * q| <= 2^30.
      const int32x4_t neg = vmovl_s16(vshr_n_s16(c, 15));
      int32x4_t level = vsubq_s32(veorq_s32(vmull_s16(c, q), neg), neg);
      level = vminq_s32(level, vmax);
      level = vshlq_s32(vaddq_s32(level, vhalf), vneg_qbits);
      nonzero[r] = vmovn_u32(vtstq_s32(level, level));
    }
    // 0xff for the positions with a non-zero level, in raster order.
    const uint8x16_t nz = vcombine_u8(vmovn_u16(vcombine_u16(nonzero[0], nonzero[1])),
                                      vmovn_u16(vcombine_u16(nonzero[2], nonzero[3])));

    if (vmaxvq_u8(nz) == 0) {
      const int16x4_t zero = vdup_n_s16(0);
      for (int r = 0; r < 4; ++r) {
        vst1_s16(dest[r], zero);
      }
      continue;
    }

    // Index of the last non-zero level in scan order.
    static const uint8_t iota[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const uint8x16_t nz_scan = vqtbl1q_u8(nz, vld1q_u8(cg_scan[scan_mode]));
    const uint8_t last = vmaxvq_u8(vandq_u8(nz_scan, vld1q_u8(iota)));

    // Clear the coefficients after it, as the generic version does.
    const int8x16_t after = vreinterpretq_s8_u8(vcgtq_u8(vld1q_u8(cg_scan_inv[scan_mode]), vdupq_n_u8(last)));
    const int16x8_t after01 = vmovl_s8(vget_low_s8(after));
    const int16x8_t after23 = vmovl_high_s8(after);
    const int16x8_t d01 = vbicq_s16(vcombine_s16(vld1_s16(dest[0]), vld1_s16(dest[1])), after01);
    const int16x8_t d23 = vbicq_s16(vcombine_s16(vld1_s16(dest[2]), vld1_s16(dest[3])), after23);
    vst1_s16(dest[0], vget_low_s16(d01));
    vst1_s16(dest[1], vget_high_s16(d01));
    vst1_s16(dest[2], vget_low_s16(d23));
    vst1_s16(dest[3], vget_high_s16(d23));

    const int32_t scanpos = scan_pos + last;
    *last_scanpos = scanpos;
    *ctx_set = (scanpos > 0 && type == 0) ? 2 : 0;
    *cg_last_scanpos = *cg_scanpos;
    sh_rates->sig_coeff_inc[scan[scanpos]] = 0;
    return;
  }
}

/*
 * Residual and reconstruction of a width x width block. For 8-bit pixels the
 * pixels are widened to 16 bits; for 10 and 12-bit pixels they are already
 * 16-bit. Like the generic code the residual and the sum of the residual and
 * the prediction wrap to int16 before the reconstruction is clipped to
 * [0, PIXEL_MAX].
 */
#if KVZ_BIT_DEPTH == 8
static INLINE uint8x8_t load4_u8(const uint8_t *p)
{
  uint32_t v;
  memcpy(&v, p, 4);
  return vreinterpret_u8_u32(vdup_n_u32(v));
}

static INLINE void store4_u8(uint8_t *p, uint8x8_t v)
{
  const uint32_t x = vget_lane_u32(vreinterpret_u32_u8(v), 0);
  memcpy(p, &x, 4);
}
#endif

static void get_residual_neon(const kvz_pixel *ref_in, const kvz_pixel *pred_in, int16_t *residual, int width, int in_stride)
{
  if (width == 4) {
    for (int y = 0; y < 4; ++y) {
      const kvz_pixel *a = ref_in + y * in_stride;
      const kvz_pixel *b = pred_in + y * in_stride;
#if KVZ_BIT_DEPTH == 8
      const uint16x8_t d = vsubl_u8(load4_u8(a), load4_u8(b));
      vst1_s16(residual + y * 4, vreinterpret_s16_u16(vget_low_u16(d)));
#else
      vst1_s16(residual + y * 4, vreinterpret_s16_u16(vsub_u16(vld1_u16(a), vld1_u16(b))));
#endif
    }
#if KVZ_BIT_DEPTH == 8
  } else if (width % 16 == 0) {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; x += 16) {
        const uint8x16_t a = vld1q_u8(ref_in + x + y * in_stride);
        const uint8x16_t b = vld1q_u8(pred_in + x + y * in_stride);
        vst1q_s16(residual + x + y * width, vreinterpretq_s16_u16(vsubl_u8(vget_low_u8(a), vget_low_u8(b))));
        vst1q_s16(residual + x + 8 + y * width, vreinterpretq_s16_u16(vsubl_high_u8(a, b)));
      }
    }
#endif
  } else if (width % 8 == 0) {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; x += 8) {
        const kvz_pixel *a = ref_in + x + y * in_stride;
        const kvz_pixel *b = pred_in + x + y * in_stride;
#if KVZ_BIT_DEPTH == 8
        const uint16x8_t d = vsubl_u8(vld1_u8(a), vld1_u8(b));
#else
        const uint16x8_t d = vsubq_u16(vld1q_u16(a), vld1q_u16(b));
#endif
        vst1q_s16(residual + x + y * width, vreinterpretq_s16_u16(d));
      }
    }
  } else {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; ++x) {
        residual[x + y * width] = (int16_t)(ref_in[x + y * in_stride] - pred_in[x + y * in_stride]);
      }
    }
  }
}

static void get_quantized_recon_neon(const int16_t *residual, const kvz_pixel *pred_in, int in_stride,
  kvz_pixel *rec_out, int out_stride, int width)
{
#if KVZ_BIT_DEPTH > 8
  const int16x8_t max = vdupq_n_s16(PIXEL_MAX);
  const int16x8_t zero = vdupq_n_s16(0);
#endif
  if (width == 4) {
    for (int y = 0; y < 4; ++y) {
      const int16x4_t res = vld1_s16(residual + y * 4);
      const kvz_pixel *pred = pred_in + y * in_stride;
      kvz_pixel *rec = rec_out + y * out_stride;
#if KVZ_BIT_DEPTH == 8
      const uint16x8_t sum = vaddw_u8(vreinterpretq_u16_s16(vcombine_s16(res, res)), load4_u8(pred));
      store4_u8(rec, vqmovun_s16(vreinterpretq_s16_u16(sum)));
#else
      int16x4_t sum = vadd_s16(res, vreinterpret_s16_u16(vld1_u16(pred)));
      sum = vmin_s16(vmax_s16(sum, vget_low_s16(zero)), vget_low_s16(max));
      vst1_u16(rec, vreinterpret_u16_s16(sum));
#endif
    }
#if KVZ_BIT_DEPTH == 8
  } else if (width % 16 == 0) {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; x += 16) {
        const uint8x16_t pred = vld1q_u8(pred_in + x + y * in_stride);
        const uint16x8_t lo = vaddw_u8(vreinterpretq_u16_s16(vld1q_s16(residual + x + y * width)), vget_low_u8(pred));
        const uint16x8_t hi = vaddw_high_u8(vreinterpretq_u16_s16(vld1q_s16(residual + x + 8 + y * width)), pred);
        vst1q_u8(rec_out + x + y * out_stride,
                 vqmovun_high_s16(vqmovun_s16(vreinterpretq_s16_u16(lo)), vreinterpretq_s16_u16(hi)));
      }
    }
#endif
  } else if (width % 8 == 0) {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; x += 8) {
        const int16x8_t res = vld1q_s16(residual + x + y * width);
        const kvz_pixel *pred = pred_in + x + y * in_stride;
        kvz_pixel *rec = rec_out + x + y * out_stride;
#if KVZ_BIT_DEPTH == 8
        const uint16x8_t sum = vaddw_u8(vreinterpretq_u16_s16(res), vld1_u8(pred));
        vst1_u8(rec, vqmovun_s16(vreinterpretq_s16_u16(sum)));
#else
        int16x8_t sum = vaddq_s16(res, vreinterpretq_s16_u16(vld1q_u16(pred)));
        sum = vminq_s16(vmaxq_s16(sum, zero), max);
        vst1q_u16(rec, vreinterpretq_u16_s16(sum));
#endif
      }
    }
  } else {
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; ++x) {
        int16_t val = residual[x + y * width] + pred_in[x + y * in_stride];
        rec_out[x + y * out_stride] = (kvz_pixel)CLIP(0, PIXEL_MAX, val);
      }
    }
  }
}

static bool has_nonzero_coeffs_neon(const coeff_t *coeff, int n)
{
  int i = 0;
  // Blocks with coefficients usually have them at the start (in raster
  // order), so check 32 coefficients at a time and stop at the first non-zero.
  for (; i + 32 <= n; i += 32) {
    const int16x8_t a = vorrq_s16(vld1q_s16(coeff + i), vld1q_s16(coeff + i + 8));
    const int16x8_t b = vorrq_s16(vld1q_s16(coeff + i + 16), vld1q_s16(coeff + i + 24));
    if (vmaxvq_u16(vreinterpretq_u16_s16(vorrq_s16(a, b))) != 0) return true;
  }
  for (; i + 8 <= n; i += 8) {
    if (vmaxvq_u16(vreinterpretq_u16_s16(vld1q_s16(coeff + i))) != 0) return true;
  }
  for (; i < n; ++i) {
    if (coeff[i] != 0) return true;
  }
  return false;
}

/*
 * Cross-component prediction (4:4:4 only). The sums of squares are computed
 * exactly in 64 bits: vmull_s16 gives exact 32-bit products of the residuals,
 * which vpadalq_s32 accumulates pairwise into 64-bit lanes.
 */
static void ccp_apply_neon(const int16_t *luma_residual, int16_t *chroma_residual, int32_t tr_width,
  int32_t luma_stride, int32_t chroma_stride, int32_t reversed_alpha, bool saturate)
{
  for (int32_t j = 0; j < tr_width; j++) {
    const int16_t *luma = luma_residual + j * luma_stride;
    int16_t *chroma = chroma_residual + j * chroma_stride;
    int32_t i = 0;
    for (; i + 4 <= tr_width; i += 4) {
      // (reversed_alpha * luma) >> 3, added to (recon) or subtracted from the chroma residual.
      const int32x4_t pred = vshrq_n_s32(vmulq_n_s32(vmovl_s16(vld1_s16(luma + i)), reversed_alpha), 3);
      const int32x4_t c = vmovl_s16(vld1_s16(chroma + i));
      if (saturate) {
        vst1_s16(chroma + i, vqmovn_s32(vaddq_s32(c, pred)));
      } else {
        vst1_s16(chroma + i, vmovn_s32(vsubq_s32(c, pred)));
      }
    }
    for (; i < tr_width; i++) {
      const int32_t pred_val = luma[i];
      if (saturate) {
        chroma[i] = (int16_t)CLIP(-32768, 32767, chroma[i] + ((reversed_alpha * pred_val) >> 3));
      } else {
        chroma[i] -= ((reversed_alpha * pred_val) >> 3);
      }
    }
  }
}

static int8_t calc_cross_component_prediction_neon(cu_info_t* const cur_cu, color_t color,
  int16_t* luma_residual, int16_t* chroma_residual, int32_t tr_width, int32_t luma_stride, int32_t chroma_stride)
{
  // Signed, because a negative correlation gives a negative alpha.
  int64_t ss_xy = 0;
  int64_t ss_xx = 0;
  int64x2_t v_xx = vdupq_n_s64(0);
  int64x2_t v_xy = vdupq_n_s64(0);
  for (int32_t j = 0; j < tr_width; j++) {
    const int16_t *luma = luma_residual + j * luma_stride;
    const int16_t *chroma = chroma_residual + j * chroma_stride;
    int32_t i = 0;
    for (; i + 4 <= tr_width; i += 4) {
      const int16x4_t l = vld1_s16(luma + i);
      const int16x4_t c = vld1_s16(chroma + i);
      v_xx = vpadalq_s32(v_xx, vmull_s16(l, l));
      v_xy = vpadalq_s32(v_xy, vmull_s16(l, c));
    }
    for (; i < tr_width; i++) {
      ss_xy += luma[i] * chroma[i];
      ss_xx += luma[i] * luma[i];
    }
  }
  ss_xx += vaddvq_s64(v_xx);
  ss_xy += vaddvq_s64(v_xy);

  int8_t alpha = 0;
  int8_t sign = 0;
  if (ss_xx != 0) {
    double falpha = (double)ss_xy / (double)ss_xx;
    // Clip before the conversion, which is undefined for values out of the int range.
    alpha = (int8_t)CLIP(-16.0, 16.0, falpha * 16);
    sign = alpha < 0 ? 1 : 0;
    alpha = alpha < 0 ? -alpha : alpha;
    const uint8_t combined_alpha_quant[17] = { 0,1,1,2,2,2,3,3,3,3,3,3,4,4,4,4,4, };

    alpha = combined_alpha_quant[alpha];

    if (alpha != 0) {
      int32_t reversed_alpha = sign ? -(1 << (alpha - 1)) : (1 << (alpha - 1));
      ccp_apply_neon(luma_residual, chroma_residual, tr_width, luma_stride, chroma_stride, reversed_alpha, false);
    }
  }
  if (color == COLOR_V) {
    cur_cu->alpha_v = alpha;
    cur_cu->alpha_v_s = sign;
  }
  else {
    cur_cu->alpha_u = alpha;
    cur_cu->alpha_u_s = sign;
  }

  return alpha;
}

static int8_t recon_cross_component_prediction_neon(cu_info_t* const cur_cu, color_t color,
  int16_t* luma_residual, int16_t* chroma_residual, int32_t tr_width, int32_t luma_stride, int32_t chroma_stride)
{
  int8_t alpha = (color == COLOR_V) ? cur_cu->alpha_v : cur_cu->alpha_u;
  uint8_t sign = (color == COLOR_V) ? cur_cu->alpha_v_s : cur_cu->alpha_u_s;
  if (alpha != 0) {
    int32_t reversed_alpha = sign ? -(1 << (alpha - 1)) : (1 << (alpha - 1));
    ccp_apply_neon(luma_residual, chroma_residual, tr_width, luma_stride, chroma_stride, reversed_alpha, true);
  }
  return alpha;
}

static bool cross_component_prediction_rdo_neon(encoder_state_t* const state,
  cu_info_t* const cur_cu, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const kvz_pixel* const pred_in,
  int16_t* const luma_residual_cross_comp[2], int16_t* const residual,
  coeff_t* const coeff_out, const int width, const kvz_pixel* const ref_in)
{
  if (color == COLOR_Y || !cbf_is_set(cur_cu->cbf, cur_cu->depth, COLOR_Y)) {
    return false;
  }

  ALIGNED(64) int16_t residual_backup[TR_MAX_WIDTH * TR_MAX_WIDTH];
  memcpy(residual_backup, residual, width * width * sizeof(int16_t));
  int8_t calculated_alpha = calc_cross_component_prediction_neon(cur_cu, color,
    luma_residual_cross_comp[0], residual, width, state->tile->frame->width, width);

  if (!calculated_alpha) {
    return false;
  }

  // Check if the cross component prediction is worth using via simple RDO
  ALIGNED(64) coeff_t test_coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) coeff_t coeff_out_temp[TR_MAX_WIDTH * TR_MAX_WIDTH];

  int16_t* test_residual[2] = { residual_backup, residual };
  double cost[2] = { 0, 0 };
  cabac_data_t cabac_copy;
  memcpy(&cabac_copy, &state->cabac, sizeof(cabac_copy));

  cabac_copy.only_count = 1;

  for (int i = 0; i < 2; i++) {
    // Transform residual. (residual -> coeff)
    if (use_trskip) {
      kvz_transformskip(state->encoder_control, test_residual[i], test_coeff, width);
    }
    else {
      kvz_transform2d(state->encoder_control, test_residual[i], test_coeff, width, color, cur_cu->type);
    }

    // Quantize coeffs. (coeff -> coeff_out)
    if (state->encoder_control->cfg.rdoq_enable && (width > 4 || !state->encoder_control->cfg.rdoq_skip)) {
      kvz_rdoq(state, test_coeff, coeff_out, width, width, 2, scan_order, cur_cu->type, 0);
    }
    else {
      kvz_quant(state, test_coeff, coeff_out, width, width, 2, scan_order, cur_cu->type);
    }
    if (i == 0) {
      memcpy(coeff_out_temp, coeff_out, width * width * sizeof(coeff_t));
    }

    double bits = 0;

    // Execute the coding function.
    // It is safe to drop the const modifier since state won't be modified
    // when cabac.only_count is set.
    if (has_nonzero_coeffs_neon(coeff_out, width * width)) {
      kvz_encode_coeff_nxn((encoder_state_t*)state, &cabac_copy, coeff_out, width, 2, scan_order, use_trskip, &bits);
      // Get quantized residual. (coeff_out -> coeff -> residual)
      kvz_dequant(state, coeff_out, test_coeff, width, width, (color == COLOR_U ? 2 : 3), cur_cu->type);
      ALIGNED(64) int16_t recon_residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
      ALIGNED(64) kvz_pixel rec_out_temp[TR_MAX_WIDTH * TR_MAX_WIDTH];

      if (use_trskip) {
        kvz_itransformskip(state->encoder_control, recon_residual, test_coeff, width);
      }
      else {
        kvz_itransform2d(state->encoder_control, recon_residual, test_coeff, width, color, cur_cu->type);
      }
      if (i == 1) {
        recon_cross_component_prediction_neon(cur_cu, color, luma_residual_cross_comp[1], recon_residual,
          width, state->tile->frame->width, width);
      }

      // Get quantized reconstruction. (residual + pred_in -> rec_out)
      get_quantized_recon_neon(recon_residual, pred_in, in_stride, rec_out_temp, width, width);

      // Calculate SATD against original block
      cost[i] = (double)(kvz_satd_any_size(width, width, ref_in, width, rec_out_temp, width));
    }
    cost[i] += bits * state->lambda;
  }
  {
    double bits = 0;
    double bits_noalpha = 0;
    int8_t alpha = calculated_alpha;
    int8_t alpha_sign = alpha < 0 ? 1 : 0;
    alpha = alpha < 0 ? -alpha : alpha;
    // Add cross-component prediction signal cost
    cabac_ctx_t* ctx = &(cabac_copy.ctx.cross_component_prediction[color == COLOR_V ? 5 : 0]);
    CABAC_FBITS_UPDATE(&cabac_copy, ctx, 1, bits, "cross_component_prediction_flag");
    if (alpha != 0) {
      alpha--;
      CABAC_FBITS_UPDATE(&cabac_copy, &ctx[1], (alpha > 0) ? 1 : 0, bits, "cross_component_prediction_alpha");
      if (alpha > 0) {
        kvz_cabac_write_unary_max_symbol(&cabac_copy, &ctx[2], alpha - 1, 1, 2, &bits);
      }
      CABAC_FBITS_UPDATE(&cabac_copy, &ctx[4], alpha_sign, bits, "cross_component_prediction_sign");
    }

    cost[1] += bits * state->lambda;

    // Add cross-component prediction signal cost for no alpha
    CABAC_FBITS_UPDATE(&cabac_copy, ctx, 0, bits_noalpha, "cross_component_prediction_flag");

    cost[0] += bits_noalpha * state->lambda;
  }

  // If the cost is not reduced, revert the cross-component prediction
  if (cost[1] >= cost[0]) {
    memcpy(residual, residual_backup, width * width * sizeof(int16_t));
    memcpy(coeff_out, coeff_out_temp, width * width * sizeof(coeff_t));
    if (color == COLOR_U) {
      cur_cu->alpha_u = 0;
    }
    else {
      cur_cu->alpha_v = 0;
    }
  }
  return true;
}

/*
 * Quantize residual and get both the reconstruction and coeffs, as
 * kvz_quantize_residual_generic. The residual, the check for non-zero
 * coefficients, the reconstruction and the cross-component prediction use
 * NEON; the transforms and the quantizers are called through the strategy
 * pointers. Instantiated with ccp_enabled true and false, so that the
 * cross-component prediction code is removed when the tool is not used.
 */
static FORCE_INLINE int quantize_residual_neon_impl(encoder_state_t *const state,
  cu_info_t *const cur_cu, const int width, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const kvz_pixel *const ref_in, const kvz_pixel *const pred_in,
  kvz_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip, int16_t* luma_residual_cross_comp[2],
  const bool ccp_enabled)
{
  // Temporary arrays to pass data to and from kvz_quant and transform functions.
  ALIGNED(64) int16_t residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) coeff_t coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

  assert(width <= TR_MAX_WIDTH);
  assert(width >= 2);

  const bool allow_cross_component_prediction = ccp_enabled && (cur_cu->tr_depth == cur_cu->depth);

  // When the transform RDO does not run for this block (tr_depth != depth)
  // the alpha fields may hold stale values from a candidate that was not
  // selected. Clear them so the encoder does not signal a cross-component
  // prediction that the reconstruction never applied.
  if (ccp_enabled && !allow_cross_component_prediction) {
    cur_cu->alpha_u = 0;
    cur_cu->alpha_v = 0;
    cur_cu->alpha_u_s = 0;
    cur_cu->alpha_v_s = 0;
  }

  // Get residual. (ref_in - pred_in -> residual)
  get_residual_neon(ref_in, pred_in, residual, width, in_stride);

  if (allow_cross_component_prediction && color == COLOR_Y) {
    // Store luma residual for cross-component prediction of chroma to the frame level buffer
    for (int yy = 0; yy < width; yy++) {
      memcpy(&luma_residual_cross_comp[0][yy * state->tile->frame->width], &residual[yy * width], width * sizeof(int16_t));
    }
  }

  if (!allow_cross_component_prediction ||
    !cross_component_prediction_rdo_neon(state, cur_cu, color, scan_order, use_trskip,
    in_stride, pred_in, luma_residual_cross_comp, residual, coeff_out, width, ref_in)) {
    // No cross-component prediction RDO done so we continue as usual, otherwise tr-quant already done
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
    }
    else {
      kvz_quant(state, coeff, coeff_out, width, width, (color == COLOR_Y ? 0 : 2),
        scan_order, cur_cu->type);
    }
  }

  const int has_coeffs = has_nonzero_coeffs_neon(coeff_out, width * width);

  // Do the inverse quantization and transformation and the reconstruction to
  // rec_out.
  if (has_coeffs && !early_skip) {
    // Get quantized residual. (coeff_out -> coeff -> residual)
    kvz_dequant(state, coeff_out, coeff, width, width, (color == COLOR_Y ? 0 : (color == COLOR_U ? 2 : 3)), cur_cu->type);
    if (use_trskip) {
      kvz_itransformskip(state->encoder_control, residual, coeff, width);
    }
    else {
      kvz_itransform2d(state->encoder_control, residual, coeff, width, color, cur_cu->type);
    }
    if (allow_cross_component_prediction) {
      if (color == COLOR_Y) {
        // Store luma residual for cross-component prediction of chroma to the frame level buffer
        for (int yy = 0; yy < width; yy++) {
          memcpy(&luma_residual_cross_comp[1][yy * state->tile->frame->width], &residual[yy * width], width * sizeof(int16_t));
        }
      }
      else if (cbf_is_set(cur_cu->cbf, cur_cu->depth, COLOR_Y)) {
        recon_cross_component_prediction_neon(cur_cu, color, luma_residual_cross_comp[1], residual, width, state->tile->frame->width, width);
      }
    }

    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    get_quantized_recon_neon(residual, pred_in, in_stride, rec_out, out_stride, width);
  }
  else {
    bool recon_done = false;
    if (allow_cross_component_prediction && cbf_is_set(cur_cu->cbf, cur_cu->depth, COLOR_Y)) {
      memset(residual, 0, width * width * sizeof(int16_t));
      if (recon_cross_component_prediction_neon(cur_cu, color, luma_residual_cross_comp[1], residual, width, state->tile->frame->width, width) != 0) {
        // Get quantized reconstruction. (residual + pred_in -> rec_out)
        get_quantized_recon_neon(residual, pred_in, in_stride, rec_out, out_stride, width);
        recon_done = true;
      }
    }
    if (!recon_done && rec_out != pred_in) {
      // With no coeffs and rec_out == pred_int we skip copying the coefficients
      // because the reconstruction is just the prediction.
      for (int y = 0; y < width; ++y) {
        memcpy(&rec_out[y * out_stride], &pred_in[y * in_stride], width * sizeof(kvz_pixel));
      }
    }
  }

  return has_coeffs;
}

static int quantize_residual_neon(encoder_state_t *const state,
  cu_info_t *const cur_cu, const int width, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const kvz_pixel *const ref_in, const kvz_pixel *const pred_in,
  kvz_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip, int16_t* luma_residual_cross_comp[2])
{
  if (state->encoder_control->cfg.enable_cross_component_prediction) {
    return quantize_residual_neon_impl(state, cur_cu, width, color, scan_order, use_trskip,
      in_stride, out_stride, ref_in, pred_in, rec_out, coeff_out,
      early_skip, luma_residual_cross_comp, true);
  }
  return quantize_residual_neon_impl(state, cur_cu, width, color, scan_order, use_trskip,
    in_stride, out_stride, ref_in, pred_in, rec_out, coeff_out,
    early_skip, luma_residual_cross_comp, false);
}

#endif // COMPILE_ARM_NEON

int kvz_strategy_register_quant_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_ARM_NEON
  // The residual and reconstruction are compiled for the pixel type of
  // KVZ_BIT_DEPTH, the bitdepth the encoder always uses.
  if (bitdepth == KVZ_BIT_DEPTH) {
    success &= kvz_strategyselector_register(opaque, "quantize_residual", "neon", 40, &quantize_residual_neon);
  }
  init_cg_scan_tables();
  // The coefficient functions do not depend on the pixel type.
  success &= kvz_strategyselector_register(opaque, "quant", "neon", 40, &kvz_quant_neon);
  success &= kvz_strategyselector_register(opaque, "dequant", "neon", 40, &kvz_dequant_neon);
  success &= kvz_strategyselector_register(opaque, "coeff_abs_sum", "neon", 40, &coeff_abs_sum_neon);
  success &= kvz_strategyselector_register(opaque, "fast_coeff_cost", "neon", 40, &fast_coeff_cost_neon);
  success &= kvz_strategyselector_register(opaque, "find_last_scanpos", "neon", 40, &find_last_scanpos_neon);
#else
  (void)opaque;
  (void)bitdepth;
#endif
  return success;
}
