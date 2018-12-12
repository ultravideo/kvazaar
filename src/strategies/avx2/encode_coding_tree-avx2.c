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

#include "strategyselector.h"

#include "cabac.h"
#include "context.h"
#include "encode_coding_tree-avx2.h"
#include "kvz_math.h"
#include <immintrin.h>

static INLINE __m256i scanord_read_vector(const int16_t *coeff, const uint32_t *scan, int8_t scan_mode, int32_t subpos, int32_t width)
{
  // For vectorized reordering of coef and q_coef
  const __m128i low128_shuffle_masks[3] = {
    _mm_setr_epi8(10,11,  4, 5, 12,13,  0, 1,  6, 7, 14,15,  8, 9,  2, 3),
    _mm_setr_epi8( 0, 1,  2, 3,  4, 5,  6, 7,  8, 9, 10,11, 12,13, 14,15),
    _mm_setr_epi8( 4, 5,  6, 7,  0, 1,  2, 3, 12,13, 14,15,  8, 9, 10,11),
  };

  const __m128i blend_masks[3] = {
    _mm_setr_epi16( 0,  0,  0, -1,  0,  0, -1, -1),
    _mm_setr_epi16( 0,  0,  0,  0,  0,  0,  0,  0),
    _mm_setr_epi16( 0,  0, -1, -1,  0,  0, -1, -1),
  };

  const __m128i invec_rearr_masks_upper[3] = {
    _mm_setr_epi8( 0, 1,  8, 9,  2, 3,  6, 7, 10,11,  4, 5, 12,13, 14,15),
    _mm_setr_epi8( 0, 1,  2, 3,  4, 5,  6, 7,  8, 9, 10,11, 12,13, 14,15),
    _mm_setr_epi8( 0, 1,  8, 9,  4, 5, 12,13,  2, 3, 10,11,  6, 7, 14,15),
  };

  const __m128i invec_rearr_masks_lower[3] = {
    _mm_setr_epi8(12,13,  6, 7,  0, 1,  2, 3, 14,15,  4, 5,  8, 9, 10,11),
    _mm_setr_epi8( 0, 1,  2, 3,  4, 5,  6, 7,  8, 9, 10,11, 12,13, 14,15),
    _mm_setr_epi8( 4, 5, 12,13,  0, 1,  8, 9,  6, 7, 14,15,  2, 3, 10,11),
  };

  const size_t row_offsets[4] = {
    scan[subpos] + width * 0,
    scan[subpos] + width * 1,
    scan[subpos] + width * 2,
    scan[subpos] + width * 3,
  };

  // NOTE: Upper means "higher in pixel order inside block", which implies
  // lower addresses (note the difference: HIGH and LOW vs UPPER and LOWER),
  // so upper 128b vector actually becomes the lower part of a 256-bit coeff
  // vector and lower vector the higher part!
  __m128d   coeffs_d_upper = _mm_castsi128_pd(_mm_set1_epi8(0));
  __m128d   coeffs_d_lower = _mm_castsi128_pd(_mm_set1_epi8(0));

  __m128i   coeffs_upper;
  __m128i   coeffs_lower;

  __m128i   coeffs_rearr1_upper;
  __m128i   coeffs_rearr1_lower;

  __m128i   coeffs_rearr2_upper;
  __m128i   coeffs_rearr2_lower;

  coeffs_d_upper = _mm_loadl_pd(coeffs_d_upper, (double *)(coeff + row_offsets[0]));
  coeffs_d_upper = _mm_loadh_pd(coeffs_d_upper, (double *)(coeff + row_offsets[1]));

  coeffs_d_lower = _mm_loadl_pd(coeffs_d_lower, (double *)(coeff + row_offsets[2]));
  coeffs_d_lower = _mm_loadh_pd(coeffs_d_lower, (double *)(coeff + row_offsets[3]));

  coeffs_upper   = _mm_castpd_si128(coeffs_d_upper);
  coeffs_lower   = _mm_castpd_si128(coeffs_d_lower);

  coeffs_lower   = _mm_shuffle_epi8(coeffs_lower, low128_shuffle_masks[scan_mode]);

  coeffs_rearr1_upper = _mm_blendv_epi8(coeffs_upper, coeffs_lower, blend_masks[scan_mode]);
  coeffs_rearr1_lower = _mm_blendv_epi8(coeffs_lower, coeffs_upper, blend_masks[scan_mode]);

  coeffs_rearr2_upper = _mm_shuffle_epi8(coeffs_rearr1_upper, invec_rearr_masks_upper[scan_mode]);
  coeffs_rearr2_lower = _mm_shuffle_epi8(coeffs_rearr1_lower, invec_rearr_masks_lower[scan_mode]);

  // Why, oh why, is there no _mm256_setr_m128i intrinsic in the header that
  // would do the exact same operation in the exact same way? :(
  return _mm256_insertf128_si256(_mm256_castsi128_si256(coeffs_rearr2_upper),
                                 coeffs_rearr2_lower,
                                 1);
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

/**
 * \brief Encode (X,Y) position of the last significant coefficient
 *
 * \param lastpos_x   X component of last coefficient
 * \param lastpos_y   Y component of last coefficient
 * \param width       Block width
 * \param height      Block height
 * \param type        plane type / luminance or chrominance
 * \param scan        scan type (diag, hor, ver)
 *
 * This method encodes the X and Y component within a block of the last
 * significant coefficient.
 */
static void encode_last_significant_xy(cabac_data_t * const cabac,
                                       uint8_t lastpos_x, uint8_t lastpos_y,
                                       uint8_t width, uint8_t height,
                                       uint8_t type, uint8_t scan)
{
  const int index = kvz_math_floor_log2(width) - 2;
  uint8_t ctx_offset = type ? 0 : (index * 3 + (index + 1) / 4);
  uint8_t shift = type ? index : (index + 3) / 4;

  cabac_ctx_t *base_ctx_x = (type ? cabac->ctx.cu_ctx_last_x_chroma : cabac->ctx.cu_ctx_last_x_luma);
  cabac_ctx_t *base_ctx_y = (type ? cabac->ctx.cu_ctx_last_y_chroma : cabac->ctx.cu_ctx_last_y_luma);

  if (scan == SCAN_VER) {
    SWAP(lastpos_x, lastpos_y, uint8_t);
  }

  const int group_idx_x = g_group_idx[lastpos_x];
  const int group_idx_y = g_group_idx[lastpos_y];

  // x prefix
  for (int last_x = 0; last_x < group_idx_x; last_x++) {
    cabac->cur_ctx = &base_ctx_x[ctx_offset + (last_x >> shift)];
    CABAC_BIN(cabac, 1, "last_sig_coeff_x_prefix");
  }
  if (group_idx_x < g_group_idx[width - 1]) {
    cabac->cur_ctx = &base_ctx_x[ctx_offset + (group_idx_x >> shift)];
    CABAC_BIN(cabac, 0, "last_sig_coeff_x_prefix");
  }

  // y prefix
  for (int last_y = 0; last_y < group_idx_y; last_y++) {
    cabac->cur_ctx = &base_ctx_y[ctx_offset + (last_y >> shift)];
    CABAC_BIN(cabac, 1, "last_sig_coeff_y_prefix");
  }
  if (group_idx_y < g_group_idx[height - 1]) {
    cabac->cur_ctx = &base_ctx_y[ctx_offset + (group_idx_y >> shift)];
    CABAC_BIN(cabac, 0, "last_sig_coeff_y_prefix");
  }

  // last_sig_coeff_x_suffix
  if (group_idx_x > 3) {
    const int suffix = lastpos_x - g_min_in_group[group_idx_x];
    const int bits = (group_idx_x - 2) / 2;
    CABAC_BINS_EP(cabac, suffix, bits, "last_sig_coeff_x_suffix");
  }

  // last_sig_coeff_y_suffix
  if (group_idx_y > 3) {
    const int suffix = lastpos_y - g_min_in_group[group_idx_y];
    const int bits = (group_idx_y - 2) / 2;
    CABAC_BINS_EP(cabac, suffix, bits, "last_sig_coeff_y_suffix");
  }
}

void kvz_encode_coeff_nxn_avx2(encoder_state_t * const state,
                               cabac_data_t * const cabac,
                               const coeff_t *coeff,
                               uint8_t width,
                               uint8_t type,
                               int8_t scan_mode,
                               int8_t tr_skip)
{
  const encoder_control_t * const encoder = state->encoder_control;
  int c1 = 1;
  uint8_t last_coeff_x = 0;
  uint8_t last_coeff_y = 0;
  int32_t i;
  uint32_t sig_coeffgroup_flag[8 * 8] = { 0 };

  int8_t be_valid = encoder->cfg.signhide_enable;
  int32_t scan_pos_sig;
  uint32_t go_rice_param = 0;
  uint32_t ctx_sig;

  // CONSTANTS
  const uint32_t num_blk_side    = width >> TR_MIN_LOG2_SIZE;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t *scan           =
    kvz_g_sig_last_scan[scan_mode][log2_block_size - 1];
  const uint32_t *scan_cg = g_sig_last_scan_cg[log2_block_size - 2][scan_mode];
  const uint32_t num_blocks = num_blk_side * num_blk_side;

  const __m256i zero = _mm256_set1_epi8(0);

  // Init base contexts according to block type
  cabac_ctx_t *base_coeff_group_ctx = &(cabac->ctx.cu_sig_coeff_group_model[type]);
  cabac_ctx_t *baseCtx           = (type == 0) ? &(cabac->ctx.cu_sig_model_luma[0]) :
                                 &(cabac->ctx.cu_sig_model_chroma[0]);

  // Scan all coeff groups to find out which of them have coeffs.
  // Populate sig_coeffgroup_flag with that info.

  // NOTE: Modified the functionality a bit, sig_coeffgroup_flag used to be
  // 1 if true and 0 if false, now it's "undefined but nonzero" if true and
  // 0 if false (not actually undefined, it's a bitmask representing the
  // significant coefficients' position in the group which in itself could
  // be useful information)
  int32_t scan_cg_last = -1;

  for (int32_t i = 0; i < num_blocks; i++) {
    const uint32_t cg_id = scan_cg[i];
    const uint32_t n_xbits = log2_block_size - 2; // How many lowest bits of scan_cg represent X coord
    const uint32_t cg_x = cg_id & ((1 << n_xbits) - 1);
    const uint32_t cg_y = cg_id >> n_xbits;

    const uint32_t cg_pos = cg_y * width * 4 + cg_x * 4;
    const uint32_t cg_pos_y = (cg_pos >> log2_block_size) >> TR_MIN_LOG2_SIZE;
    const uint32_t cg_pos_x = (cg_pos & (width - 1)) >> TR_MIN_LOG2_SIZE;
    const uint32_t addr = cg_pos_x + cg_pos_y * num_blk_side;

    __m128d coeffs_d_upper;
    __m128d coeffs_d_lower;
    __m128i coeffs_upper;
    __m128i coeffs_lower;
    __m256i cur_coeffs;

    coeffs_d_upper = _mm_loadl_pd(coeffs_d_upper, (double *)(coeff + cg_pos + 0 * width));
    coeffs_d_upper = _mm_loadh_pd(coeffs_d_upper, (double *)(coeff + cg_pos + 1 * width));
    coeffs_d_lower = _mm_loadl_pd(coeffs_d_lower, (double *)(coeff + cg_pos + 2 * width));
    coeffs_d_lower = _mm_loadh_pd(coeffs_d_lower, (double *)(coeff + cg_pos + 3 * width));

    coeffs_upper = _mm_castpd_si128(coeffs_d_upper);
    coeffs_lower = _mm_castpd_si128(coeffs_d_lower);

    cur_coeffs = _mm256_insertf128_si256(_mm256_castsi128_si256(coeffs_upper),
                                         coeffs_lower,
                                         1);

    __m256i coeffs_zero = _mm256_cmpeq_epi16(cur_coeffs, zero);

    uint32_t nz_coeffs_2b = ~((uint32_t)_mm256_movemask_epi8(coeffs_zero));
    sig_coeffgroup_flag[addr] = nz_coeffs_2b;

    if (nz_coeffs_2b)
      scan_cg_last = i;
  }
  // Rest of the code assumes at least one non-zero coeff.
  assert(scan_cg_last >= 0);

  ALIGNED(64) int16_t coeff_reord[LCU_WIDTH * LCU_WIDTH];
  for (int32_t i = scan_cg_last; i >= 0; i--) {
    int32_t subpos = i * 16;
    __m256i coeffs_r = scanord_read_vector(coeff, scan, scan_mode, subpos, width);
    _mm256_store_si256((__m256i *)(coeff_reord + subpos), coeffs_r);
  }

  // Find the last coeff by going backwards in scan order.
  uint32_t scan_pos_last;
  uint32_t baseaddr = scan_cg_last * 16;
  __m256i cur_coeffs = _mm256_loadu_si256((__m256i *)(coeff_reord + baseaddr));
  __m256i cur_coeffs_zeros = _mm256_cmpeq_epi16(cur_coeffs, zero);
  uint32_t nz_bytes = ~(_mm256_movemask_epi8(cur_coeffs_zeros));
  scan_pos_last = baseaddr + ((31 - _lzcnt_u32(nz_bytes)) >> 1);

  int pos_last = scan[scan_pos_last];

  // transform skip flag
  if(width == 4 && encoder->cfg.trskip_enable) {
    cabac->cur_ctx = (type == 0) ? &(cabac->ctx.transform_skip_model_luma) : &(cabac->ctx.transform_skip_model_chroma);
    CABAC_BIN(cabac, tr_skip, "transform_skip_flag");
  }

  last_coeff_x = pos_last & (width - 1);
  last_coeff_y = (uint8_t)(pos_last >> log2_block_size);

  // Code last_coeff_x and last_coeff_y
  encode_last_significant_xy(cabac,
                             last_coeff_x,
                             last_coeff_y,
                             width,
                             width,
                             type,
                             scan_mode);

  scan_pos_sig = scan_pos_last;

  int32_t abs_coeff[16];
  abs_coeff[0] = abs(coeff[pos_last]);
  uint32_t coeff_signs  = (coeff[pos_last] < 0);
  int32_t num_non_zero = 1;
  int32_t last_nz_pos_in_cg  = scan_pos_sig;
  int32_t first_nz_pos_in_cg = scan_pos_sig;
  scan_pos_sig--;

  // significant_coeff_flag
  for (i = scan_cg_last; i >= 0; i--) {
    int32_t sub_pos        = i << 4; // LOG2_SCAN_SET_SIZE;
    int32_t cg_blk_pos     = scan_cg[i];
    int32_t cg_pos_y       = cg_blk_pos / num_blk_side;
    int32_t cg_pos_x       = cg_blk_pos - (cg_pos_y * num_blk_side);

    go_rice_param = 0;

    if (i == scan_cg_last || i == 0) {
      sig_coeffgroup_flag[cg_blk_pos] = 1;
    } else {
      uint32_t sig_coeff_group   = (sig_coeffgroup_flag[cg_blk_pos] != 0);
      uint32_t ctx_sig  = kvz_context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                      cg_pos_y, width);
      cabac->cur_ctx = &base_coeff_group_ctx[ctx_sig];
      CABAC_BIN(cabac, sig_coeff_group, "coded_sub_block_flag");
    }

    if (sig_coeffgroup_flag[cg_blk_pos]) {
      int32_t pattern_sig_ctx = kvz_context_calc_pattern_sig_ctx(sig_coeffgroup_flag,
                                                             cg_pos_x, cg_pos_y, width);

      ALIGNED(64) int16_t abs_coeff_buf[16];
      ALIGNED(32) int16_t pos_ys_buf[16];
      ALIGNED(32) int16_t pos_xs_buf[16];

      const __m256i coeff_pos_zero = _mm256_castsi128_si256(_mm_cvtsi32_si128(0xffff));
      const __m128i log2_block_size_128 = _mm_cvtsi32_si128(log2_block_size);

      __m256i coeffs = _mm256_load_si256((__m256i *)(coeff_reord + sub_pos));
      __m256i sigs_inv = _mm256_cmpeq_epi16(coeffs, zero);
      __m256i is = _mm256_set1_epi16(i);
      __m256i is_zero = _mm256_cmpeq_epi16(is, zero);
      __m256i coeffs_subzero = _mm256_cmpgt_epi16(zero, coeffs);

      __m256i masked_coeffs = _mm256_andnot_si256(sigs_inv, coeffs);
      __m256i abs_coeffs = _mm256_abs_epi16(masked_coeffs);

      // TODO: obtain 16-bit block positions, maybe? :P
      __m256i blk_poses_hi = _mm256_loadu_si256((__m256i *)(scan + sub_pos + 8));
      __m256i blk_poses_lo = _mm256_loadu_si256((__m256i *)(scan + sub_pos + 0));
      __m256i blk_poses_tmp = _mm256_packs_epi32(blk_poses_lo, blk_poses_hi);
      __m256i blk_poses = _mm256_permute4x64_epi64(blk_poses_tmp, _MM_SHUFFLE(3, 1, 2, 0));

      __m256i pos_ys = _mm256_srl_epi16(blk_poses, log2_block_size_128);
      __m256i pos_xs = _mm256_sub_epi16(blk_poses, _mm256_sll_epi16(pos_ys, log2_block_size_128));

      _mm256_store_si256((__m256i *)pos_ys_buf, pos_ys);
      _mm256_store_si256((__m256i *)pos_xs_buf, pos_xs);

      __m256i encode_sig_coeff_flags_inv = _mm256_andnot_si256(is_zero, coeff_pos_zero);

      get_first_last_nz_int16(masked_coeffs, &first_nz_pos_in_cg, &last_nz_pos_in_cg);
      _mm256_store_si256((__m256i *)abs_coeff_buf, abs_coeffs);

      uint32_t esc_flags = ~(_mm256_movemask_epi8(encode_sig_coeff_flags_inv));
      uint32_t sigs = ~(_mm256_movemask_epi8(sigs_inv));
      uint32_t coeff_sign_buf = _mm256_movemask_epi8(coeffs_subzero);

      for (; scan_pos_sig >= sub_pos; scan_pos_sig--) {
        uint32_t id = scan_pos_sig - sub_pos;
        uint32_t shamt = (id << 1) + 1;

        uint32_t curr_sig = (sigs >> shamt) & 1;
        uint32_t curr_esc_flag = (esc_flags >> shamt) & 1;
        uint32_t curr_coeff_sign = (coeff_sign_buf >> shamt) & 1;

        uint32_t curr_pos_x = pos_xs_buf[id];
        uint32_t curr_pos_y = pos_ys_buf[id];

        if (curr_esc_flag | num_non_zero) {
          ctx_sig  = kvz_context_get_sig_ctx_inc(pattern_sig_ctx, scan_mode, curr_pos_x, curr_pos_y,
                                             log2_block_size, type);
          cabac->cur_ctx = &baseCtx[ctx_sig];
          CABAC_BIN(cabac, curr_sig, "sig_coeff_flag");
        }

        if (curr_sig) {
          abs_coeff[num_non_zero] = abs_coeff_buf[id];
          coeff_signs              = 2 * coeff_signs + curr_coeff_sign;
          num_non_zero++;
        }
      }
    } else {
      scan_pos_sig = sub_pos - 1;
    }

    if (num_non_zero > 0) {
      bool sign_hidden = last_nz_pos_in_cg - first_nz_pos_in_cg >= 4 /* SBH_THRESHOLD */
                         && !encoder->cfg.lossless;
      uint32_t ctx_set  = (i > 0 && type == 0) ? 2 : 0;
      cabac_ctx_t *base_ctx_mod;
      int32_t num_c1_flag, first_c2_flag_idx, idx, first_coeff2;

      if (c1 == 0) {
        ctx_set++;
      }

      c1 = 1;

      base_ctx_mod     = (type == 0) ? &(cabac->ctx.cu_one_model_luma[4 * ctx_set]) :
                         &(cabac->ctx.cu_one_model_chroma[4 * ctx_set]);
      num_c1_flag      = MIN(num_non_zero, C1FLAG_NUMBER);
      first_c2_flag_idx = -1;

      for (idx = 0; idx < num_c1_flag; idx++) {
        uint32_t symbol = (abs_coeff[idx] > 1) ? 1 : 0;
        cabac->cur_ctx = &base_ctx_mod[c1];
        CABAC_BIN(cabac, symbol, "coeff_abs_level_greater1_flag");

        if (symbol) {
          c1 = 0;

          if (first_c2_flag_idx == -1) {
            first_c2_flag_idx = idx;
          }
        } else if ((c1 < 3) && (c1 > 0)) {
          c1++;
        }
      }

      if (c1 == 0) {
        base_ctx_mod = (type == 0) ? &(cabac->ctx.cu_abs_model_luma[ctx_set]) :
                       &(cabac->ctx.cu_abs_model_chroma[ctx_set]);

        if (first_c2_flag_idx != -1) {
          uint8_t symbol = (abs_coeff[first_c2_flag_idx] > 2) ? 1 : 0;
          cabac->cur_ctx      = &base_ctx_mod[0];
          CABAC_BIN(cabac, symbol, "coeff_abs_level_greater2_flag");
        }
      }
      int32_t shiftamt = (be_valid && sign_hidden) ? 1 : 0;
      int32_t nnz = num_non_zero - shiftamt;
      coeff_signs >>= shiftamt;
      if (!cabac->only_count) {
        if (encoder->cfg.crypto_features & KVZ_CRYPTO_TRANSF_COEFF_SIGNS) {
          coeff_signs ^= kvz_crypto_get_key(state->crypto_hdl, nnz);
        }
      }
      CABAC_BINS_EP(cabac, coeff_signs, nnz, "coeff_sign_flag");

      if (c1 == 0 || num_non_zero > C1FLAG_NUMBER) {
        first_coeff2 = 1;

        for (idx = 0; idx < num_non_zero; idx++) {
          int32_t base_level  = (idx < C1FLAG_NUMBER) ? (2 + first_coeff2) : 1;

          if (abs_coeff[idx] >= base_level) {
            if (!cabac->only_count && (encoder->cfg.crypto_features & KVZ_CRYPTO_TRANSF_COEFFS)) {
              kvz_cabac_write_coeff_remain_encry(state, cabac, abs_coeff[idx] - base_level, go_rice_param, base_level);
            } else {
              kvz_cabac_write_coeff_remain(cabac, abs_coeff[idx] - base_level, go_rice_param);
            }

            if (abs_coeff[idx] > 3 * (1 << go_rice_param)) {
              go_rice_param = MIN(go_rice_param + 1, 4);
            }
          }

          if (abs_coeff[idx] >= 2) {
            first_coeff2 = 0;
          }
        }
      }
    }
    last_nz_pos_in_cg = -1;
    first_nz_pos_in_cg = 16;
    num_non_zero = 0;
    coeff_signs = 0;
  }
}

int kvz_strategy_register_encode_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "encode_coeff_nxn", "avx2", 40, &kvz_encode_coeff_nxn_avx2);

  return success;
}
