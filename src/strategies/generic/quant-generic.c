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

#include "strategies/generic/quant-generic.h"

#include <stdlib.h>

#include "encoder.h"
#include "rdo.h"
#include "scalinglist.h"
#include "strategies/strategies-quant.h"
#include "strategyselector.h"
#include "transform.h"
#include "fast_coeff_cost.h"

#define QUANT_SHIFT 14
/**
* \brief quantize transformed coefficents
*
*/
void kvz_quant_generic(const encoder_state_t * const state, coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t * const scan = kvz_g_sig_last_scan[scan_idx][log2_block_size - 1];

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6);
  const uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  const coeff_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  const int32_t q_bits = QUANT_SHIFT + qp_scaled / 6 + transform_shift;
  const int32_t add = ((state->frame->slicetype == KVZ_SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;

  uint32_t ac_sum = 0;

  for (int32_t n = 0; n < width * height; n++) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t  sign;

    sign = (level < 0 ? -1 : 1);

    int32_t curr_quant_coeff = quant_coeff[n];
    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    ac_sum += level;

    level *= sign;
    q_coef[n] = (coeff_t)(CLIP(-32768, 32767, level));

  }

  if (!encoder->cfg.signhide_enable || ac_sum < 2) return;

  int32_t delta_u[LCU_WIDTH*LCU_WIDTH >> 2];

  for (int32_t n = 0; n < width * height; n++) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t curr_quant_coeff = quant_coeff[n];

    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    delta_u[n] = (int32_t)((abs_level * curr_quant_coeff - (level << q_bits)) >> q_bits8);
  }

  if (ac_sum >= 2) {
#define SCAN_SET_SIZE 16
#define LOG2_SCAN_SET_SIZE 4
    int32_t n, last_cg = -1, abssum = 0, subset, subpos;
    for (subset = (width*height - 1) >> LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg = -1;
      subpos = subset << LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n = SCAN_SET_SIZE - 1; n >= 0; n--)  {
        if (q_coef[scan[n + subpos]])  {
          last_nz_pos_in_cg = n;
          break;
        }
      }

      // First coeff pos
      for (n = 0; n <SCAN_SET_SIZE; n++) {
        if (q_coef[scan[n + subpos]]) {
          first_nz_pos_in_cg = n;
          break;
        }
      }

      // Sum all kvz_quant coeffs between first and last
      for (n = first_nz_pos_in_cg; n <= last_nz_pos_in_cg; n++) {
        abssum += q_coef[scan[n + subpos]];
      }

      if (last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1);
        if (signbit != (abssum & 0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff, min_pos = -1, cur_cost = 0x7fffffff;
          int16_t final_change = 0, cur_change = 0;
          for (n = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n >= 0; n--) {
            uint32_t blkPos = scan[n + subpos];
            if (q_coef[blkPos] != 0) {
              if (delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change = 1;
              }
              else if (n == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost = 0x7fffffff;
              }
              else {
                cur_cost = delta_u[blkPos];
                cur_change = -1;
              }
            }
            else if (n < first_nz_pos_in_cg && ((coef[blkPos] >= 0) ? 0 : 1) != signbit) {
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

#undef SCAN_SET_SIZE
#undef LOG2_SCAN_SET_SIZE
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
* \param early_skip if this is used for early skip, bypass IT and IQ
*
* \returns  Whether coeff_out contains any non-zero coefficients.
*/
int kvz_quantize_residual_generic(encoder_state_t *const state,
  const cu_info_t *const cur_cu, const int width, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const kvz_pixel *const ref_in, const kvz_pixel *const pred_in,
  kvz_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip)
{
  // Temporary arrays to pass data to and from kvz_quant and transform functions.
  ALIGNED(64) int16_t residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  ALIGNED(64) coeff_t coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

  int has_coeffs = 0;

  assert(width <= TR_MAX_WIDTH);
  assert(width >= TR_MIN_WIDTH);

  // Get residual. (ref_in - pred_in -> residual)
  {
    int y, x;
    for (y = 0; y < width; ++y) {
      for (x = 0; x < width; ++x) {
        residual[x + y * width] = (int16_t)(ref_in[x + y * in_stride] - pred_in[x + y * in_stride]);
      }
    }
  }

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
  {
    int i;
    for (i = 0; i < width * width; ++i) {
      if (coeff_out[i] != 0) {
        has_coeffs = 1;
        break;
      }
    }
  }

  // Do the inverse quantization and transformation and the reconstruction to
  // rec_out.
  if (has_coeffs && !early_skip) {
    int y, x;

    // Get quantized residual. (coeff_out -> coeff -> residual)
    kvz_dequant(state, coeff_out, coeff, width, width, (color == COLOR_Y ? 0 : (color == COLOR_U ? 2 : 3)), cur_cu->type);
    if (use_trskip) {
      kvz_itransformskip(state->encoder_control, residual, coeff, width);
    }
    else {
      kvz_itransform2d(state->encoder_control, residual, coeff, width, color, cur_cu->type);
    }

    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    for (y = 0; y < width; ++y) {
      for (x = 0; x < width; ++x) {
        int16_t val = residual[x + y * width] + pred_in[x + y * in_stride];
        rec_out[x + y * out_stride] = (kvz_pixel)CLIP(0, PIXEL_MAX, val);
      }
    }
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

/**
 * \brief inverse quantize transformed and quantized coefficents
 *
 */
void kvz_dequant_generic(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type)
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

    const coeff_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
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

    for (n = 0; n < width*height; n++) {
      coeff_q   = (q_coef[n] * scale + add) >> shift;
      coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}

static uint32_t coeff_abs_sum_generic(const coeff_t *coeffs, size_t length)
{
  uint32_t sum = 0;
  for (int i = 0; i < length; i++) {
    sum += abs(coeffs[i]);
  }
  return sum;
}

static INLINE void get_coeff_weights(uint64_t wts_packed, uint16_t *weights)
{
  weights[0] = (wts_packed >>  0) & 0xffff;
  weights[1] = (wts_packed >> 16) & 0xffff;
  weights[2] = (wts_packed >> 32) & 0xffff;
  weights[3] = (wts_packed >> 48) & 0xffff;
}

static double fast_coeff_cost_generic(const coeff_t *coeff, int32_t width, uint64_t weights)
{
  uint32_t sum = 0;
  uint16_t weights_unpacked[4];

  get_coeff_weights(weights, weights_unpacked);

  for (int32_t i = 0; i < width * width; i++) {
     int16_t curr = coeff[i];
    uint32_t curr_abs = abs(curr);
    if (curr_abs > 3) {
      curr_abs = 3;
    }
    sum += weights_unpacked[curr_abs];
  }
  return (double) sum  / 256.0;
}



static void find_last_scanpos_generic(coeff_t* coef, coeff_t* dest_coeff, int8_t type, int32_t q_bits, const coeff_t* quant_coeff, struct kvz_sh_rates_t* sh_rates, const uint32_t cg_size, uint16_t* ctx_set, const uint32_t* scan, int32_t* cg_last_scanpos, int32_t* last_scanpos, uint32_t cg_num, int32_t* cg_scanpos, int32_t width, int8_t scan_mode) {
  for (*cg_scanpos = (cg_num - 1); *cg_scanpos >= 0; (*cg_scanpos)--) {
    for (int32_t scanpos_in_cg = (cg_size - 1); scanpos_in_cg >= 0; scanpos_in_cg--) {
      int32_t  scanpos = *cg_scanpos * cg_size + scanpos_in_cg;
      uint32_t blkpos = scan[scanpos];
      int32_t q = quant_coeff[blkpos];
      int32_t level_double = coef[blkpos];
      level_double = MIN(abs(level_double) * q, MAX_INT - (1 << (q_bits - 1)));
      uint32_t max_abs_level = (level_double + (1 << (q_bits - 1))) >> q_bits;

      if (max_abs_level > 0) {
        *last_scanpos = scanpos;
        *ctx_set = (scanpos > 0 && type == 0) ? 2 : 0;
        *cg_last_scanpos = *cg_scanpos;
        sh_rates->sig_coeff_inc[blkpos] = 0;
        return;
      }
      dest_coeff[blkpos] = 0;
    }
  }
}


int kvz_strategy_register_quant_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "quant", "generic", 0, &kvz_quant_generic);
  success &= kvz_strategyselector_register(opaque, "quantize_residual", "generic", 0, &kvz_quantize_residual_generic);
  success &= kvz_strategyselector_register(opaque, "dequant", "generic", 0, &kvz_dequant_generic);
  success &= kvz_strategyselector_register(opaque, "coeff_abs_sum", "generic", 0, &coeff_abs_sum_generic);
  success &= kvz_strategyselector_register(opaque, "fast_coeff_cost", "generic", 0, &fast_coeff_cost_generic);
  success &= kvz_strategyselector_register(opaque, "find_last_scanpos", "generic", 0, &find_last_scanpos_generic);

  return success;
}
