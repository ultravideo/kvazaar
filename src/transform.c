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

#include "transform.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "config.h"
#include "nal.h"
#include "rdo.h"
#include "strategies/strategies-dct.h"

//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
//


const uint8_t g_chroma_scale[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//

/**
 * \brief Get scaled QP used in quantization
 *
 */
int32_t get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset)
{
  int32_t qp_scaled = 0;
  if(type == 0) {
    qp_scaled = qp + qp_offset;
  } else {
    qp_scaled = CLIP(-qp_offset, 57, qp);
    if(qp_scaled < 0) {
      qp_scaled = qp_scaled + qp_offset;
    } else {
      qp_scaled = g_chroma_scale[qp_scaled] + qp_offset;
    }
  }
  return qp_scaled;
}

/**
 * \brief NxN inverse transform (2D)
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size input data (width of transform)
 */
void transformskip(const encoder_control * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  g_convert_to_bit[block_size] + 2;
  int32_t  shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size;
  int32_t  j,k;
  for (j = 0; j < block_size; j++) {
    for(k = 0; k < block_size; k ++) {
      coeff[j * block_size + k] = block[j * block_size + k] << shift;
    }
  }
}

/**
 * \brief inverse transform skip
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size width of transform
 */
void itransformskip(const encoder_control * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  g_convert_to_bit[block_size] + 2;
  int32_t  shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size;
  int32_t  j,k;
  int32_t offset;
  offset = (1 << (shift -1)); // For rounding
  for ( j = 0; j < block_size; j++ ) {
    for(k = 0; k < block_size; k ++) {
      block[j * block_size + k] =  (coeff[j * block_size + k] + offset) >> shift;
    }
  }
}

/**
 * \brief forward transform (2D)
 * \param block input residual
 * \param coeff transform coefficients
 * \param block_size width of transform
 */
void transform2d(const encoder_control * const encoder, int16_t *block, int16_t *coeff, int8_t block_size, int32_t mode)
{
  dct_func *dct_func = get_dct_func(block_size, mode);  
  dct_func(encoder->bitdepth, block, coeff);
}

void itransform2d(const encoder_control * const encoder, int16_t *block, int16_t *coeff, int8_t block_size, int32_t mode)
{
  dct_func *idct_func = get_idct_func(block_size, mode);
  idct_func(encoder->bitdepth, coeff, block);
}


#define QUANT_SHIFT 14
/**
 * \brief quantize transformed coefficents
 *
 */
void quant(const encoder_state_t * const encoder_state, int16_t *coef, int16_t *q_coef, int32_t width,
           int32_t height, int8_t type, int8_t scan_idx, int8_t block_type )
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  const uint32_t log2_block_size = g_convert_to_bit[ width ] + 2;
  const uint32_t * const scan = g_sig_last_scan[ scan_idx ][ log2_block_size - 1 ];

  int32_t qp_scaled = get_scaled_qp(type, encoder_state->global->QP, 0);

  const uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  const int32_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  const int32_t q_bits = QUANT_SHIFT + qp_scaled/6 + transform_shift;
  const int32_t add = ((encoder_state->global->slicetype == SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;

  uint32_t ac_sum = 0;

  for (int32_t n = 0; n < width * height; n++) {
    int32_t level;
    int32_t  sign;

    level = coef[n];
    sign  = (level < 0 ? -1: 1);

    level = ((int64_t)abs(level) * quant_coeff[n] + add) >> q_bits;
    ac_sum += level;

    level *= sign;
    q_coef[n] = (int16_t)(CLIP( -32768, 32767, level));
  }

  if (!(encoder->sign_hiding && ac_sum >= 2)) return;

  int32_t delta_u[LCU_WIDTH*LCU_WIDTH >> 2];

  for (int32_t n = 0; n < width * height; n++) {
    int32_t level;
    level = coef[n];
    level = ((int64_t)abs(level) * quant_coeff[n] + add) >> q_bits;
    delta_u[n] = (int32_t)(((int64_t)abs(coef[n]) * quant_coeff[n] - (level << q_bits)) >> q_bits8);
  }

  if(ac_sum >= 2) {
    #define SCAN_SET_SIZE 16
    #define LOG2_SCAN_SET_SIZE 4
    int32_t n,last_cg = -1, abssum = 0, subset, subpos;
    for(subset = (width*height - 1)>>LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg=-1;
      subpos = subset<<LOG2_SCAN_SET_SIZE;
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

      // Sum all quant coeffs between first and last
      for(n = first_nz_pos_in_cg; n <= last_nz_pos_in_cg; n++) {
        abssum += q_coef[scan[n + subpos]];
      }

      if(last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if(last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1) ;
        if(signbit != (abssum&0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff,  min_pos =-1, cur_cost=0x7fffffff;
          int16_t final_change = 0, cur_change=0;
          for(n = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n >= 0; n--) {
            uint32_t blkPos  = scan[n + subpos];
            if(q_coef[blkPos] != 0) {
              if(delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change=1;
              } else if(n == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost=0x7fffffff;
              } else {
                cur_cost = delta_u[blkPos];
                cur_change =-1;
              }
            } else if(n < first_nz_pos_in_cg && ((coef[blkPos] >= 0)?0:1) != signbit) {
              cur_cost = 0x7fffffff;
            } else {
              cur_cost   = -delta_u[blkPos];
              cur_change = 1;
            }

            if(cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos      = blkPos;
            }
          } // CG loop

          if(q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if(coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;
        } // Hide
      }
      if (last_cg == 1) last_cg=0;
    }

    #undef SCAN_SET_SIZE
    #undef LOG2_SCAN_SET_SIZE
  }
}

/**
 * \brief inverse quantize transformed and quantized coefficents
 *
 */
void dequant(const encoder_state_t * const encoder_state, int16_t *q_coef, int16_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  int32_t shift,add,coeff_q;
  int32_t n;
  int32_t transform_shift = 15 - encoder->bitdepth - (g_convert_to_bit[ width ] + 2);

  int32_t qp_scaled = get_scaled_qp(type, encoder_state->global->QP, 0);

  shift = 20 - QUANT_SHIFT - transform_shift;

  if (encoder->scaling_list.enable)
  {
    uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    const int32_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
    shift += 4;

    if (shift >qp_scaled / 6) {
      add = 1 << (shift - qp_scaled/6 - 1);

      for (n = 0; n < width * height; n++) {
        coeff_q = ((q_coef[n] * dequant_coef[n]) + add ) >> (shift -  qp_scaled/6);
        coef[n] = (int16_t)CLIP(-32768,32767,coeff_q);
      }
    } else {
      for (n = 0; n < width * height; n++) {
        // Clip to avoid possible overflow in following shift left operation
        coeff_q   = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q << (qp_scaled/6 - shift));
      }
    }
  } else {
    int32_t scale = g_inv_quant_scales[qp_scaled%6] << (qp_scaled/6);
    add = 1 << (shift-1);

    for (n = 0; n < width*height; n++) {
      coeff_q   = (q_coef[n] * scale + add) >> shift;
      coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}


/**
 * \brief Quantize residual and get both the reconstruction and coeffs.
 * 
 * \param width  Transform width.
 * \param color  Color.
 * \param scan_order  Coefficient scan order.
 * \param use_trskip  Whether transform skip is used.
 * \param stride  Stride for ref_in, pred_in rec_out and coeff_out.
 * \param ref_in  Reference pixels.
 * \param pred_in  Predicted pixels.
 * \param rec_out  Reconstructed pixels.
 * \param coeff_out  Coefficients used for reconstruction of rec_out.
 *
 * \returns  Whether coeff_out contains any non-zero coefficients.
 */
int quantize_residual(encoder_state_t *const encoder_state,
                      const cu_info *const cur_cu, const int width, const color_index color,
                      const coeff_scan_order_t scan_order, const int use_trskip, 
                      const int in_stride, const int out_stride,
                      const pixel *const ref_in, const pixel *const pred_in, 
                      pixel *rec_out, coefficient *coeff_out)
{
  // Temporary arrays to pass data to and from quant and transform functions.
  int16_t residual[TR_MAX_WIDTH * TR_MAX_WIDTH];
  coefficient quant_coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];
  coefficient coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

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
    transformskip(encoder_state->encoder_control, residual, coeff, width);
  } else {
    transform2d(encoder_state->encoder_control, residual, coeff, width, (color == COLOR_Y ? 0 : 65535));
  }

  // Quantize coeffs. (coeff -> quant_coeff)
  if (encoder_state->encoder_control->rdoq_enable) {
    int8_t tr_depth = cur_cu->tr_depth - cur_cu->depth;
    tr_depth += (cur_cu->part_size == SIZE_NxN ? 1 : 0);
    rdoq(encoder_state, coeff, quant_coeff, width, width, (color == COLOR_Y ? 0 : 2),
         scan_order, cur_cu->type, tr_depth);
  } else {
    quant(encoder_state, coeff, quant_coeff, width, width, (color == COLOR_Y ? 0 : 2),
          scan_order, cur_cu->type);
  }

  // Check if there are any non-zero coefficients.
  {
    int i;
    for (i = 0; i < width * width; ++i) {
      if (quant_coeff[i] != 0) {
        has_coeffs = 1;
        break;
      }
    }
  }

  // Copy coefficients to coeff_out.
  coefficients_blit(quant_coeff, coeff_out, width, width, width, out_stride);

  // Do the inverse quantization and transformation and the reconstruction to
  // rec_out.
  if (has_coeffs) {
    int y, x;

    // Get quantized residual. (quant_coeff -> coeff -> residual)
    dequant(encoder_state, quant_coeff, coeff, width, width, (color == COLOR_Y ? 0 : (color == COLOR_U ? 2 : 3)), cur_cu->type);
    if (use_trskip) {
      itransformskip(encoder_state->encoder_control, residual, coeff, width);
    } else {
      itransform2d(encoder_state->encoder_control, residual, coeff, width, (color == COLOR_Y ? 0 : 65535));
    }

    // Get quantized reconstruction. (residual + pred_in -> rec_out)
    for (y = 0; y < width; ++y) {
      for (x = 0; x < width; ++x) {
        int16_t val = residual[x + y * width] + pred_in[x + y * in_stride];
        rec_out[x + y * out_stride] = (uint8_t)CLIP(0, 255, val);
      }
    }
  } else if (rec_out != pred_in) {
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
 * \brief Like quantize_residual except that this uses trskip if that is better.
 *
 * Using this function saves one step of quantization and inverse quantization
 * compared to doing the decision separately from the actual operation.
 *
 * \param width  Transform width.
 * \param color  Color.
 * \param scan_order  Coefficient scan order.
 * \param trskip_out  Whether transform skip is used.
 * \param stride  Stride for ref_in, pred_in rec_out and coeff_out.
 * \param ref_in  Reference pixels.
 * \param pred_in  Predicted pixels.
 * \param rec_out  Reconstructed pixels.
 * \param coeff_out  Coefficients used for reconstruction of rec_out.
 *
 * \returns  Whether coeff_out contains any non-zero coefficients.
 */
int quantize_residual_trskip(
    encoder_state_t *const encoder_state,
    const cu_info *const cur_cu, const int width, const color_index color,
    const coeff_scan_order_t scan_order, int8_t *trskip_out, 
    const int in_stride, const int out_stride,
    const pixel *const ref_in, const pixel *const pred_in, 
    pixel *rec_out, coefficient *coeff_out)
{
  struct {
    pixel rec[4*4];
    coefficient coeff[4*4];
    uint32_t cost;
    int has_coeffs;
  } skip, noskip, *best;

  const int bit_cost = (int)(encoder_state->global->cur_lambda_cost+0.5);
  
  noskip.has_coeffs = quantize_residual(
      encoder_state, cur_cu, width, color, scan_order,
      0, in_stride, 4,
      ref_in, pred_in, noskip.rec, noskip.coeff);
  noskip.cost = pixels_calc_ssd(ref_in, noskip.rec, in_stride, 4, 4);
  noskip.cost += get_coeff_cost(encoder_state, noskip.coeff, 4, 0, scan_order) * bit_cost;

  skip.cost += get_coeff_cost(encoder_state, skip.coeff, 4, 0, scan_order) * bit_cost;

  if (noskip.cost <= skip.cost) {
    *trskip_out = 0;
    best = &noskip;
  } else {
    *trskip_out = 1;
    best = &skip;
  }

  if (best->has_coeffs || rec_out != pred_in) {
    // If there is no residual and reconstruction is already in rec_out, 
    // we can skip this.
    pixels_blit(best->rec, rec_out, width, width, 4, out_stride);
  }
  coefficients_blit(best->coeff, coeff_out, width, width, 4, out_stride);

  return best->has_coeffs;
}


/**
 * This function calculates the residual coefficients for a region of the LCU
 * (defined by x, y and depth) and updates the reconstruction with the
 * kvantized residual.
 *
 * It handles recursion for transform split, but that is currently only work
 * for 64x64 inter to 32x32 transform blocks.
 *
 * Inputs are:
 * - lcu->rec  pixels after prediction for the area
 * - lcu->ref  reference pixels for the area
 * - lcu->cu   for the area
 *
 * Outputs are:
 * - lcu->rec  reconstruction after quantized residual
 * - lcu->coeff  quantized coefficients for the area
 * - lcu->cbf  coded block flags for the area
 * - lcu->cu.intra[].tr_skip  for the area
 */
void quantize_lcu_luma_residual(encoder_state_t * const encoder_state, int32_t x, int32_t y, const uint8_t depth, cu_info *cur_cu, lcu_t* lcu)
{
  // we have 64>>depth transform size
  const vector2d lcu_px = {x & 0x3f, y & 0x3f};
  const int pu_index = PU_INDEX(lcu_px.x / 4, lcu_px.y / 4);
  if (cur_cu == NULL) {
    cur_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];
  }
  const int8_t width = LCU_WIDTH>>depth;
  
  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);

  // Split transform and increase depth
  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;
    quantize_lcu_luma_residual(encoder_state, x,          y,          depth+1, NULL, lcu);
    quantize_lcu_luma_residual(encoder_state, x + offset, y,          depth+1, NULL, lcu);
    quantize_lcu_luma_residual(encoder_state, x,          y + offset, depth+1, NULL, lcu);
    quantize_lcu_luma_residual(encoder_state, x + offset, y + offset, depth+1, NULL, lcu);

    // Propagate coded block flags from child CUs to parent CU.
    if (depth < MAX_DEPTH) {
      cu_info *cu_a =  &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset)>>3) +  (lcu_px.y>>3)        *LCU_T_CU_WIDTH];
      cu_info *cu_b =  &lcu->cu[LCU_CU_OFFSET +  (lcu_px.x>>3)           + ((lcu_px.y+offset)>>3)*LCU_T_CU_WIDTH];
      cu_info *cu_c =  &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset)>>3) + ((lcu_px.y+offset)>>3)*LCU_T_CU_WIDTH];
      if (cbf_is_set(cu_a->cbf.y, depth+1) || cbf_is_set(cu_b->cbf.y, depth+1) || cbf_is_set(cu_c->cbf.y, depth+1)) {
        cbf_set(&cur_cu->cbf.y, depth);
      }
    }

    return;
  }

  {
    const int luma_offset = lcu_px.x + lcu_px.y * LCU_WIDTH;

    // Pointers to current location in arrays with prediction.
    pixel *recbase_y = &lcu->rec.y[luma_offset];
    // Pointers to current location in arrays with reference.
    const pixel *base_y = &lcu->ref.y[luma_offset];
    // Pointers to current location in arrays with kvantized coefficients.
    coefficient *orig_coeff_y = &lcu->coeff.y[luma_offset];

    coeff_scan_order_t scan_idx_luma = get_scan_order(cur_cu->type, cur_cu->intra[pu_index].mode, depth);

    #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
    uint32_t residual_sum = 0;
    #endif

    // Clear coded block flag structures for depths lower than current depth.
    // This should ensure that the CBF data doesn't get corrupted if this function
    // is called more than once.
    cbf_clear(&cur_cu->cbf.y, depth + pu_index);

    if (width == 4 && 
        encoder_state->encoder_control->trskip_enable)
    {
      // Try quantization with trskip and use it if it's better.
      int has_coeffs = quantize_residual_trskip(
          encoder_state, cur_cu, width, COLOR_Y, scan_idx_luma,
          &cur_cu->intra[pu_index].tr_skip,
          LCU_WIDTH, LCU_WIDTH,
          base_y, recbase_y, recbase_y, orig_coeff_y
      );
      if (has_coeffs) {
        cbf_set(&cur_cu->cbf.y, depth + pu_index);
      }
    } else {
      int has_coeffs = quantize_residual(
          encoder_state, cur_cu, width, COLOR_Y, scan_idx_luma,
          0,
          LCU_WIDTH, LCU_WIDTH,
          base_y, recbase_y, recbase_y, orig_coeff_y
      );
      if (has_coeffs) {
        cbf_set(&cur_cu->cbf.y, depth + pu_index);
      }
    }
  }
}


void quantize_lcu_chroma_residual(encoder_state_t * const encoder_state, int32_t x, int32_t y, const uint8_t depth, cu_info *cur_cu, lcu_t* lcu)
{
  // we have 64>>depth transform size
  const vector2d lcu_px = {x & 0x3f, y & 0x3f};
  const int pu_index = PU_INDEX(lcu_px.x / 4, lcu_px.y / 4);
  const int8_t width = LCU_WIDTH>>depth;
  if (cur_cu == NULL) {
    cur_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];
  }
  
  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);

  // Split transform and increase depth
  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;
    quantize_lcu_chroma_residual(encoder_state, x,          y,          depth+1, NULL, lcu);
    quantize_lcu_chroma_residual(encoder_state, x + offset, y,          depth+1, NULL, lcu);
    quantize_lcu_chroma_residual(encoder_state, x,          y + offset, depth+1, NULL, lcu);
    quantize_lcu_chroma_residual(encoder_state, x + offset, y + offset, depth+1, NULL, lcu);

    // Propagate coded block flags from child CUs to parent CU.
    if (depth < MAX_DEPTH) {
      cu_info *cu_a =  &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset)>>3) +  (lcu_px.y>>3)        *LCU_T_CU_WIDTH];
      cu_info *cu_b =  &lcu->cu[LCU_CU_OFFSET +  (lcu_px.x>>3)           + ((lcu_px.y+offset)>>3)*LCU_T_CU_WIDTH];
      cu_info *cu_c =  &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset)>>3) + ((lcu_px.y+offset)>>3)*LCU_T_CU_WIDTH];
      if (cbf_is_set(cu_a->cbf.u, depth+1) || cbf_is_set(cu_b->cbf.u, depth+1) || cbf_is_set(cu_c->cbf.u, depth+1)) {
        cbf_set(&cur_cu->cbf.u, depth);
      }
      if (cbf_is_set(cu_a->cbf.v, depth+1) || cbf_is_set(cu_b->cbf.v, depth+1) || cbf_is_set(cu_c->cbf.v, depth+1)) {
        cbf_set(&cur_cu->cbf.v, depth);
      }
    }

    return;
  }

  // If luma is 4x4, do chroma for the 8x8 luma area when handling the top
  // left PU because the coordinates are correct.
  if (depth <= MAX_DEPTH || pu_index == 0) {
    cbf_clear(&cur_cu->cbf.u, depth);
    cbf_clear(&cur_cu->cbf.v, depth);

    const int chroma_offset = lcu_px.x / 2 + lcu_px.y / 2 * LCU_WIDTH_C;
    pixel *recbase_u = &lcu->rec.u[chroma_offset];
    pixel *recbase_v = &lcu->rec.v[chroma_offset];
    const pixel *base_u = &lcu->ref.u[chroma_offset];
    const pixel *base_v = &lcu->ref.v[chroma_offset];
    coefficient *orig_coeff_u = &lcu->coeff.u[chroma_offset];
    coefficient *orig_coeff_v = &lcu->coeff.v[chroma_offset];
    coeff_scan_order_t scan_idx_chroma;
    int tr_skip = 0;
    int chroma_depth = (depth == MAX_PU_DEPTH ? depth - 1 : depth);
    int chroma_width = LCU_WIDTH_C >> chroma_depth;

    scan_idx_chroma = get_scan_order(cur_cu->type, cur_cu->intra[0].mode_chroma, depth);
    if (quantize_residual(encoder_state, cur_cu, chroma_width, COLOR_U, scan_idx_chroma, tr_skip, LCU_WIDTH_C, LCU_WIDTH_C, base_u, recbase_u, recbase_u, orig_coeff_u)) {
      cbf_set(&cur_cu->cbf.u, depth);
    }
    if (quantize_residual(encoder_state, cur_cu, chroma_width, COLOR_V, scan_idx_chroma, tr_skip, LCU_WIDTH_C, LCU_WIDTH_C, base_v, recbase_v, recbase_v, orig_coeff_v)) {
      cbf_set(&cur_cu->cbf.v, depth);
    }
  }
}

