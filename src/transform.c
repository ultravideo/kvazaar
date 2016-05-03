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

#include "transform.h"

#include "image.h"
#include "kvazaar.h"
#include "rdo.h"
#include "strategies/strategies-dct.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-quant.h"
#include "tables.h"


//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
//


const uint8_t kvz_g_chroma_scale[58]=
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
int32_t kvz_get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset)
{
  int32_t qp_scaled = 0;
  if(type == 0) {
    qp_scaled = qp + qp_offset;
  } else {
    qp_scaled = CLIP(-qp_offset, 57, qp);
    if(qp_scaled < 0) {
      qp_scaled = qp_scaled + qp_offset;
    } else {
      qp_scaled = kvz_g_chroma_scale[qp_scaled] + qp_offset;
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
void kvz_transformskip(const encoder_control_t * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  kvz_g_convert_to_bit[block_size] + 2;
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
void kvz_itransformskip(const encoder_control_t * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  kvz_g_convert_to_bit[block_size] + 2;
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
void kvz_transform2d(const encoder_control_t * const encoder, int16_t *block, int16_t *coeff, int8_t block_size, int32_t mode)
{
  dct_func *dct_func = kvz_get_dct_func(block_size, mode);  
  dct_func(encoder->bitdepth, block, coeff);
}

void kvz_itransform2d(const encoder_control_t * const encoder, int16_t *block, int16_t *coeff, int8_t block_size, int32_t mode)
{
  dct_func *idct_func = kvz_get_idct_func(block_size, mode);
  idct_func(encoder->bitdepth, coeff, block);
}

/**
 * \brief Like kvz_quantize_residual except that this uses trskip if that is better.
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
int kvz_quantize_residual_trskip(
    encoder_state_t *const state,
    const cu_info_t *const cur_cu, const int width, const color_t color,
    const coeff_scan_order_t scan_order, int8_t *trskip_out, 
    const int in_stride, const int out_stride,
    const kvz_pixel *const ref_in, const kvz_pixel *const pred_in, 
    kvz_pixel *rec_out, coeff_t *coeff_out)
{
  struct {
    kvz_pixel rec[4*4];
    coeff_t coeff[4*4];
    uint32_t cost;
    int has_coeffs;
  } skip, noskip, *best;

  const int bit_cost = (int)(state->global->cur_lambda_cost+0.5);
  
  noskip.has_coeffs = kvz_quantize_residual(
      state, cur_cu, width, color, scan_order,
      0, in_stride, 4,
      ref_in, pred_in, noskip.rec, noskip.coeff);
  noskip.cost = kvz_pixels_calc_ssd(ref_in, noskip.rec, in_stride, 4, 4);
  noskip.cost += kvz_get_coeff_cost(state, noskip.coeff, 4, 0, scan_order) * bit_cost;

  skip.has_coeffs = kvz_quantize_residual(
    state, cur_cu, width, color, scan_order,
    1, in_stride, 4,
    ref_in, pred_in, skip.rec, skip.coeff);
  skip.cost = kvz_pixels_calc_ssd(ref_in, skip.rec, in_stride, 4, 4);
  skip.cost += kvz_get_coeff_cost(state, skip.coeff, 4, 0, scan_order) * bit_cost;

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
    kvz_pixels_blit(best->rec, rec_out, width, width, 4, out_stride);
  }
  kvz_coefficients_blit(best->coeff, coeff_out, width, width, 4, out_stride);

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
void kvz_quantize_lcu_luma_residual(encoder_state_t * const state, int32_t x, int32_t y, const uint8_t depth, cu_info_t *cur_cu, lcu_t* lcu)
{
  // we have 64>>depth transform size
  const vector2d_t lcu_px = { SUB_SCU(x), SUB_SCU(y) };
  const int pu_index = PU_INDEX(lcu_px.x / 4, lcu_px.y / 4);
  if (cur_cu == NULL) {
    cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  }
  const int8_t width = LCU_WIDTH>>depth;
  
  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);

  // Split transform and increase depth
  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;
    kvz_quantize_lcu_luma_residual(state, x,          y,          depth+1, NULL, lcu);
    kvz_quantize_lcu_luma_residual(state, x + offset, y,          depth+1, NULL, lcu);
    kvz_quantize_lcu_luma_residual(state, x,          y + offset, depth+1, NULL, lcu);
    kvz_quantize_lcu_luma_residual(state, x + offset, y + offset, depth+1, NULL, lcu);

    // Propagate coded block flags from child CUs to parent CU.
    if (depth < MAX_DEPTH) {
      cu_info_t *cu_a = LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y);
      cu_info_t *cu_b = LCU_GET_CU_AT_PX(lcu, lcu_px.x,          lcu_px.y + offset);
      cu_info_t *cu_c = LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y + offset);
      if (cbf_is_set(cu_a->cbf.y, depth+1) || cbf_is_set(cu_b->cbf.y, depth+1) || cbf_is_set(cu_c->cbf.y, depth+1)) {
        cbf_set(&cur_cu->cbf.y, depth);
      }
    }

    return;
  }

  {
    const int luma_offset = lcu_px.x + lcu_px.y * LCU_WIDTH;

    // Pointers to current location in arrays with prediction.
    kvz_pixel *recbase_y = &lcu->rec.y[luma_offset];
    // Pointers to current location in arrays with reference.
    const kvz_pixel *base_y = &lcu->ref.y[luma_offset];
    // Pointers to current location in arrays with kvantized coefficients.
    coeff_t *orig_coeff_y = &lcu->coeff.y[luma_offset];

    coeff_scan_order_t scan_idx_luma = kvz_get_scan_order(cur_cu->type, cur_cu->intra[pu_index].mode, depth);

    #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
    uint32_t residual_sum = 0;
    #endif

    // Clear coded block flag structures for depths lower than current depth.
    // This should ensure that the CBF data doesn't get corrupted if this function
    // is called more than once.
    cbf_clear(&cur_cu->cbf.y, depth + pu_index);

    if (width == 4 && 
        state->encoder_control->trskip_enable)
    {
      // Try quantization with trskip and use it if it's better.
      int has_coeffs = kvz_quantize_residual_trskip(
          state, cur_cu, width, COLOR_Y, scan_idx_luma,
          &cur_cu->intra[pu_index].tr_skip,
          LCU_WIDTH, LCU_WIDTH,
          base_y, recbase_y, recbase_y, orig_coeff_y
      );
      if (has_coeffs) {
        cbf_set(&cur_cu->cbf.y, depth + pu_index);
      }
    } else {
      int has_coeffs = kvz_quantize_residual(
          state, cur_cu, width, COLOR_Y, scan_idx_luma,
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


void kvz_quantize_lcu_chroma_residual(encoder_state_t * const state, int32_t x, int32_t y, const uint8_t depth, cu_info_t *cur_cu, lcu_t* lcu)
{
  // we have 64>>depth transform size
  const vector2d_t lcu_px = { SUB_SCU(x), SUB_SCU(y) };
  const int pu_index = PU_INDEX(lcu_px.x / 4, lcu_px.y / 4);
  const int8_t width = LCU_WIDTH>>depth;
  if (cur_cu == NULL) {
    cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  }
  
  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);

  // Split transform and increase depth
  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;
    kvz_quantize_lcu_chroma_residual(state, x,          y,          depth+1, NULL, lcu);
    kvz_quantize_lcu_chroma_residual(state, x + offset, y,          depth+1, NULL, lcu);
    kvz_quantize_lcu_chroma_residual(state, x,          y + offset, depth+1, NULL, lcu);
    kvz_quantize_lcu_chroma_residual(state, x + offset, y + offset, depth+1, NULL, lcu);

    // Propagate coded block flags from child CUs to parent CU.
    if (depth < MAX_DEPTH) {
      cu_info_t *cu_a = LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y);
      cu_info_t *cu_b = LCU_GET_CU_AT_PX(lcu, lcu_px.x,          lcu_px.y + offset);
      cu_info_t *cu_c = LCU_GET_CU_AT_PX(lcu, lcu_px.x + offset, lcu_px.y + offset);
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
    kvz_pixel *recbase_u = &lcu->rec.u[chroma_offset];
    kvz_pixel *recbase_v = &lcu->rec.v[chroma_offset];
    const kvz_pixel *base_u = &lcu->ref.u[chroma_offset];
    const kvz_pixel *base_v = &lcu->ref.v[chroma_offset];
    coeff_t *orig_coeff_u = &lcu->coeff.u[chroma_offset];
    coeff_t *orig_coeff_v = &lcu->coeff.v[chroma_offset];
    coeff_scan_order_t scan_idx_chroma;
    int tr_skip = 0;
    int chroma_depth = (depth == MAX_PU_DEPTH ? depth - 1 : depth);
    int chroma_width = LCU_WIDTH_C >> chroma_depth;

    scan_idx_chroma = kvz_get_scan_order(cur_cu->type, cur_cu->intra[0].mode_chroma, depth);
    if (kvz_quantize_residual(state, cur_cu, chroma_width, COLOR_U, scan_idx_chroma, tr_skip, LCU_WIDTH_C, LCU_WIDTH_C, base_u, recbase_u, recbase_u, orig_coeff_u)) {
      cbf_set(&cur_cu->cbf.u, depth);
    }
    if (kvz_quantize_residual(state, cur_cu, chroma_width, COLOR_V, scan_idx_chroma, tr_skip, LCU_WIDTH_C, LCU_WIDTH_C, base_v, recbase_v, recbase_v, orig_coeff_v)) {
      cbf_set(&cur_cu->cbf.v, depth);
    }
  }
}

