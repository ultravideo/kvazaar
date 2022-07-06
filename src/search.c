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

#include "search.h"

#include <limits.h>
#include <string.h>

#include "cabac.h"
#include "encoder.h"
#include "encode_coding_tree.h"
#include "imagelist.h"
#include "inter.h"
#include "intra.h"
#include "kvazaar.h"
#include "rdo.h"
#include "search_inter.h"
#include "search_intra.h"
#include "threadqueue.h"
#include "transform.h"
#include "videoframe.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-quant.h"


#define IN_FRAME(x, y, width, height, block_width, block_height) \
  ((x) >= 0 && (y) >= 0 \
  && (x) + (block_width) <= (width) \
  && (y) + (block_height) <= (height))

// Cost threshold for doing intra search in inter frames with --rd=0.
static const int INTRA_THRESHOLD = 8;


static INLINE void copy_cu_info(int x_local, int y_local, int width, lcu_t *from, lcu_t *to)
{
  for   (int y = y_local; y < y_local + width; y += SCU_WIDTH) {
    for (int x = x_local; x < x_local + width; x += SCU_WIDTH) {
      *LCU_GET_CU_AT_PX(to, x, y) = *LCU_GET_CU_AT_PX(from, x, y);
    }
  }
}

static INLINE void copy_cu_pixels(int x_local, int y_local, int width, lcu_t *from, lcu_t *to)
{
  const int luma_index = x_local + y_local * LCU_WIDTH;
  const int chroma_index = (x_local / 2) + (y_local / 2) * (LCU_WIDTH / 2);

  kvz_pixels_blit(&from->rec.y[luma_index], &to->rec.y[luma_index],
                  width, width, LCU_WIDTH, LCU_WIDTH);
  if (from->rec.chroma_format != KVZ_CSP_400) {
    kvz_pixels_blit(&from->rec.u[chroma_index], &to->rec.u[chroma_index],
                    width / 2, width / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    kvz_pixels_blit(&from->rec.v[chroma_index], &to->rec.v[chroma_index],
                    width / 2, width / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
  }
}

static INLINE void copy_cu_coeffs(int x_local, int y_local, int width, lcu_t *from, lcu_t *to)
{
  const int luma_z = xy_to_zorder(LCU_WIDTH, x_local, y_local);
  copy_coeffs(&from->coeff.y[luma_z], &to->coeff.y[luma_z], width);

  if (from->rec.chroma_format != KVZ_CSP_400) {
    const int chroma_z = xy_to_zorder(LCU_WIDTH_C, x_local >> 1, y_local >> 1);
    copy_coeffs(&from->coeff.u[chroma_z], &to->coeff.u[chroma_z], width >> 1);
    copy_coeffs(&from->coeff.v[chroma_z], &to->coeff.v[chroma_z], width >> 1);
  }
}

/**
 * Copy all non-reference CU data from next level to current level.
 */
static void work_tree_copy_up(int x_local, int y_local, int depth, lcu_t *work_tree)
{
  const int width = LCU_WIDTH >> depth;
  copy_cu_info  (x_local, y_local, width, &work_tree[depth + 1], &work_tree[depth]);
  copy_cu_pixels(x_local, y_local, width, &work_tree[depth + 1], &work_tree[depth]);
  copy_cu_coeffs(x_local, y_local, width, &work_tree[depth + 1], &work_tree[depth]);
}


/**
 * Copy all non-reference CU data from current level to all lower levels.
 */
static void work_tree_copy_down(int x_local, int y_local, int depth, lcu_t *work_tree)
{
  const int width = LCU_WIDTH >> depth;
  for (int i = depth + 1; i <= MAX_PU_DEPTH; i++) {
    copy_cu_info  (x_local, y_local, width, &work_tree[depth], &work_tree[i]);
    copy_cu_pixels(x_local, y_local, width, &work_tree[depth], &work_tree[i]);
  }
}

void kvz_lcu_fill_trdepth(lcu_t *lcu, int x_px, int y_px, int depth, int tr_depth)
{
  const int x_local = SUB_SCU(x_px);
  const int y_local = SUB_SCU(y_px);
  const int width = LCU_WIDTH >> depth;

  for (unsigned y = 0; y < width; y += SCU_WIDTH) {
    for (unsigned x = 0; x < width; x += SCU_WIDTH) {
      LCU_GET_CU_AT_PX(lcu, x_local + x, y_local + y)->tr_depth = tr_depth;
    }
  }
}

static void lcu_fill_cu_info(lcu_t *lcu, int x_local, int y_local, int width, int height, cu_info_t *cu)
{
  // Set mode in every CU covered by part_mode in this depth.
  for (int y = y_local; y < y_local + height; y += SCU_WIDTH) {
    for (int x = x_local; x < x_local + width; x += SCU_WIDTH) {
      cu_info_t *to = LCU_GET_CU_AT_PX(lcu, x, y);
      to->type      = cu->type;
      to->depth     = cu->depth;
      to->part_size = cu->part_size;
      to->qp        = cu->qp;

      if (cu->type == CU_INTRA) {
        to->intra.mode        = cu->intra.mode;
        to->intra.mode_chroma = cu->intra.mode_chroma;
      } else {
        to->skipped   = cu->skipped;
        to->merged    = cu->merged;
        to->merge_idx = cu->merge_idx;
        to->inter     = cu->inter;
      }
    }
  }
}

static void lcu_fill_inter(lcu_t *lcu, int x_local, int y_local, int cu_width)
{
  const part_mode_t part_mode = LCU_GET_CU_AT_PX(lcu, x_local, y_local)->part_size;
  const int num_pu = kvz_part_mode_num_parts[part_mode];

  for (int i = 0; i < num_pu; ++i) {
    const int x_pu      = PU_GET_X(part_mode, cu_width, x_local, i);
    const int y_pu      = PU_GET_Y(part_mode, cu_width, y_local, i);
    const int width_pu  = PU_GET_W(part_mode, cu_width, i);
    const int height_pu = PU_GET_H(part_mode, cu_width, i);

    cu_info_t *pu  = LCU_GET_CU_AT_PX(lcu, x_pu, y_pu);
    pu->type = CU_INTER;
    lcu_fill_cu_info(lcu, x_pu, y_pu, width_pu, height_pu, pu);
  }
}

static void lcu_fill_cbf(lcu_t *lcu, int x_local, int y_local, int width, cu_info_t *cur_cu)
{
  const uint32_t tr_split = cur_cu->tr_depth - cur_cu->depth;
  const uint32_t mask = ~((width >> tr_split)-1);

  // Set coeff flags in every CU covered by part_mode in this depth.
  for (uint32_t y = y_local; y < y_local + width; y += SCU_WIDTH) {
    for (uint32_t x = x_local; x < x_local + width; x += SCU_WIDTH) {
      // Use TU top-left CU to propagate coeff flags
      cu_info_t *cu_from = LCU_GET_CU_AT_PX(lcu, x & mask, y & mask);
      cu_info_t *cu_to   = LCU_GET_CU_AT_PX(lcu, x, y);
      if (cu_from != cu_to) {
        // Chroma coeff data is not used, luma is needed for deblocking
        cbf_copy(&cu_to->cbf, cu_from->cbf, COLOR_Y);
      }
    }
  }
}


//Calculates cost for all zero coeffs
static double cu_zero_coeff_cost(const encoder_state_t *state, lcu_t *work_tree, const int x, const int y,
  const int depth)
{
  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);
  int cu_width = LCU_WIDTH >> depth;
  lcu_t *const lcu = &work_tree[depth];

  const int luma_index = y_local * LCU_WIDTH + x_local;
  const int chroma_index = (y_local / 2) * LCU_WIDTH_C + (x_local / 2);

  double ssd = 0.0;
  ssd += KVZ_LUMA_MULT * kvz_pixels_calc_ssd(
    &lcu->ref.y[luma_index], &lcu->rec.y[luma_index],
    LCU_WIDTH, LCU_WIDTH, cu_width
    );
  if (x % 8 == 0 && y % 8 == 0 && state->encoder_control->chroma_format != KVZ_CSP_400) {
    ssd += KVZ_CHROMA_MULT * kvz_pixels_calc_ssd(
      &lcu->ref.u[chroma_index], &lcu->rec.u[chroma_index],
      LCU_WIDTH_C, LCU_WIDTH_C, cu_width / 2
      );
    ssd += KVZ_CHROMA_MULT * kvz_pixels_calc_ssd(
      &lcu->ref.v[chroma_index], &lcu->rec.v[chroma_index],
      LCU_WIDTH_C, LCU_WIDTH_C, cu_width / 2
      );
  }
  // Save the pixels at a lower level of the working tree.
  copy_cu_pixels(x_local, y_local, cu_width, lcu, &work_tree[depth + 1]);

  return ssd;
}


/**
* Calculate RD cost for a Coding Unit.
* \return Cost of block
* \param ref_cu  CU used for prediction parameters.
*
* Calculates the RDO cost of a single CU that will not be split further.
* Takes into account SSD of reconstruction and the cost of encoding whatever
* prediction unit data needs to be coded.
*/
double kvz_cu_rd_cost_luma(const encoder_state_t *const state,
                           const int x_px, const int y_px, const int depth,
                           const cu_info_t *const pred_cu,
                           lcu_t *const lcu)
{
  const int width = LCU_WIDTH >> depth;
  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type == CU_INTER && pred_cu->cbf == 0);

  // cur_cu is used for TU parameters.
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, x_px, y_px);

  double coeff_bits = 0;
  double tr_tree_bits = 0;

  // Check that lcu is not in 
  assert(x_px >= 0 && x_px < LCU_WIDTH);
  assert(y_px >= 0 && y_px < LCU_WIDTH);

  const uint8_t tr_depth = tr_cu->tr_depth - depth;

  cabac_data_t* cabac = (cabac_data_t *)&state->search_cabac;

  // Add transform_tree split_transform_flag bit cost.
  bool intra_split_flag = pred_cu->type == CU_INTRA && pred_cu->part_size == SIZE_NxN && depth == 3;
  int max_tr_depth;
  if (pred_cu->type == CU_INTRA) {
    max_tr_depth = state->encoder_control->cfg.tr_depth_intra + intra_split_flag;
  }
  else {
    max_tr_depth = state->encoder_control->tr_depth_inter;
  }
  if (width <= TR_MAX_WIDTH
      && width > TR_MIN_WIDTH
      && !intra_split_flag
      && MIN(tr_cu->tr_depth, depth) - tr_cu->depth < max_tr_depth
      && !skip_residual_coding)
  {
    cabac_ctx_t *ctx = &(cabac->ctx.trans_subdiv_model[5 - (6 - depth)]);
    CABAC_FBITS_UPDATE(cabac, ctx, tr_depth > 0, tr_tree_bits, "tr_split_search");
  }

  if (tr_depth > 0) {
    int offset = width / 2;
    double sum = 0;

    sum += kvz_cu_rd_cost_luma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * state->lambda;
  }


  if (cabac->update && tr_cu->tr_depth == tr_cu->depth && !skip_residual_coding) {
    // Because these need to be coded before the luma cbf they also need to be counted
    // before the cabac state changes. However, since this branch is only executed when
    // calculating the last RD cost it is not problem to include the chroma cbf costs in
    // luma, because the chroma cost is calculated right after the luma cost.
    // However, if we have different tr_depth, the bits cannot be written in correct
    // order anyways so do not touch the chroma cbf here.
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      cabac_ctx_t* cr_ctx = &(cabac->ctx.qt_cbf_model_chroma[depth - tr_cu->depth]);
      cabac->cur_ctx = cr_ctx;
      int u_is_set = cbf_is_set(pred_cu->cbf, depth, COLOR_U);
      int v_is_set = cbf_is_set(pred_cu->cbf, depth, COLOR_V);
      CABAC_FBITS_UPDATE(cabac, cr_ctx, u_is_set, tr_tree_bits, "cbf_cb_search");
      CABAC_FBITS_UPDATE(cabac, cr_ctx, v_is_set, tr_tree_bits, "cbf_cb_search");
    }
  }

  // Add transform_tree cbf_luma bit cost.
  const int is_tr_split = tr_cu->tr_depth - tr_cu->depth;
  int is_set = cbf_is_set(pred_cu->cbf, depth, COLOR_Y);
  if (pred_cu->type == CU_INTRA ||
      is_tr_split ||
      cbf_is_set(tr_cu->cbf, depth, COLOR_U) ||
      cbf_is_set(tr_cu->cbf, depth, COLOR_V))
  {
    cabac_ctx_t *ctx = &(cabac->ctx.qt_cbf_model_luma[!is_tr_split]);

    CABAC_FBITS_UPDATE(cabac, ctx, is_set, tr_tree_bits, "cbf_y_search");
  }

  // SSD between reconstruction and original
  int ssd = 0;
  if (!state->encoder_control->cfg.lossless) {
    int index = y_px * LCU_WIDTH + x_px;
    ssd = kvz_pixels_calc_ssd(&lcu->ref.y[index], &lcu->rec.y[index],
                                        LCU_WIDTH,          LCU_WIDTH,
                                        width);
  }


  if (!skip_residual_coding) {
    int8_t luma_scan_mode = kvz_get_scan_order(pred_cu->type, pred_cu->intra.mode, depth);
    const coeff_t *coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, x_px, y_px)];

    if(is_set)
      coeff_bits += kvz_get_coeff_cost(state, coeffs, width, 0, luma_scan_mode);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * KVZ_LUMA_MULT + bits * state->lambda;
}


double kvz_cu_rd_cost_chroma(const encoder_state_t *const state,
                             const int x_px, const int y_px, const int depth,
                             const cu_info_t *const pred_cu,
                             lcu_t *const lcu)
{
  const vector2d_t lcu_px = { x_px / 2, y_px / 2 };
  const int width = (depth <= MAX_DEPTH) ? LCU_WIDTH >> (depth + 1) : LCU_WIDTH >> depth;
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, x_px, y_px);
  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type == CU_INTER && pred_cu->cbf == 0);

  double tr_tree_bits = 0;
  double coeff_bits = 0;

  assert(x_px >= 0 && x_px < LCU_WIDTH);
  assert(y_px >= 0 && y_px < LCU_WIDTH);

  if (x_px % 8 != 0 || y_px % 8 != 0) {
    // For MAX_PU_DEPTH calculate chroma for previous depth for the first
    // block and return 0 cost for all others.
    return 0;
  }

  int u_is_set = cbf_is_set(pred_cu->cbf, depth, COLOR_U);
  int v_is_set = cbf_is_set(pred_cu->cbf, depth, COLOR_V);
  // See luma for why the second condition
  if (depth < MAX_PU_DEPTH && (!state->search_cabac.update || tr_cu->tr_depth != tr_cu->depth) && !skip_residual_coding) {
    const int tr_depth = depth - pred_cu->depth;
    cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;
    cabac_ctx_t *ctx = &(cabac->ctx.qt_cbf_model_chroma[tr_depth]);
    cabac->cur_ctx = ctx;
    if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_U)) {
      CABAC_FBITS_UPDATE(cabac, ctx, u_is_set, tr_tree_bits, "cbf_cb_search");
    }
    if (tr_depth == 0 || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_V)) {
      CABAC_FBITS_UPDATE(cabac, ctx, v_is_set, tr_tree_bits, "cbf_cb_search");
    }
  }

  if (tr_cu->tr_depth > depth) {
    int offset = LCU_WIDTH >> (depth + 1);
    double sum = 0;

    sum += kvz_cu_rd_cost_chroma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * state->lambda;
  }

  // Chroma SSD
  int ssd = 0;
  if (!state->encoder_control->cfg.lossless) {
    int index = lcu_px.y * LCU_WIDTH_C + lcu_px.x;
    int ssd_u = kvz_pixels_calc_ssd(&lcu->ref.u[index], &lcu->rec.u[index],
                                    LCU_WIDTH_C,         LCU_WIDTH_C,
                                    width);
    int ssd_v = kvz_pixels_calc_ssd(&lcu->ref.v[index], &lcu->rec.v[index],
                                    LCU_WIDTH_C,        LCU_WIDTH_C,
                                    width);
    ssd = ssd_u + ssd_v;
  }

  if (!skip_residual_coding)
  {
    int8_t scan_order = kvz_get_scan_order(pred_cu->type, pred_cu->intra.mode_chroma, depth);
    const int index = xy_to_zorder(LCU_WIDTH_C, lcu_px.x, lcu_px.y);

    if(u_is_set)coeff_bits += kvz_get_coeff_cost(state, &lcu->coeff.u[index], width, 2, scan_order);
    if(v_is_set)coeff_bits += kvz_get_coeff_cost(state, &lcu->coeff.v[index], width, 2, scan_order);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * KVZ_CHROMA_MULT + bits * state->lambda;
}

static double cu_rd_cost_tr_split_accurate(const encoder_state_t* const state,
                                           const int x_px, const int y_px, const int depth,
                                           const cu_info_t* const pred_cu,
                                           lcu_t* const lcu) {
  const int width = LCU_WIDTH >> depth;

  const int skip_residual_coding = pred_cu->skipped || (pred_cu->type == CU_INTER && pred_cu->cbf == 0);
  // cur_cu is used for TU parameters.
  cu_info_t* const tr_cu = LCU_GET_CU_AT_PX(lcu, x_px, y_px);

  double coeff_bits = 0;
  double tr_tree_bits = 0;

  // Check that lcu is not in 
  assert(x_px >= 0 && x_px < LCU_WIDTH);
  assert(y_px >= 0 && y_px < LCU_WIDTH);

  const uint8_t tr_depth = tr_cu->tr_depth - depth;

  const int cb_flag_u = cbf_is_set(tr_cu->cbf, depth, COLOR_U);
  const int cb_flag_v = cbf_is_set(tr_cu->cbf, depth, COLOR_V);

  cabac_data_t* cabac = (cabac_data_t*)&state->search_cabac;

  {
    int cbf = cbf_is_set_any(pred_cu->cbf, depth);
    // Only need to signal coded block flag if not skipped or merged
    // skip = no coded residual, merge = coded residual
    if (pred_cu->type == CU_INTER && (pred_cu->part_size != SIZE_2Nx2N || !pred_cu->merged)) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_qt_root_cbf_model), cbf, tr_tree_bits, "rqt_root_cbf");
    }

  }
  // Add transform_tree split_transform_flag bit cost.
  bool intra_split_flag = pred_cu->type == CU_INTRA && pred_cu->part_size == SIZE_NxN && depth == 3;
  int max_tr_depth;
  if (pred_cu->type == CU_INTRA) {
    max_tr_depth = state->encoder_control->cfg.tr_depth_intra + intra_split_flag;
  }
  else {
    max_tr_depth = state->encoder_control->tr_depth_inter;
  }
  if (width <= TR_MAX_WIDTH
    && width > TR_MIN_WIDTH
    && !intra_split_flag
    && MIN(tr_cu->tr_depth, depth) - tr_cu->depth < max_tr_depth
    && !skip_residual_coding)
  {
    cabac_ctx_t* ctx = &(cabac->ctx.trans_subdiv_model[5 - (6 - depth)]);
    CABAC_FBITS_UPDATE(cabac, ctx, tr_depth > 0, tr_tree_bits, "tr_split_search");
  }

  if(state->encoder_control->chroma_format != KVZ_CSP_400 && !skip_residual_coding) {
    if(tr_cu->depth == depth || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_U)) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_chroma[depth - tr_cu->depth]), cb_flag_u, tr_tree_bits, "cbf_cb");
    } 
    if(tr_cu->depth == depth || cbf_is_set(pred_cu->cbf, depth - 1, COLOR_V)) {
      CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.qt_cbf_model_chroma[depth - tr_cu->depth]), cb_flag_v, tr_tree_bits, "cbf_cr");
    } 
  }

  if (tr_depth > 0) {
    int offset = LCU_WIDTH >> (depth + 1);
    double sum = 0;

    sum += cu_rd_cost_tr_split_accurate(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_tr_split_accurate(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_tr_split_accurate(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_tr_split_accurate(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);
    return sum + tr_tree_bits * state->lambda;
  }
  const int cb_flag_y = cbf_is_set(tr_cu->cbf, depth, COLOR_Y) ;

  // Add transform_tree cbf_luma bit cost.
  const int is_tr_split = depth - tr_cu->depth;
  if ((pred_cu->type == CU_INTRA ||
    is_tr_split ||
    cb_flag_u ||
    cb_flag_v) 
      && !skip_residual_coding)
  {
    cabac_ctx_t* ctx = &(cabac->ctx.qt_cbf_model_luma[!is_tr_split]);

    CABAC_FBITS_UPDATE(cabac, ctx, cb_flag_y, tr_tree_bits, "cbf_y_search");
  }
  // SSD between reconstruction and original
  unsigned luma_ssd = 0;
  if (!state->encoder_control->cfg.lossless) {
    int index = y_px * LCU_WIDTH + x_px;
    luma_ssd = kvz_pixels_calc_ssd(&lcu->ref.y[index], &lcu->rec.y[index],
      LCU_WIDTH, LCU_WIDTH,
      width);
  }

  {
    int8_t luma_scan_mode = kvz_get_scan_order(pred_cu->type, pred_cu->intra.mode, depth);
    const coeff_t* coeffs = &lcu->coeff.y[xy_to_zorder(LCU_WIDTH, x_px, y_px)];

    if(cb_flag_y)
      coeff_bits += kvz_get_coeff_cost(state, coeffs, width, 0, luma_scan_mode);
  }

  unsigned chroma_ssd = 0;
  if(state->encoder_control->chroma_format != KVZ_CSP_400 && x_px % 8 == 0 && y_px % 8 == 0) {
    const vector2d_t lcu_px = { x_px / 2, y_px / 2 };
    const int chroma_width = (depth <= MAX_DEPTH) ? LCU_WIDTH >> (depth + 1) : LCU_WIDTH >> depth;
    if (!state->encoder_control->cfg.lossless) {
      int index = lcu_px.y * LCU_WIDTH_C + lcu_px.x;
      unsigned ssd_u = kvz_pixels_calc_ssd(&lcu->ref.u[index], &lcu->rec.u[index],
        LCU_WIDTH_C, LCU_WIDTH_C,
        chroma_width);
      unsigned ssd_v = kvz_pixels_calc_ssd(&lcu->ref.v[index], &lcu->rec.v[index],
        LCU_WIDTH_C, LCU_WIDTH_C,
        chroma_width);
      chroma_ssd = ssd_u + ssd_v;
    }

     {
      int8_t scan_order = kvz_get_scan_order(pred_cu->type, pred_cu->intra.mode_chroma, depth);
      const unsigned index = xy_to_zorder(LCU_WIDTH_C, lcu_px.x, lcu_px.y);

      if(cb_flag_u)coeff_bits += kvz_get_coeff_cost(state, &lcu->coeff.u[index], chroma_width, 2, scan_order);
      if (cb_flag_v)coeff_bits += kvz_get_coeff_cost(state, &lcu->coeff.v[index], chroma_width, 2, scan_order);
    }
  }

  double bits = tr_tree_bits + coeff_bits;
  return luma_ssd * KVZ_LUMA_MULT + chroma_ssd * KVZ_CHROMA_MULT + bits * state->lambda;
}


// Return estimate of bits used to code prediction mode of cur_cu.
static double calc_mode_bits(const encoder_state_t *state,
                             const lcu_t *lcu,
                             const cu_info_t * cur_cu,
                             int x, int y)
{
  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);

  assert(cur_cu->type == CU_INTRA);

  int8_t candidate_modes[3];
  {
    const cu_info_t *left_cu  = ((x >= SCU_WIDTH) ? LCU_GET_CU_AT_PX(lcu, x_local - SCU_WIDTH, y_local) : NULL);
    const cu_info_t *above_cu = ((y >= SCU_WIDTH) ? LCU_GET_CU_AT_PX(lcu, x_local, y_local - SCU_WIDTH) : NULL);
    kvz_intra_get_dir_luma_predictor(x, y, candidate_modes, cur_cu, left_cu, above_cu);
  }

  double mode_bits = kvz_luma_mode_bits(state, cur_cu->intra.mode, candidate_modes);

  if (x % 8 == 0 && y % 8 == 0 && state->encoder_control->chroma_format != KVZ_CSP_400) {
    mode_bits += kvz_chroma_mode_bits(state, cur_cu->intra.mode_chroma, cur_cu->intra.mode);
  }

  return mode_bits;
}


// TODO: replace usages of this by the kvz_sort_indices_by_cost function.
/**
 * \brief Sort modes and costs to ascending order according to costs.
 */
void kvz_sort_modes(int8_t *__restrict modes, double *__restrict costs, uint8_t length)
{
  // Length for intra is always between 5 and 23, and is either 21, 17, 9 or 8 about
  // 60% of the time, so there should be no need for anything more complex
  // than insertion sort.
  // Length for merge is 5 or less.
  for (uint8_t i = 1; i < length; ++i) {
    const double cur_cost = costs[i];
    const int8_t cur_mode = modes[i];
    uint8_t j = i;
    while (j > 0 && cur_cost < costs[j - 1]) {
      costs[j] = costs[j - 1];
      modes[j] = modes[j - 1];
      --j;
    }
    costs[j] = cur_cost;
    modes[j] = cur_mode;
  }
}


/**
 * \brief Sort keys (indices) to ascending order according to costs.
 */
void kvz_sort_keys_by_cost(unit_stats_map_t *__restrict map)
{
  // Size of sorted arrays is expected to be "small". No need for faster algorithm.
  for (uint8_t i = 1; i < map->size; ++i) {
    const int8_t cur_indx = map->keys[i];
    const double cur_cost = map->cost[cur_indx];
    uint8_t j = i;
    while (j > 0 && cur_cost < map->cost[map->keys[j - 1]]) {
      map->keys[j] = map->keys[j - 1];
      --j;
    }
    map->keys[j] = cur_indx;
  }
}


static uint8_t get_ctx_cu_split_model(const lcu_t *lcu, int x, int y, int depth)
{
  vector2d_t lcu_cu = { SUB_SCU(x), SUB_SCU(y) };
  bool condA = x >= 8 && LCU_GET_CU_AT_PX(lcu, lcu_cu.x - 1, lcu_cu.y    )->depth > depth;
  bool condL = y >= 8 && LCU_GET_CU_AT_PX(lcu, lcu_cu.x,     lcu_cu.y - 1)->depth > depth;
  return condA + condL;
}

/**
 * Search every mode from 0 to MAX_PU_DEPTH and return cost of best mode.
 * - The recursion is started at depth 0 and goes in Z-order to MAX_PU_DEPTH.
 * - Data structure work_tree is maintained such that the neighbouring SCUs
 *   and pixels to the left and up of current CU are the final CUs decided
 *   via the search. This is done by copying the relevant data to all
 *   relevant levels whenever a decision is made whether to split or not.
 * - All the final data for the LCU gets eventually copied to depth 0, which
 *   will be the final output of the recursion.
 */
static double search_cu(encoder_state_t * const state, int x, int y, int depth, lcu_t *work_tree)
{
  const encoder_control_t* ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  int cu_width = LCU_WIDTH >> depth;
  double cost = MAX_DOUBLE;
  double inter_zero_coeff_cost = MAX_DOUBLE;
  double inter_bitcost = MAX_INT;
  cu_info_t *cur_cu;
  cabac_data_t pre_search_cabac;
  memcpy(&pre_search_cabac, &state->search_cabac, sizeof(pre_search_cabac));

  struct {
    int32_t min;
    int32_t max;
  } pu_depth_inter, pu_depth_intra;

  lcu_t *const lcu = &work_tree[depth];

  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);

  // Stop recursion if the CU is completely outside the frame.
  if (x >= frame->width || y >= frame->height) {
    // Return zero cost because this CU does not have to be coded.
    return 0;
  }

  int gop_layer = ctrl->cfg.gop_len != 0 ? ctrl->cfg.gop[state->frame->gop_offset].layer - 1 : 0;

  // Assign correct depth limit
  constraint_t* constr = state->constraint;
  if(constr->ml_intra_depth_ctu) {
    pu_depth_intra.min = constr->ml_intra_depth_ctu->_mat_upper_depth[(x_local >> 3) + (y_local >> 3) * 8];
    pu_depth_intra.max = constr->ml_intra_depth_ctu->_mat_lower_depth[(x_local >> 3) + (y_local >> 3) * 8];
  }
  else {
    pu_depth_intra.min = ctrl->cfg.pu_depth_intra.min[gop_layer] >= 0 ? ctrl->cfg.pu_depth_intra.min[gop_layer] : ctrl->cfg.pu_depth_intra.min[0];
    pu_depth_intra.max = ctrl->cfg.pu_depth_intra.max[gop_layer] >= 0 ? ctrl->cfg.pu_depth_intra.max[gop_layer] : ctrl->cfg.pu_depth_intra.max[0];
  }
  pu_depth_inter.min = ctrl->cfg.pu_depth_inter.min[gop_layer] >= 0 ? ctrl->cfg.pu_depth_inter.min[gop_layer] : ctrl->cfg.pu_depth_inter.min[0];
  pu_depth_inter.max = ctrl->cfg.pu_depth_inter.max[gop_layer] >= 0 ? ctrl->cfg.pu_depth_inter.max[gop_layer] : ctrl->cfg.pu_depth_inter.max[0];

  cur_cu = LCU_GET_CU_AT_PX(lcu, x_local, y_local);
  // Assign correct depth
  cur_cu->depth = depth > MAX_DEPTH ? MAX_DEPTH : depth;
  cur_cu->tr_depth = depth > 0 ? depth : 1;
  cur_cu->type = CU_NOTSET;
  cur_cu->part_size = SIZE_2Nx2N;
  cur_cu->qp = state->qp;

  // If the CU is completely inside the frame at this depth, search for
  // prediction modes at this depth.
  if (x + cu_width <= frame->width &&
      y + cu_width <= frame->height)
  {
    int cu_width_inter_min = LCU_WIDTH >> pu_depth_inter.max;
    bool can_use_inter =
      state->frame->slicetype != KVZ_SLICE_I &&
      depth <= MAX_DEPTH &&
      (
        WITHIN(depth, pu_depth_inter.min, pu_depth_inter.max) ||
        // When the split was forced because the CTU is partially outside the
        // frame, we permit inter coding even if pu_depth_inter would
        // otherwise forbid it.
        (x & ~(cu_width_inter_min - 1)) + cu_width_inter_min > frame->width ||
        (y & ~(cu_width_inter_min - 1)) + cu_width_inter_min > frame->height
      );

    if (can_use_inter) {
      double mode_cost;
      double mode_bitcost;
      kvz_search_cu_inter(state,
                          x, y,
                          depth,
                          lcu,
                          &mode_cost, &mode_bitcost);
      if (mode_cost < cost) {
        cost = mode_cost;
        inter_bitcost = mode_bitcost;
        cur_cu->type = CU_INTER;
      }

      if (!(ctrl->cfg.early_skip && cur_cu->skipped)) {
        // Try SMP and AMP partitioning.
        static const part_mode_t mp_modes[] = {
          // SMP
          SIZE_2NxN, SIZE_Nx2N,
          // AMP
          SIZE_2NxnU, SIZE_2NxnD,
          SIZE_nLx2N, SIZE_nRx2N,
        };

        const int first_mode = ctrl->cfg.smp_enable ? 0 : 2;
        const int last_mode = (ctrl->cfg.amp_enable && cu_width >= 16) ? 5 : 1;
        for (int i = first_mode; i <= last_mode; ++i) {
          kvz_search_cu_smp(state,
                            x, y,
                            depth,
                            mp_modes[i],
                            &work_tree[depth + 1],
                            &mode_cost, &mode_bitcost);
          if (mode_cost < cost) {
            cost = mode_cost;
            inter_bitcost = mode_bitcost;
            // Copy inter prediction info to current level.
            copy_cu_info(x_local, y_local, cu_width, &work_tree[depth + 1], lcu);
          }
        }
      }
    }

    // Try to skip intra search in rd==0 mode.
    // This can be quite severe on bdrate. It might be better to do this
    // decision after reconstructing the inter frame.
    bool skip_intra = (state->encoder_control->cfg.rdo == 0
                      && cur_cu->type != CU_NOTSET
                      && cost / (cu_width * cu_width) < INTRA_THRESHOLD)
                      || (ctrl->cfg.early_skip && cur_cu->skipped);

    int32_t cu_width_intra_min = LCU_WIDTH >> pu_depth_intra.max;
    bool can_use_intra =
      (WITHIN(depth, pu_depth_intra.min, pu_depth_intra.max) ||
        // When the split was forced because the CTU is partially outside
        // the frame, we permit intra coding even if pu_depth_intra would
        // otherwise forbid it.
        (x & ~(cu_width_intra_min - 1)) + cu_width_intra_min > frame->width ||
        (y & ~(cu_width_intra_min - 1)) + cu_width_intra_min > frame->height) &&
      !(state->encoder_control->cfg.force_inter && state->frame->slicetype != KVZ_SLICE_I);

    if (can_use_intra && !skip_intra) {
      int8_t intra_mode;
      double intra_cost;
      kvz_search_cu_intra(state, x, y, depth, lcu,
                          &intra_mode, &intra_cost);
#ifdef COMPLETE_PRED_MODE_BITS
      // Technically counting these bits would be correct, however counting
      // them universally degrades quality so this block is disabled by default
      if(state->frame->slicetype != KVZ_SLICE_I) {
        double pred_mode_type_bits = 0;
        CABAC_FBITS_UPDATE(&state->search_cabac, &state->search_cabac.ctx.cu_pred_mode_model, 1, pred_mode_type_bits, "pred_mode_flag");
        CABAC_FBITS_UPDATE(&state->search_cabac, &state->search_cabac.ctx.cu_skip_flag_model[kvz_get_skip_context(x, y, lcu, NULL)], 0, pred_mode_type_bits, "skip_flag");
        intra_cost += pred_mode_type_bits * state->lambda;
      }
#endif
      if (intra_cost < cost) {
        cost = intra_cost;
        cur_cu->type = CU_INTRA;
        cur_cu->part_size = depth > MAX_DEPTH ? SIZE_NxN : SIZE_2Nx2N;
        cur_cu->intra.mode = intra_mode;
        cur_cu->skipped = 0;
        cur_cu->merged = 0;
      }
    }

    // Reconstruct best mode because we need the reconstructed pixels for
    // mode search of adjacent CUs.
    if (cur_cu->type == CU_INTRA) {
      assert(cur_cu->part_size == SIZE_2Nx2N || cur_cu->part_size == SIZE_NxN);
      cur_cu->intra.mode_chroma = cur_cu->intra.mode;
      lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_width, cur_cu);
      kvz_intra_recon_cu(state,
                         x, y,
                         depth,
                         cur_cu->intra.mode, -1, // skip chroma
                         NULL, lcu);

      if (x % 8 == 0 && y % 8 == 0 && state->encoder_control->chroma_format != KVZ_CSP_400) {
        // There is almost no benefit to doing the chroma mode search for
        // rd2. Possibly because the luma mode search already takes chroma
        // into account, so there is less of a chanse of luma mode being
        // really bad for chroma.
        if (ctrl->cfg.rdo >= 2 && ctrl->cfg.intra_chroma_search) {
          cur_cu->intra.mode_chroma = kvz_search_cu_intra_chroma(state, x, y, depth, lcu);
          lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_width, cur_cu);
        }

        kvz_intra_recon_cu(state,
                           x, y,
                           depth,
                           -1, cur_cu->intra.mode_chroma, // skip luma
                           NULL, lcu);
      }
    } else if (cur_cu->type == CU_INTER) {

      if (!cur_cu->skipped) {
        // Reset transform depth because intra messes with them.
        // This will no longer be necessary if the transform depths are not shared.
        int tr_depth = MAX(1, depth);
        if (cur_cu->part_size != SIZE_2Nx2N) {
          tr_depth = depth + 1;
        }
        kvz_lcu_fill_trdepth(lcu, x, y, depth, tr_depth);

        const bool has_chroma = state->encoder_control->chroma_format != KVZ_CSP_400;
        kvz_inter_recon_cu(state, lcu, x, y, cu_width, true, has_chroma);

        if (ctrl->cfg.zero_coeff_rdo && !ctrl->cfg.lossless && !ctrl->cfg.rdoq_enable) {
          //Calculate cost for zero coeffs
          inter_zero_coeff_cost = cu_zero_coeff_cost(state, work_tree, x, y, depth) + inter_bitcost * state->lambda;

        }

        kvz_quantize_lcu_residual(state,
          true, has_chroma,
          x, y, depth,
          NULL,
          lcu,
          false);

        int cbf = cbf_is_set_any(cur_cu->cbf, depth);

        if (cur_cu->merged && !cbf && cur_cu->part_size == SIZE_2Nx2N) {
          cur_cu->merged = 0;
          cur_cu->skipped = 1;
          // Selecting skip reduces bits needed to code the CU
          int skip_ctx = kvz_get_skip_context(x, y, lcu, NULL);
          inter_bitcost = CTX_ENTROPY_FBITS(&state->search_cabac.ctx.cu_skip_flag_model[skip_ctx], 1);
          inter_bitcost += CTX_ENTROPY_FBITS(&(state->search_cabac.ctx.cu_merge_idx_ext_model), cur_cu->merge_idx != 0);
          inter_bitcost += cur_cu->merge_idx;        
        }
      }
      lcu_fill_inter(lcu, x_local, y_local, cu_width);
      lcu_fill_cbf(lcu, x_local, y_local, cu_width, cur_cu);
    }
  }

  if (cur_cu->type == CU_INTRA || cur_cu->type == CU_INTER) {
    double bits = 0;
    cabac_data_t* cabac  = &state->search_cabac;
    cabac->update = 1;

    if(cur_cu->type != CU_INTRA || cur_cu->part_size == SIZE_2Nx2N) {
      bits += kvz_mock_encode_coding_unit(
        state,
        cabac,
        x, y, depth,
        lcu,
        cur_cu);
    }
    else {
      // Intra 4×4 PUs
      if (state->frame->slicetype != KVZ_SLICE_I) {
        cabac_ctx_t* ctx = &(cabac->ctx.cu_pred_mode_model);
        CABAC_FBITS_UPDATE(cabac, ctx, 1, bits, "pred_mode_flag");
      }
      bits += calc_mode_bits(state, lcu, cur_cu, x, y);
    }
    
    cost = bits * state->lambda;

    cost += cu_rd_cost_tr_split_accurate(state, x_local, y_local, depth, cur_cu, lcu);
    
    if (ctrl->cfg.zero_coeff_rdo && inter_zero_coeff_cost <= cost) {
      cost = inter_zero_coeff_cost;

      // Restore saved pixels from lower level of the working tree.
      copy_cu_pixels(x_local, y_local, cu_width, &work_tree[depth + 1], lcu);

      if (cur_cu->merged && cur_cu->part_size == SIZE_2Nx2N) {
        cur_cu->merged = 0;
        cur_cu->skipped = 1;
        lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_width, cur_cu);
      }

      if (cur_cu->tr_depth != depth) {
        // Reset transform depth since there are no coefficients. This
        // ensures that CBF is cleared for the whole area of the CU.
        kvz_lcu_fill_trdepth(lcu, x, y, depth, depth);
      }

      cur_cu->cbf = 0;
      lcu_fill_cbf(lcu, x_local, y_local, cu_width, cur_cu);
    }
    cabac->update = 0;
  } 

  bool can_split_cu =
    // If the CU is partially outside the frame, we need to split it even
    // if pu_depth_intra and pu_depth_inter would not permit it.
    cur_cu->type == CU_NOTSET ||
    (depth < pu_depth_intra.max && !(state->encoder_control->cfg.force_inter&& state->frame->slicetype != KVZ_SLICE_I)) ||
    (state->frame->slicetype != KVZ_SLICE_I &&
      depth < pu_depth_inter.max);

  // Recursively split all the way to max search depth.
  if (can_split_cu) {
    int half_cu = cu_width / 2;
    double split_cost = 0.0;
    int cbf = cbf_is_set_any(cur_cu->cbf, depth);
    cabac_data_t post_seach_cabac;
    memcpy(&post_seach_cabac, &state->search_cabac, sizeof(post_seach_cabac));
    memcpy(&state->search_cabac, &pre_search_cabac, sizeof(post_seach_cabac));
    state->search_cabac.update = 1;

    double split_bits = 0;

    if (depth < MAX_DEPTH) {
      // Add cost of cu_split_flag.
      uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);
      cabac_ctx_t *ctx = &(state->search_cabac.ctx.split_flag_model[split_model]);
      CABAC_FBITS_UPDATE(&state->search_cabac, ctx, 1, split_bits, "split_search");
    }

    if (cur_cu->type == CU_INTRA && depth == MAX_DEPTH) {
      // Add cost of intra part_size.
      cabac_ctx_t *ctx = &(state->search_cabac.ctx.part_size_model[0]);
      CABAC_FBITS_UPDATE(&state->search_cabac, ctx, 0, split_bits, "split_search");
    }
    state->search_cabac.update = 0;
    split_cost += split_bits * state->lambda;

    // If skip mode was selected for the block, skip further search.
    // Skip mode means there's no coefficients in the block, so splitting
    // might not give any better results but takes more time to do.
    // It is ok to interrupt the search as soon as it is known that
    // the split costs at least as much as not splitting.
    if (cur_cu->type == CU_NOTSET || cbf || state->encoder_control->cfg.cu_split_termination == KVZ_CU_SPLIT_TERMINATION_OFF) {
      if (split_cost < cost) split_cost += search_cu(state, x,           y,           depth + 1, work_tree);
      if (split_cost < cost) split_cost += search_cu(state, x + half_cu, y,           depth + 1, work_tree);
      if (split_cost < cost) split_cost += search_cu(state, x,           y + half_cu, depth + 1, work_tree);
      if (split_cost < cost) split_cost += search_cu(state, x + half_cu, y + half_cu, depth + 1, work_tree);
    } else {
      split_cost = INT_MAX;
    }

    // If no search is not performed for this depth, try just the best mode
    // of the top left CU from the next depth. This should ensure that 64x64
    // gets used, at least in the most obvious cases, while avoiding any
    // searching.
    if (cur_cu->type == CU_NOTSET && depth < MAX_PU_DEPTH
        && x + cu_width <= frame->width && y + cu_width <= frame->height 
        && state->encoder_control->cfg.combine_intra_cus)
    {

      cu_info_t *cu_d1 = LCU_GET_CU_AT_PX(&work_tree[depth + 1], x_local, y_local);

      // If the best CU in depth+1 is intra and the biggest it can be, try it.
      if (cu_d1->type == CU_INTRA && cu_d1->depth == depth + 1) {
        cabac_data_t temp_cabac;
        memcpy(&temp_cabac, &state->search_cabac, sizeof(temp_cabac));
        memcpy(&state->search_cabac, &pre_search_cabac, sizeof(pre_search_cabac));
        cost = 0;
        double bits = 0;
        if (depth < MAX_DEPTH) {
          uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);
          cabac_ctx_t* ctx = &(state->search_cabac.ctx.split_flag_model[split_model]);
          CABAC_FBITS_UPDATE(&state->search_cabac, ctx, 0, bits, "no_split_search");
        }
        else if (depth == MAX_DEPTH && cur_cu->type == CU_INTRA) {
          // Add cost of intra part_size.
          cabac_ctx_t* ctx = &(state->search_cabac.ctx.part_size_model[0]);
          CABAC_FBITS_UPDATE(&state->search_cabac, ctx, 1, bits, "no_split_search");
        }

        cur_cu->intra = cu_d1->intra;
        cur_cu->type = CU_INTRA;
        cur_cu->part_size = SIZE_2Nx2N;

        kvz_lcu_fill_trdepth(lcu, x, y, depth, cur_cu->tr_depth);
        lcu_fill_cu_info(lcu, x_local, y_local, cu_width, cu_width, cur_cu);

        const bool has_chroma = state->encoder_control->chroma_format != KVZ_CSP_400;
        const int8_t mode_chroma = has_chroma ? cur_cu->intra.mode_chroma : -1;
        kvz_intra_recon_cu(state,
                           x, y,
                           depth,
                           cur_cu->intra.mode, mode_chroma,
                           NULL, lcu);

        double mode_bits = calc_mode_bits(state, lcu, cur_cu, x, y) + bits;
        cost += mode_bits * state->lambda;

        cost += cu_rd_cost_tr_split_accurate(state, x_local, y_local, depth, cur_cu, lcu);

        memcpy(&post_seach_cabac, &state->search_cabac, sizeof(post_seach_cabac));
        memcpy(&state->search_cabac, &temp_cabac, sizeof(temp_cabac));
      }
    }

    if (split_cost < cost) {
      // Copy split modes to this depth.
      cost = split_cost;
      work_tree_copy_up(x_local, y_local, depth, work_tree);
#if KVZ_DEBUG
      debug_split = 1;
#endif
    } else if (depth > 0) {
      // Copy this CU's mode all the way down for use in adjacent CUs mode
      // search.
      memcpy(&state->search_cabac, &post_seach_cabac, sizeof(post_seach_cabac));
      work_tree_copy_down(x_local, y_local, depth, work_tree);
    }
  } else if (depth >= 0 && depth < MAX_PU_DEPTH) {
    // Need to copy modes down since the lower level of the work tree is used
    // when searching SMP and AMP blocks.
    work_tree_copy_down(x_local, y_local, depth, work_tree);
  }

  assert(cur_cu->type != CU_NOTSET);

  return cost;
}


/**
 * Initialize lcu_t for search.
 * - Copy reference CUs.
 * - Copy reference pixels from neighbouring LCUs.
 * - Copy reference pixels from this LCU.
 */
static void init_lcu_t(const encoder_state_t * const state, const int x, const int y, lcu_t *lcu, const yuv_t *hor_buf, const yuv_t *ver_buf)
{
  const videoframe_t * const frame = state->tile->frame;

  FILL(*lcu, 0);
  
  lcu->rec.chroma_format = state->encoder_control->chroma_format;
  lcu->ref.chroma_format = state->encoder_control->chroma_format;

  // Copy reference cu_info structs from neighbouring LCUs.

  // Copy top CU row.
  if (y > 0) {
    for (int i = 0; i < LCU_WIDTH; i += SCU_WIDTH) {
      const cu_info_t *from_cu = kvz_cu_array_at_const(frame->cu_array, x + i, y - 1);
      cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, i, -1);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }
  }
  // Copy left CU column.
  if (x > 0) {
    for (int i = 0; i < LCU_WIDTH; i += SCU_WIDTH) {
      const cu_info_t *from_cu = kvz_cu_array_at_const(frame->cu_array, x - 1, y + i);
      cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, -1, i);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }
  }
  // Copy top-left CU.
  if (x > 0 && y > 0) {
    const cu_info_t *from_cu = kvz_cu_array_at_const(frame->cu_array, x - 1, y - 1);
    cu_info_t *to_cu = LCU_GET_CU_AT_PX(lcu, -1, -1);
    memcpy(to_cu, from_cu, sizeof(*to_cu));
  }

  // Copy top-right CU.
  if (y > 0 && x + LCU_WIDTH < frame->width) {
    const cu_info_t *from_cu = kvz_cu_array_at_const(frame->cu_array, x + LCU_WIDTH, y - 1);
    cu_info_t *to_cu = LCU_GET_TOP_RIGHT_CU(lcu);
    memcpy(to_cu, from_cu, sizeof(*to_cu));
  }

  // Copy reference pixels.
  {
    const int pic_width = frame->width;
    // Copy top reference pixels.
    if (y > 0) {
      // hor_buf is of size pic_width so there might not be LCU_REF_PX_WIDTH
      // number of allocated pixels left.
      int x_max = MIN(LCU_REF_PX_WIDTH, pic_width - x);
      int x_min_in_lcu = (x>0) ? 0 : 1;
      int luma_offset = OFFSET_HOR_BUF(x, y, frame, x_min_in_lcu - 1);
      int chroma_offset = OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu - 1);
      int luma_bytes = (x_max + (1 - x_min_in_lcu))*sizeof(kvz_pixel);
      int chroma_bytes = (x_max / 2 + (1 - x_min_in_lcu))*sizeof(kvz_pixel);

      memcpy(&lcu->top_ref.y[x_min_in_lcu], &hor_buf->y[luma_offset], luma_bytes);
      if (state->encoder_control->chroma_format != KVZ_CSP_400) {
        memcpy(&lcu->top_ref.u[x_min_in_lcu], &hor_buf->u[chroma_offset], chroma_bytes);
        memcpy(&lcu->top_ref.v[x_min_in_lcu], &hor_buf->v[chroma_offset], chroma_bytes);
      }
    }
    // Copy left reference pixels.
    if (x > 0) {
      int y_min_in_lcu = (y>0) ? 0 : 1;
      int luma_offset = OFFSET_VER_BUF(x, y, frame, y_min_in_lcu - 1);
      int chroma_offset = OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu - 1);
      int luma_bytes = (LCU_WIDTH + (1 - y_min_in_lcu)) * sizeof(kvz_pixel);
      int chroma_bytes = (LCU_WIDTH / 2 + (1 - y_min_in_lcu)) * sizeof(kvz_pixel);

      memcpy(&lcu->left_ref.y[y_min_in_lcu], &ver_buf->y[luma_offset], luma_bytes);
      if (state->encoder_control->chroma_format != KVZ_CSP_400) {
        memcpy(&lcu->left_ref.u[y_min_in_lcu], &ver_buf->u[chroma_offset], chroma_bytes);
        memcpy(&lcu->left_ref.v[y_min_in_lcu], &ver_buf->v[chroma_offset], chroma_bytes);
      }
    }
  }

  // Copy LCU pixels.
  {
    const videoframe_t * const frame = state->tile->frame;
    int x_max = MIN(x + LCU_WIDTH, frame->width) - x;
    int y_max = MIN(y + LCU_WIDTH, frame->height) - y;

    int x_c = x / 2;
    int y_c = y / 2;
    int x_max_c = x_max / 2;
    int y_max_c = y_max / 2;

    kvz_pixels_blit(&frame->source->y[x + y * frame->source->stride], lcu->ref.y,
                        x_max, y_max, frame->source->stride, LCU_WIDTH);
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      kvz_pixels_blit(&frame->source->u[x_c + y_c * frame->source->stride / 2], lcu->ref.u,
                      x_max_c, y_max_c, frame->source->stride / 2, LCU_WIDTH / 2);
      kvz_pixels_blit(&frame->source->v[x_c + y_c * frame->source->stride / 2], lcu->ref.v,
                      x_max_c, y_max_c, frame->source->stride / 2, LCU_WIDTH / 2);
    }
  }
}


/**
 * Copy CU and pixel data to it's place in picture datastructure.
 */
static void copy_lcu_to_cu_data(const encoder_state_t * const state, int x_px, int y_px, const lcu_t *lcu)
{
  // Copy non-reference CUs to picture.
  kvz_cu_array_copy_from_lcu(state->tile->frame->cu_array, x_px, y_px, lcu);

  // Copy pixels to picture.
  {
    videoframe_t * const pic = state->tile->frame;
    const int pic_width = pic->width;
    const int x_max = MIN(x_px + LCU_WIDTH, pic_width) - x_px;
    const int y_max = MIN(y_px + LCU_WIDTH, pic->height) - y_px;

    kvz_pixels_blit(lcu->rec.y, &pic->rec->y[x_px + y_px * pic->rec->stride],
                        x_max, y_max, LCU_WIDTH, pic->rec->stride);

    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      kvz_pixels_blit(lcu->rec.u, &pic->rec->u[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
      kvz_pixels_blit(lcu->rec.v, &pic->rec->v[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    }
  }
}


/**
 * Search LCU for modes.
 * - Best mode gets copied to current picture.
 */
void kvz_search_lcu(encoder_state_t * const state, const int x, const int y, const yuv_t * const hor_buf, const yuv_t * const ver_buf)
{
  memcpy(&state->search_cabac, &state->cabac, sizeof(cabac_data_t));
  state->search_cabac.only_count = 1;
  assert(x % LCU_WIDTH == 0);
  assert(y % LCU_WIDTH == 0);

  // Initialize the same starting state to every depth. The search process
  // will use these as temporary storage for predictions before making
  // a decision on which to use, and they get updated during the search
  // process.
  lcu_t work_tree[MAX_PU_DEPTH + 1];
  init_lcu_t(state, x, y, &work_tree[0], hor_buf, ver_buf);
  for (int depth = 1; depth <= MAX_PU_DEPTH; ++depth) {
    work_tree[depth] = work_tree[0];
  }

  // If the ML depth prediction is enabled, 
  // generate the depth prediction interval 
  // for the current lcu
  constraint_t* constr = state->constraint;
  if (constr->ml_intra_depth_ctu) {
    kvz_lcu_luma_depth_pred(constr->ml_intra_depth_ctu, work_tree[0].ref.y, state->qp);
  }

  // Start search from depth 0.
  double cost = search_cu(state, x, y, 0, work_tree);

  // Save squared cost for rate control.
  if(state->encoder_control->cfg.rc_algorithm == KVZ_LAMBDA) {
    kvz_get_lcu_stats(state, x / LCU_WIDTH, y / LCU_WIDTH)->weight = cost * cost;
  }

  // The best decisions through out the LCU got propagated back to depth 0,
  // so copy those back to the frame.
  copy_lcu_to_cu_data(state, x, y, &work_tree[0]);

  // Copy coeffs to encoder state.
  copy_coeffs(work_tree[0].coeff.y, state->coeff->y, LCU_WIDTH);
  copy_coeffs(work_tree[0].coeff.u, state->coeff->u, LCU_WIDTH_C);
  copy_coeffs(work_tree[0].coeff.v, state->coeff->v, LCU_WIDTH_C);
}
