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

#include "search.h"

#include <limits.h>
#include <string.h>
#include <math.h>

#include "cabac.h"
#include "encoder.h"
#include "imagelist.h"
#include "inter.h"
#include "intra.h"
#include "kvazaar.h"
#include "rdo.h"
#include "search_inter.h"
#include "search_intra.h"
#include "strategies/strategies-picture.h"
#include "threadqueue.h"
#include "threads.h"
#include "transform.h"
#include "videoframe.h"
#include "visualization.h"


#define IN_FRAME(x, y, width, height, block_width, block_height) \
  ((x) >= 0 && (y) >= 0 \
  && (x) + (block_width) <= (width) \
  && (y) + (block_height) <= (height))

// Cost treshold for doing intra search in inter frames with --rd=0.
#ifndef INTRA_TRESHOLD
# define INTRA_TRESHOLD 20
#endif

// Disable early cu-split pruning.
#ifndef FULL_CU_SPLIT_SEARCH
#  define FULL_CU_SPLIT_SEARCH false
#endif
// Modify weight of luma SSD.
#ifndef LUMA_MULT
# define LUMA_MULT 0.8
#endif
// Modify weight of chroma SSD.
#ifndef CHROMA_MULT
# define CHROMA_MULT 1.5
#endif


/**
 * Copy all non-reference CU data from depth+1 to depth.
 */
static void work_tree_copy_up(int x_px, int y_px, int depth, lcu_t work_tree[MAX_PU_DEPTH + 1])
{
  assert(depth >= 0 && depth < MAX_PU_DEPTH);

  // Copy non-reference CUs.
  {
    const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
    const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
    const int width_cu = LCU_WIDTH >> MAX_DEPTH >> depth;
    int x, y;
    for (y = y_cu; y < y_cu + width_cu; ++y) {
      for (x = x_cu; x < x_cu + width_cu; ++x) {
        const cu_info_t *from_cu = LCU_GET_CU(&work_tree[depth + 1], x, y);
        cu_info_t *to_cu = LCU_GET_CU(&work_tree[depth], x, y);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
  }

  // Copy reconstructed pixels.
  {
    const int x = SUB_SCU(x_px);
    const int y = SUB_SCU(y_px);
    const int width_px = LCU_WIDTH >> depth;
    const int luma_index = x + y * LCU_WIDTH;
    const int chroma_index = (x / 2) + (y / 2) * (LCU_WIDTH / 2);

    const lcu_yuv_t *from = &work_tree[depth + 1].rec;
    lcu_yuv_t *to = &work_tree[depth].rec;

    const lcu_coeff_t *from_coeff = &work_tree[depth + 1].coeff;
    lcu_coeff_t *to_coeff = &work_tree[depth].coeff;

    kvz_pixels_blit(&from->y[luma_index], &to->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    kvz_pixels_blit(&from->u[chroma_index], &to->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    kvz_pixels_blit(&from->v[chroma_index], &to->v[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);

    // Copy coefficients up. They do not have to be copied down because they
    // are not used for the search.
    kvz_coefficients_blit(&from_coeff->y[luma_index], &to_coeff->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    kvz_coefficients_blit(&from_coeff->u[chroma_index], &to_coeff->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    kvz_coefficients_blit(&from_coeff->v[chroma_index], &to_coeff->v[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
  }
}


/**
 * Copy all non-reference CU data from depth to depth+1..MAX_PU_DEPTH.
 */
static void work_tree_copy_down(int x_px, int y_px, int depth, lcu_t work_tree[MAX_PU_DEPTH + 1])
{
  assert(depth >= 0 && depth < MAX_PU_DEPTH);

  // TODO: clean up to remove the copy pasta
  const int width_px = LCU_WIDTH >> depth;

  int d;

  for (d = depth + 1; d < MAX_PU_DEPTH + 1; ++d) {
    const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
    const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
    const int width_cu = width_px >> MAX_DEPTH;

    int x, y;
    for (y = y_cu; y < y_cu + width_cu; ++y) {
      for (x = x_cu; x < x_cu + width_cu; ++x) {
        const cu_info_t *from_cu = LCU_GET_CU(&work_tree[depth], x, y);
        cu_info_t *to_cu = LCU_GET_CU(&work_tree[d], x, y);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
  }

  // Copy reconstructed pixels.
  for (d = depth + 1; d < MAX_PU_DEPTH + 1; ++d) {
    const int x = SUB_SCU(x_px);
    const int y = SUB_SCU(y_px);

    const int luma_index = x + y * LCU_WIDTH;
    const int chroma_index = (x / 2) + (y / 2) * (LCU_WIDTH / 2);

    lcu_yuv_t *from = &work_tree[depth].rec;
    lcu_yuv_t *to = &work_tree[d].rec;

    kvz_pixels_blit(&from->y[luma_index], &to->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    kvz_pixels_blit(&from->u[chroma_index], &to->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    kvz_pixels_blit(&from->v[chroma_index], &to->v[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
  }
}


void kvz_lcu_set_trdepth(lcu_t *lcu, int x_px, int y_px, int depth, int tr_depth)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const vector2d_t lcu_cu = { SUB_SCU(x_px) / 8, SUB_SCU(y_px) / 8 };
  int x, y;

  // Depth 4 doesn't go inside the loop. Set the top-left CU.
  LCU_GET_CU(lcu, lcu_cu.x, lcu_cu.y)->tr_depth = tr_depth;

  for (y = 0; y < width_cu; ++y) {
    for (x = 0; x < width_cu; ++x) {
      cu_info_t *cu = LCU_GET_CU(lcu, lcu_cu.x + x, lcu_cu.y + y);
      cu->tr_depth = tr_depth;
    }
  }
}


static void lcu_set_intra_mode(lcu_t *lcu, int x_px, int y_px, int depth, int pred_mode, int chroma_mode, int part_mode)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  int x, y;

  // NxN can only be applied to a single CU at a time.
  if (part_mode == SIZE_NxN) {
    cu_info_t *cu = LCU_GET_CU(lcu, x_cu, y_cu);
    cu->depth = MAX_DEPTH;
    cu->type = CU_INTRA;
    cu->intra[PU_INDEX(x_px / 4, y_px / 4)].mode = pred_mode;
    cu->intra[PU_INDEX(x_px / 4, y_px / 4)].mode_chroma = chroma_mode;
    cu->part_size = part_mode;
    return;
  }

  // Set mode in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info_t *cu = LCU_GET_CU(lcu, x, y);
      cu->depth = depth;
      cu->type = CU_INTRA;
      cu->intra[0].mode = pred_mode;
      cu->intra[1].mode = pred_mode;
      cu->intra[2].mode = pred_mode;
      cu->intra[3].mode = pred_mode;
      cu->intra[0].mode_chroma = chroma_mode;
      cu->part_size = part_mode;
      cu->coded = 1;
    }
  }
}


static void lcu_set_inter_pu(lcu_t *lcu, int x_pu, int y_pu, int width_pu, int height_pu, cu_info_t *cur_pu)
{
  // Set mode in every CU covered by part_mode in this depth.
  for (int y = y_pu; y < y_pu + height_pu; ++y) {
    for (int x = x_pu; x < x_pu + width_pu; ++x) {
      cu_info_t *cu = LCU_GET_CU(lcu, x, y);
      //Check if this could be moved inside the if
      cu->coded    = 1;
      if (cu != cur_pu) {
        cu->depth     = cur_pu->depth;
        cu->part_size = cur_pu->part_size;
        cu->type     = CU_INTER;
        cu->tr_depth  = cur_pu->tr_depth;
        cu->merged    = cur_pu->merged;
        cu->skipped   = cur_pu->skipped;
        memcpy(&cu->inter, &cur_pu->inter, sizeof(cur_pu->inter));
      }
    }
  }
}


static void lcu_set_inter(lcu_t *lcu, int x_px, int y_px, int depth, cu_info_t *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  const int num_pu = kvz_part_mode_num_parts[cur_cu->part_size];

  for (int i = 0; i < num_pu; ++i) {
    const int x_pu      = PU_GET_X(cur_cu->part_size, width_cu, x_cu, i);
    const int y_pu      = PU_GET_Y(cur_cu->part_size, width_cu, y_cu, i);
    const int width_pu  = PU_GET_W(cur_cu->part_size, width_cu, i);
    const int height_pu = PU_GET_H(cur_cu->part_size, width_cu, i);
    cu_info_t *cur_pu   = LCU_GET_CU(lcu, x_pu, y_pu);
    lcu_set_inter_pu(lcu, x_pu, y_pu, width_pu, height_pu, cur_pu);
  }
}


static void lcu_set_coeff(lcu_t *lcu, int x_px, int y_px, int depth, cu_info_t *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  int x, y;
  int tr_split = cur_cu->tr_depth-cur_cu->depth;

  // Set coeff flags in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info_t *cu = LCU_GET_CU(lcu, x, y);
      // Use TU top-left CU to propagate coeff flags
      uint32_t mask = ~((width_cu>>tr_split)-1);
      cu_info_t *cu_from = LCU_GET_CU(lcu, x & mask, y & mask);
      if (cu != cu_from) {
        // Chroma coeff data is not used, luma is needed for deblocking
        cu->cbf.y = cu_from->cbf.y;
      }
    }
  }
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
  const uint8_t pu_index = PU_INDEX(x_px / 4, y_px / 4);

  // cur_cu is used for TU parameters.
  cu_info_t *const tr_cu = LCU_GET_CU_AT_PX(lcu, x_px, y_px);

  double coeff_bits = 0;
  double tr_tree_bits = 0;

  // Check that lcu is not in 
  assert(x_px >= 0 && x_px < LCU_WIDTH);
  assert(y_px >= 0 && y_px < LCU_WIDTH);

  const uint8_t tr_depth = tr_cu->tr_depth - depth;

  // Add transform_tree split_transform_flag bit cost.
  bool intra_split_flag = pred_cu->type == CU_INTRA && pred_cu->part_size == SIZE_NxN && depth == 3;
  if (width <= TR_MAX_WIDTH
      && width > TR_MIN_WIDTH
      && !intra_split_flag)
  {
    const cabac_ctx_t *ctx = &(state->cabac.ctx.trans_subdiv_model[5 - (6 - depth)]);
    tr_tree_bits += CTX_ENTROPY_FBITS(ctx, tr_depth > 0);
  }

  if (tr_depth > 0) {
    int offset = width / 2;
    double sum = 0;

    sum += kvz_cu_rd_cost_luma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_luma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * state->global->cur_lambda_cost;
  }

  // Add transform_tree cbf_luma bit cost.
  if (pred_cu->type == CU_INTRA ||
      tr_depth > 0 ||
      cbf_is_set(tr_cu->cbf.u, depth) ||
      cbf_is_set(tr_cu->cbf.v, depth))
  {
    const cabac_ctx_t *ctx = &(state->cabac.ctx.qt_cbf_model_luma[!tr_depth]);
    tr_tree_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.y, depth + pu_index));
  }

  unsigned ssd = 0;
  // SSD between reconstruction and original
  for (int y = y_px; y < y_px + width; ++y) {
    for (int x = x_px; x < x_px + width; ++x) {
      int diff = (int)lcu->rec.y[y * LCU_WIDTH + x] - (int)lcu->ref.y[y * LCU_WIDTH + x];
      ssd += diff*diff;
    }
  }

  {
    coeff_t coeff_temp[32 * 32];
    int8_t luma_scan_mode = kvz_get_scan_order(pred_cu->type, pred_cu->intra[PU_INDEX(x_px / 4, y_px / 4)].mode, depth);

    // Code coeffs using cabac to get a better estimate of real coding costs.
    kvz_coefficients_blit(&lcu->coeff.y[(y_px*LCU_WIDTH) + x_px], coeff_temp, width, width, LCU_WIDTH, width);
    coeff_bits += kvz_get_coeff_cost(state, coeff_temp, width, 0, luma_scan_mode);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * LUMA_MULT + bits * state->global->cur_lambda_cost;
}


double kvz_cu_rd_cost_chroma(const encoder_state_t *const state,
                         const int x_px, const int y_px, const int depth,
                         const cu_info_t *const pred_cu,
                         lcu_t *const lcu)
{
  const vector2d_t lcu_px = { x_px / 2, y_px / 2 };
  const int width = (depth <= MAX_DEPTH) ? LCU_WIDTH >> (depth + 1) : LCU_WIDTH >> depth;
  cu_info_t *const tr_cu = LCU_GET_CU(lcu, lcu_px.x / 4, lcu_px.y / 4);

  double tr_tree_bits = 0;
  double coeff_bits = 0;

  assert(x_px >= 0 && x_px < LCU_WIDTH);
  assert(y_px >= 0 && y_px < LCU_WIDTH);

  if (PU_INDEX(x_px / 4, y_px / 4) != 0) {
    // For MAX_PU_DEPTH calculate chroma for previous depth for the first
    // block and return 0 cost for all others.
    return 0;
  }

  if (depth < MAX_PU_DEPTH) {
    const int tr_depth = depth - pred_cu->depth;
    const cabac_ctx_t *ctx = &(state->cabac.ctx.qt_cbf_model_chroma[tr_depth]);
    if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.u, depth - 1)) {
      tr_tree_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.u, depth));
    }
    if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.v, depth - 1)) {
      tr_tree_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.v, depth));
    }
  }

  if (tr_cu->tr_depth > depth) {
    int offset = LCU_WIDTH >> (depth + 1);
    int sum = 0;

    sum += kvz_cu_rd_cost_chroma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += kvz_cu_rd_cost_chroma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * state->global->cur_lambda_cost;
  }

  // Chroma SSD
  int ssd = 0;
  for (int y = lcu_px.y; y < lcu_px.y + width; ++y) {
    for (int x = lcu_px.x; x < lcu_px.x + width; ++x) {
      int diff = (int)lcu->rec.u[y * LCU_WIDTH_C + x] - (int)lcu->ref.u[y * LCU_WIDTH_C + x];
      ssd += diff * diff;
    }
  }
  for (int y = lcu_px.y; y < lcu_px.y + width; ++y) {
    for (int x = lcu_px.x; x < lcu_px.x + width; ++x) {
      int diff = (int)lcu->rec.v[y * LCU_WIDTH_C + x] - (int)lcu->ref.v[y * LCU_WIDTH_C + x];
      ssd += diff * diff;
    }
  }

  {
    coeff_t coeff_temp[16 * 16];
    int8_t scan_order = kvz_get_scan_order(pred_cu->type, pred_cu->intra[0].mode_chroma, depth);
    
    kvz_coefficients_blit(&lcu->coeff.u[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += kvz_get_coeff_cost(state, coeff_temp, width, 2, scan_order);

    kvz_coefficients_blit(&lcu->coeff.v[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += kvz_get_coeff_cost(state, coeff_temp, width, 2, scan_order);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * CHROMA_MULT + bits * state->global->cur_lambda_cost;
}


// Return estimate of bits used to code prediction mode of cur_cu.
static double calc_mode_bits(const encoder_state_t *state,
                             const cu_info_t * cur_cu,
                             int x, int y)
{
  double mode_bits;
  
  if (cur_cu->type == CU_INTER) {
    mode_bits = cur_cu->inter.bitcost;
  } else {
    int8_t candidate_modes[3];
    {
      const cu_info_t *left_cu  = ((x >= 8) ? CU_GET_CU(cur_cu, -1,  0) : NULL);
      const cu_info_t *above_cu = ((y >= 8) ? CU_GET_CU(cur_cu,  0, -1) : NULL);
      kvz_intra_get_dir_luma_predictor(x, y, candidate_modes, cur_cu, left_cu, above_cu);
    }

    mode_bits = kvz_luma_mode_bits(state, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode, candidate_modes);
    if (PU_INDEX(x >> 2, y >> 2) == 0) {
      mode_bits += kvz_chroma_mode_bits(state, cur_cu->intra[0].mode_chroma, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode);
    }
  }

  return mode_bits;
}


static uint8_t get_ctx_cu_split_model(const lcu_t *lcu, int x, int y, int depth)
{
  vector2d_t lcu_cu = { SUB_SCU(x) / 8, SUB_SCU(y) / 8 };
  bool condA = x >= 8 && LCU_GET_CU(lcu, lcu_cu.x - 1, lcu_cu.y    )->depth > depth;
  bool condL = y >= 8 && LCU_GET_CU(lcu, lcu_cu.x,     lcu_cu.y - 1)->depth > depth;
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
static double search_cu(encoder_state_t * const state, int x, int y, int depth, lcu_t work_tree[MAX_PU_DEPTH + 1])
{
  const encoder_control_t* ctrl = state->encoder_control;
  const videoframe_t * const frame = state->tile->frame;
  int cu_width = LCU_WIDTH >> depth;
  double cost = MAX_INT;
  cu_info_t *cur_cu;

  lcu_t *const lcu = &work_tree[depth];

  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);
#ifdef KVZ_DEBUG
  int debug_split = 0;
#endif
  PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHCU);

  // Stop recursion if the CU is completely outside the frame.
  if (x >= frame->width || y >= frame->height) {
    // Return zero cost because this CU does not have to be coded.
    return 0;
  }

  cur_cu = LCU_GET_CU_AT_PX(&work_tree[depth], x_local, y_local);
  // Assign correct depth
  cur_cu->depth = depth > MAX_DEPTH ? MAX_DEPTH : depth;
  cur_cu->tr_depth = depth > 0 ? depth : 1;
  cur_cu->type = CU_NOTSET;
  cur_cu->part_size = SIZE_2Nx2N;
  // If the CU is completely inside the frame at this depth, search for
  // prediction modes at this depth.
  if (x + cu_width <= frame->width &&
      y + cu_width <= frame->height)
  {

    bool can_use_inter =
        state->global->slicetype != KVZ_SLICE_I
        && WITHIN(depth, ctrl->pu_depth_inter.min, ctrl->pu_depth_inter.max);

    if (can_use_inter) {
      double mode_cost = kvz_search_cu_inter(state, x, y, depth, &work_tree[depth]);
      if (mode_cost < cost) {
        cost = mode_cost;
        cur_cu->type = CU_INTER;
      }

      if (depth < MAX_DEPTH) {
        // Try SMP and AMP partitioning.
        static const part_mode_t mp_modes[] = {
          // SMP
          SIZE_2NxN, SIZE_Nx2N,
          // AMP
          SIZE_2NxnU, SIZE_2NxnD,
          SIZE_nLx2N, SIZE_nRx2N,
        };

        const int first_mode = ctrl->cfg->smp_enable ? 0 : 2;
        const int last_mode  = (ctrl->cfg->amp_enable && cu_width >= 32) ? 5 : 1;
        for (int i = first_mode; i <= last_mode; ++i) {
          mode_cost = kvz_search_cu_smp(state,
                                        x, y,
                                        depth,
                                        mp_modes[i],
                                        &work_tree[depth + 1]);
          // TODO: take cost of coding part mode into account
          if (mode_cost < cost) {
            cost = mode_cost;
            // TODO: only copy inter prediction info, not pixels
            work_tree_copy_up(x, y, depth, work_tree);
    }
        }
      }
    }

    // Try to skip intra search in rd==0 mode.
    // This can be quite severe on bdrate. It might be better to do this
    // decision after reconstructing the inter frame.
    bool skip_intra = state->encoder_control->rdo == 0
                      && cur_cu->type != CU_NOTSET
                      && cost / (cu_width * cu_width) < INTRA_TRESHOLD;
    if (!skip_intra 
        && WITHIN(depth, ctrl->pu_depth_intra.min, ctrl->pu_depth_intra.max))
    {
      double mode_cost = kvz_search_cu_intra(state, x, y, depth, &work_tree[depth]);
      if (mode_cost < cost) {
        cost = mode_cost;
        cur_cu->type = CU_INTRA;
        cur_cu->part_size = depth > MAX_DEPTH ? SIZE_NxN : SIZE_2Nx2N;
      }
    }

    // Reconstruct best mode because we need the reconstructed pixels for
    // mode search of adjacent CUs.
    if (cur_cu->type == CU_INTRA) {
      assert(cur_cu->part_size == SIZE_2Nx2N || cur_cu->part_size == SIZE_NxN);
      int8_t intra_mode = cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode;
      lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                         intra_mode,
                         intra_mode,
                         cur_cu->part_size);
      kvz_intra_recon_lcu_luma(state, x, y, depth, intra_mode, NULL, &work_tree[depth]);

      if (PU_INDEX(x >> 2, y >> 2) == 0) {
        int8_t intra_mode_chroma = intra_mode;

        // There is almost no benefit to doing the chroma mode search for
        // rd2. Possibly because the luma mode search already takes chroma
        // into account, so there is less of a chanse of luma mode being
        // really bad for chroma.
        if (state->encoder_control->rdo == 3) {
          intra_mode_chroma = kvz_search_cu_intra_chroma(state, x, y, depth, &work_tree[depth]);
          lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                             intra_mode, intra_mode_chroma,
                             cur_cu->part_size);
        }

        kvz_intra_recon_lcu_chroma(state, x, y, depth, intra_mode_chroma, NULL, &work_tree[depth]);
      }
    } else if (cur_cu->type == CU_INTER) {
      // Reset transform depth because intra messes with them.
      // This will no longer be necessary if the transform depths are not shared.
      int tr_depth = depth > 0 ? depth : 1;
      kvz_lcu_set_trdepth(&work_tree[depth], x, y, depth, tr_depth);

      const int cu_width = LCU_WIDTH >> depth;
      const int num_pu = kvz_part_mode_num_parts[cur_cu->part_size];

      for (int i = 0; i < num_pu; ++i) {
        const int pu_x = PU_GET_X(cur_cu->part_size, cu_width, x, i);
        const int pu_y = PU_GET_Y(cur_cu->part_size, cu_width, y, i);
        const int pu_w = PU_GET_W(cur_cu->part_size, cu_width, i);
        const int pu_h = PU_GET_H(cur_cu->part_size, cu_width, i);

        cu_info_t *cur_pu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(pu_x), SUB_SCU(pu_y));

        if (cur_pu->inter.mv_dir == 3) {
          const kvz_picture *const refs[2] = {
            state->global->ref->images[cur_pu->inter.mv_ref[0]],
            state->global->ref->images[cur_pu->inter.mv_ref[1]],
          };
          kvz_inter_recon_lcu_bipred(state,
                                     refs[0], refs[1],
                                     pu_x, pu_y,
                                     pu_w, pu_h,
                                     cur_pu->inter.mv,
                                     &work_tree[depth]);
      } else {
          const int mv_idx = cur_pu->inter.mv_dir - 1;
          const kvz_picture *const ref =
              state->global->ref->images[cur_pu->inter.mv_ref[mv_idx]];
          kvz_inter_recon_lcu(state,
                              ref,
                              pu_x, pu_y,
                              pu_w, pu_h,
                              cur_pu->inter.mv[mv_idx],
                              &work_tree[depth],
                              0);
      }
      }

      kvz_quantize_lcu_luma_residual(state, x, y, depth, NULL, &work_tree[depth]);
      kvz_quantize_lcu_chroma_residual(state, x, y, depth, NULL, &work_tree[depth]);

      int cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);

      if(cur_cu->merged && !cbf && cur_cu->part_size == SIZE_2Nx2N) {
        cur_cu->merged = 0;
        cur_cu->skipped = 1;
        // Selecting skip reduces bits needed to code the CU
        if (cur_cu->inter.bitcost > 1) {
          cur_cu->inter.bitcost -= 1;
        }
      }
      lcu_set_inter(&work_tree[depth], x, y, depth, cur_cu);
      lcu_set_coeff(&work_tree[depth], x, y, depth, cur_cu);
    }
  }
  if (cur_cu->type == CU_INTRA || cur_cu->type == CU_INTER) {
    cost = kvz_cu_rd_cost_luma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    cost += kvz_cu_rd_cost_chroma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    double mode_bits = calc_mode_bits(state, cur_cu, x, y);
    cost += mode_bits * state->global->cur_lambda_cost;
  }
  
  // Recursively split all the way to max search depth.
  if (depth < ctrl->pu_depth_intra.max || (depth < ctrl->pu_depth_inter.max && state->global->slicetype != KVZ_SLICE_I)) {
    int half_cu = cu_width / 2;
    double split_cost = 0.0;
    int cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);
        
    if (depth < MAX_DEPTH) {
      // Add cost of cu_split_flag.
      uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);
      const cabac_ctx_t *ctx = &(state->cabac.ctx.split_flag_model[split_model]);
      cost += CTX_ENTROPY_FBITS(ctx, 0) * state->global->cur_lambda_cost;
      split_cost += CTX_ENTROPY_FBITS(ctx, 1) * state->global->cur_lambda_cost;
    }

    if (cur_cu->type == CU_INTRA && depth == MAX_DEPTH) {
      // Add cost of intra part_size.
      const cabac_ctx_t *ctx = &(state->cabac.ctx.part_size_model[0]);
      cost += CTX_ENTROPY_FBITS(ctx, 1) * state->global->cur_lambda_cost;  // 2Nx2N
      split_cost += CTX_ENTROPY_FBITS(ctx, 0) * state->global->cur_lambda_cost;  // NxN
    }

    // If skip mode was selected for the block, skip further search.
    // Skip mode means there's no coefficients in the block, so splitting
    // might not give any better results but takes more time to do.
    if (cur_cu->type == CU_NOTSET || cbf || FULL_CU_SPLIT_SEARCH) {
      split_cost += search_cu(state, x,           y,           depth + 1, work_tree);
      split_cost += search_cu(state, x + half_cu, y,           depth + 1, work_tree);
      split_cost += search_cu(state, x,           y + half_cu, depth + 1, work_tree);
      split_cost += search_cu(state, x + half_cu, y + half_cu, depth + 1, work_tree);
    } else {
      split_cost = INT_MAX;
    }

    // If no search is not performed for this depth, try just the best mode
    // of the top left CU from the next depth. This should ensure that 64x64
    // gets used, at least in the most obvious cases, while avoiding any
    // searching.
    if (cur_cu->type == CU_NOTSET && depth < MAX_PU_DEPTH
        && x + cu_width <= frame->width && y + cu_width <= frame->height)
    {
      vector2d_t lcu_cu = { x_local / 8, y_local / 8 };
      cu_info_t *cu_d1 = LCU_GET_CU(&work_tree[depth + 1], lcu_cu.x, lcu_cu.y);

      // If the best CU in depth+1 is intra and the biggest it can be, try it.
      if (cu_d1->type == CU_INTRA && cu_d1->depth == depth + 1) {
        cost = 0;

        cur_cu->intra[0] = cu_d1->intra[0];
        cur_cu->type = CU_INTRA;
        cur_cu->part_size = depth > MAX_DEPTH ? SIZE_NxN : SIZE_2Nx2N;

        kvz_lcu_set_trdepth(&work_tree[depth], x, y, depth, cur_cu->tr_depth);
        lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                           cur_cu->intra[0].mode, cur_cu->intra[0].mode_chroma,
                           cur_cu->part_size);
        kvz_intra_recon_lcu_luma(state, x, y, depth, cur_cu->intra[0].mode, NULL, &work_tree[depth]);
        kvz_intra_recon_lcu_chroma(state, x, y, depth, cur_cu->intra[0].mode_chroma, NULL, &work_tree[depth]);
        cost += kvz_cu_rd_cost_luma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
        cost += kvz_cu_rd_cost_chroma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);

        // Add the cost of coding no-split.
        uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);
        const cabac_ctx_t *ctx = &(state->cabac.ctx.split_flag_model[split_model]);
        cost += CTX_ENTROPY_FBITS(ctx, 0) * state->global->cur_lambda_cost;

        // Add the cost of coding intra mode only once.
        double mode_bits = calc_mode_bits(state, cur_cu, x, y);
        cost += mode_bits * state->global->cur_lambda_cost;
      }
    }

    if (split_cost < cost) {
      // Copy split modes to this depth.
      cost = split_cost;
      work_tree_copy_up(x, y, depth, work_tree);
#if KVZ_DEBUG
      debug_split = 1;
#endif
    } else if (depth > 0) {
      // Copy this CU's mode all the way down for use in adjacent CUs mode
      // search.
      work_tree_copy_down(x, y, depth, work_tree);
    }
  } else if (depth >= 0 && depth < MAX_PU_DEPTH) {
    // Need to copy modes down since the lower level of the work tree is used
    // when searching SMP and AMP blocks.
    work_tree_copy_down(x, y, depth, work_tree);
  }

#if KVZ_VISUALIZATION == 1
  kvz_visualization_draw_block_with_delay(state, &work_tree[depth], cur_cu, x, y, depth);
#endif
  
  PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHCU, state->encoder_control->threadqueue, "type=search_cu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d,depth=%d,split=%d,cur_cu_is_intra=%d", state->global->frame, state->tile->id, state->slice->id,
                          (state->tile->lcu_offset_x * LCU_WIDTH) + x,
                          (state->tile->lcu_offset_x * LCU_WIDTH) + x + (LCU_WIDTH >> depth), 
                          (state->tile->lcu_offset_y * LCU_WIDTH) + y,
                          (state->tile->lcu_offset_y * LCU_WIDTH) + y + (LCU_WIDTH >> depth), 
                          depth, debug_split, (cur_cu->type==CU_INTRA)?1:0);

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
  
  // Copy reference cu_info structs from neighbouring LCUs.
  {
    const int x_cu = x >> MAX_DEPTH;
    const int y_cu = y >> MAX_DEPTH;

    // Copy top CU row.
    if (y_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info_t *from_cu = kvz_videoframe_get_cu_const(frame, x_cu + i, y_cu - 1);
        cu_info_t *to_cu = LCU_GET_CU(lcu, i, -1);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy left CU column.
    if (x_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info_t *from_cu = kvz_videoframe_get_cu_const(frame, x_cu - 1, y_cu + i);
        cu_info_t *to_cu = LCU_GET_CU(lcu, -1, i);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy top-left CU.
    if (x_cu > 0 && y_cu > 0) {
      const cu_info_t *from_cu = kvz_videoframe_get_cu_const(frame, x_cu - 1, y_cu - 1);
      cu_info_t *to_cu = LCU_GET_CU(lcu, -1, -1);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }

    // Copy top-right CU.
    if (y_cu > 0 && x + LCU_WIDTH < frame->width) {
      const cu_info_t *from_cu = kvz_videoframe_get_cu_const(frame, x_cu + LCU_CU_WIDTH, y_cu - 1);
      cu_info_t *to_cu = LCU_GET_TOP_RIGHT_CU(lcu);
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }
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
      memcpy(&lcu->top_ref.y[x_min_in_lcu], &hor_buf->y[OFFSET_HOR_BUF(x, y, frame, x_min_in_lcu-1)], (x_max + (1-x_min_in_lcu))*sizeof(kvz_pixel));
      memcpy(&lcu->top_ref.u[x_min_in_lcu], &hor_buf->u[OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu - 1)], (x_max / 2 + (1 - x_min_in_lcu))*sizeof(kvz_pixel));
      memcpy(&lcu->top_ref.v[x_min_in_lcu], &hor_buf->v[OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu - 1)], (x_max / 2 + (1 - x_min_in_lcu))*sizeof(kvz_pixel));
    }
    // Copy left reference pixels.
    if (x > 0) {
      int y_min_in_lcu = (y>0) ? 0 : 1;
      memcpy(&lcu->left_ref.y[y_min_in_lcu], &ver_buf->y[OFFSET_VER_BUF(x, y, frame, y_min_in_lcu - 1)], (LCU_WIDTH + (1 - y_min_in_lcu))*sizeof(kvz_pixel));
      memcpy(&lcu->left_ref.u[y_min_in_lcu], &ver_buf->u[OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu - 1)], (LCU_WIDTH / 2 + (1 - y_min_in_lcu))*sizeof(kvz_pixel));
      memcpy(&lcu->left_ref.v[y_min_in_lcu], &ver_buf->v[OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu - 1)], (LCU_WIDTH / 2 + (1 - y_min_in_lcu))*sizeof(kvz_pixel));
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
    kvz_pixels_blit(&frame->source->u[x_c + y_c * frame->source->stride/2], lcu->ref.u,
                        x_max_c, y_max_c, frame->source->stride/2, LCU_WIDTH / 2);
    kvz_pixels_blit(&frame->source->v[x_c + y_c * frame->source->stride/2], lcu->ref.v,
                        x_max_c, y_max_c, frame->source->stride/2, LCU_WIDTH / 2);
  }
}


/**
 * Copy CU and pixel data to it's place in picture datastructure.
 */
static void copy_lcu_to_cu_data(const encoder_state_t * const state, int x_px, int y_px, const lcu_t *lcu)
{
  // Copy non-reference CUs to picture.
  {
    const int x_cu = x_px >> MAX_DEPTH;
    const int y_cu = y_px >> MAX_DEPTH;
    videoframe_t * const frame = state->tile->frame;

    int x, y;
    for (y = 0; y < LCU_CU_WIDTH; ++y) {
      for (x = 0; x < LCU_CU_WIDTH; ++x) {
        const cu_info_t *from_cu = LCU_GET_CU(lcu, x, y);
        cu_info_t *to_cu = kvz_videoframe_get_cu(frame, x_cu + x, y_cu + y);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
  }

  // Copy pixels to picture.
  {
    videoframe_t * const pic = state->tile->frame;
    const int pic_width = pic->width;
    const int x_max = MIN(x_px + LCU_WIDTH, pic_width) - x_px;
    const int y_max = MIN(y_px + LCU_WIDTH, pic->height) - y_px;
    const int luma_index = x_px + y_px * pic_width;
    const int chroma_index = (x_px / 2) + (y_px / 2) * (pic_width / 2);

    kvz_pixels_blit(lcu->rec.y, &pic->rec->y[x_px + y_px * pic->rec->stride],
                        x_max, y_max, LCU_WIDTH, pic->rec->stride);
    kvz_coefficients_blit(lcu->coeff.y, &pic->coeff_y[luma_index],
                        x_max, y_max, LCU_WIDTH, pic_width);

    kvz_pixels_blit(lcu->rec.u, &pic->rec->u[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    kvz_pixels_blit(lcu->rec.v, &pic->rec->v[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    kvz_coefficients_blit(lcu->coeff.u, &pic->coeff_u[chroma_index],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);
    kvz_coefficients_blit(lcu->coeff.v, &pic->coeff_v[chroma_index],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);
  }
}


/**
 * Search LCU for modes.
 * - Best mode gets copied to current picture.
 */
void kvz_search_lcu(encoder_state_t * const state, const int x, const int y, const yuv_t * const hor_buf, const yuv_t * const ver_buf)
{
  // Initialize the same starting state to every depth. The search process
  // will use these as temporary storage for predictions before making
  // a decision on which to use, and they get updated during the search
  // process.
  lcu_t work_tree[MAX_PU_DEPTH + 1];
  init_lcu_t(state, x, y, &work_tree[0], hor_buf, ver_buf);
  for (int depth = 1; depth <= MAX_PU_DEPTH; ++depth) {
    work_tree[depth] = work_tree[0];
  }

#if KVZ_VISUALIZATION == 1
  kvz_visualization_mv_clear_lcu(state, x, y);
#endif

  // Start search from depth 0.
  search_cu(state, x, y, 0, work_tree);

#if KVZ_VISUALIZATION == 1
  kvz_visualization_mv_draw_lcu(state, x, y, &work_tree[0]);
#endif

  // The best decisions through out the LCU got propagated back to depth 0,
  // so copy those back to the frame.
  copy_lcu_to_cu_data(state, x, y, &work_tree[0]);
}
