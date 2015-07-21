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

#include "search.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "bitstream.h"
#include "image.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-ipol.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "rdo.h"
#include "transform.h"
#include "encoder.h"
#include "search_inter.h"

#define IN_FRAME(x, y, width, height, block_width, block_height) \
  ((x) >= 0 && (y) >= 0 \
  && (x) + (block_width) <= (width) \
  && (y) + (block_height) <= (height))

// Cost treshold for doing intra search in inter frames with --rd=0.
#ifndef INTRA_TRESHOLD
# define INTRA_TRESHOLD 20
#endif

// Extra cost for CU split.
// Compensates for missing or incorrect bit costs. Must be recalculated if
// bits are added or removed from cu-tree search.
#ifndef CU_COST
#  define CU_COST 3
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
// Normalize SAD for comparison against SATD to estimate transform skip
// for 4x4 blocks.
#ifndef TRSKIP_RATIO
# define TRSKIP_RATIO 1.7
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
        const cu_info_t *from_cu = &work_tree[depth + 1].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
        cu_info_t *to_cu = &work_tree[depth].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
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

    pixels_blit(&from->y[luma_index], &to->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    pixels_blit(&from->u[chroma_index], &to->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    pixels_blit(&from->v[chroma_index], &to->v[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);

    // Copy coefficients up. They do not have to be copied down because they
    // are not used for the search.
    coefficients_blit(&from_coeff->y[luma_index], &to_coeff->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    coefficients_blit(&from_coeff->u[chroma_index], &to_coeff->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    coefficients_blit(&from_coeff->v[chroma_index], &to_coeff->v[chroma_index],
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
        const cu_info_t *from_cu = &work_tree[depth].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
        cu_info_t *to_cu = &work_tree[d].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
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

    pixels_blit(&from->y[luma_index], &to->y[luma_index],
                        width_px, width_px, LCU_WIDTH, LCU_WIDTH);
    pixels_blit(&from->u[chroma_index], &to->u[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
    pixels_blit(&from->v[chroma_index], &to->v[chroma_index],
                        width_px / 2, width_px / 2, LCU_WIDTH / 2, LCU_WIDTH / 2);
  }
}


static void lcu_set_trdepth(lcu_t *lcu, int x_px, int y_px, int depth, int tr_depth)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const vector2d_t lcu_cu = { (x_px & (LCU_WIDTH - 1)) / 8, (y_px & (LCU_WIDTH - 1)) / 8 };
  cu_info_t *const cur_cu = &lcu->cu[lcu_cu.x + lcu_cu.y * LCU_T_CU_WIDTH + LCU_CU_OFFSET];
  int x, y;

  // Depth 4 doesn't go inside the loop. Set the top-left CU.
  cur_cu->tr_depth = tr_depth;

  for (y = 0; y < width_cu; ++y) {
    for (x = 0; x < width_cu; ++x) {
      cu_info_t *cu = &cur_cu[x + y * LCU_T_CU_WIDTH];
      cu->tr_depth = tr_depth;
    }
  }
}


static void lcu_set_intra_mode(lcu_t *lcu, int x_px, int y_px, int depth, int pred_mode, int chroma_mode, int part_mode)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info_t *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;

  // NxN can only be applied to a single CU at a time.
  if (part_mode == SIZE_NxN) {
    cu_info_t *cu = &lcu_cu[x_cu + y_cu * LCU_T_CU_WIDTH];
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
      cu_info_t *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
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


static void lcu_set_inter(lcu_t *lcu, int x_px, int y_px, int depth, cu_info_t *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info_t *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;
  // Set mode in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info_t *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
      //Check if this could be moved inside the if
      cu->coded    = 1;
      if (cu != cur_cu) {
        cu->depth    = cur_cu->depth;
        cu->type     = CU_INTER;
        cu->tr_depth = cur_cu->tr_depth;
        cu->merged   = cur_cu->merged;
        cu->skipped  = cur_cu->skipped;
        memcpy(&cu->inter, &cur_cu->inter, sizeof(cur_cu->inter));
      }
    }
  }
}


static void lcu_set_coeff(lcu_t *lcu, int x_px, int y_px, int depth, cu_info_t *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info_t *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;
  int tr_split = cur_cu->tr_depth-cur_cu->depth;

  // Set coeff flags in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info_t *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
      // Use TU top-left CU to propagate coeff flags
      uint32_t mask = ~((width_cu>>tr_split)-1);
      cu_info_t *cu_from = &lcu_cu[(x & mask) + (y & mask) * LCU_T_CU_WIDTH];
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
static double cu_rd_cost_luma(const encoder_state_t *const state,
                              const int x_px, const int y_px, const int depth,
                              const cu_info_t *const pred_cu,
                              lcu_t *const lcu)
{
  const int width = LCU_WIDTH >> depth;
  const uint8_t pu_index = PU_INDEX(x_px / 4, y_px / 4);

  // cur_cu is used for TU parameters.
  cu_info_t *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (x_px / 8) + (y_px / 8) * LCU_T_CU_WIDTH];

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

    sum += cu_rd_cost_luma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

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
    int8_t luma_scan_mode = get_scan_order(pred_cu->type, pred_cu->intra[PU_INDEX(x_px / 4, y_px / 4)].mode, depth);

    // Code coeffs using cabac to get a better estimate of real coding costs.
    coefficients_blit(&lcu->coeff.y[(y_px*LCU_WIDTH) + x_px], coeff_temp, width, width, LCU_WIDTH, width);
    coeff_bits += get_coeff_cost(state, coeff_temp, width, 0, luma_scan_mode);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * LUMA_MULT + bits * state->global->cur_lambda_cost;
}


static double cu_rd_cost_chroma(const encoder_state_t *const state,
                                const int x_px, const int y_px, const int depth,
                                const cu_info_t *const pred_cu,
                                lcu_t *const lcu)
{
  const vector2d_t lcu_px = { x_px / 2, y_px / 2 };
  const int width = (depth <= MAX_DEPTH) ? LCU_WIDTH >> (depth + 1) : LCU_WIDTH >> depth;
  cu_info_t *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x / 4) + (lcu_px.y / 4)*LCU_T_CU_WIDTH];

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

    sum += cu_rd_cost_chroma(state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

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
    int8_t scan_order = get_scan_order(pred_cu->type, pred_cu->intra[0].mode_chroma, depth);
    
    coefficients_blit(&lcu->coeff.u[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += get_coeff_cost(state, coeff_temp, width, 2, scan_order);

    coefficients_blit(&lcu->coeff.v[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += get_coeff_cost(state, coeff_temp, width, 2, scan_order);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * CHROMA_MULT + bits * state->global->cur_lambda_cost;
}


/**
* \brief Perform search for best intra transform split configuration.
*
* This function does a recursive search for the best intra transform split
* configuration for a given intra prediction mode.
*
* \return RD cost of best transform split configuration. Splits in lcu->cu.
* \param depth  Current transform depth.
* \param max_depth  Depth to which TR split will be tried.
* \param intra_mode  Intra prediction mode.
* \param cost_treshold  RD cost at which search can be stopped.
*/
static double search_intra_trdepth(encoder_state_t * const state,
                                   int x_px, int y_px, int depth, int max_depth,
                                   int intra_mode, int cost_treshold,
                                   cu_info_t *const pred_cu,
                                   lcu_t *const lcu)
{
  assert(depth >= 0 && depth <= MAX_PU_DEPTH);

  const int width = LCU_WIDTH >> depth;
  const int width_c = width > TR_MIN_WIDTH ? width / 2 : width;

  const int offset = width / 2;
  const vector2d_t lcu_px = { x_px & 0x3f, y_px & 0x3f };
  cu_info_t *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];

  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);

  struct {
    kvz_pixel y[TR_MAX_WIDTH*TR_MAX_WIDTH];
    kvz_pixel u[TR_MAX_WIDTH*TR_MAX_WIDTH];
    kvz_pixel v[TR_MAX_WIDTH*TR_MAX_WIDTH];
  } nosplit_pixels;
  cu_cbf_t nosplit_cbf;

  double split_cost = INT32_MAX;
  double nosplit_cost = INT32_MAX;

  if (depth > 0) {
    tr_cu->tr_depth = depth;
    pred_cu->tr_depth = depth;

    nosplit_cost = 0.0;

    cbf_clear(&pred_cu->cbf.y, depth + PU_INDEX(x_px / 4, y_px / 4));

    intra_recon_lcu_luma(state, x_px, y_px, depth, intra_mode, pred_cu, lcu);
    nosplit_cost += cu_rd_cost_luma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);

    if (reconstruct_chroma) {
      cbf_clear(&pred_cu->cbf.u, depth);
      cbf_clear(&pred_cu->cbf.v, depth);

      intra_recon_lcu_chroma(state, x_px, y_px, depth, intra_mode, pred_cu, lcu);
      nosplit_cost += cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
    }

    // Early stop codition for the recursive search.
    // If the cost of any 1/4th of the transform is already larger than the
    // whole transform, assume that splitting further is a bad idea.
    if (nosplit_cost >= cost_treshold) {
      return nosplit_cost;
    }

    nosplit_cbf = pred_cu->cbf;

    pixels_blit(lcu->rec.y, nosplit_pixels.y, width, width, LCU_WIDTH, width);
    if (reconstruct_chroma) {
      pixels_blit(lcu->rec.u, nosplit_pixels.u, width_c, width_c, LCU_WIDTH_C, width_c);
      pixels_blit(lcu->rec.v, nosplit_pixels.v, width_c, width_c, LCU_WIDTH_C, width_c);
    }
  }

  // Recurse further if all of the following:
  // - Current depth is less than maximum depth of the search (max_depth).
  //   - Maximum transform hierarchy depth is constrained by clipping
  //     max_depth.
  // - Min transform size hasn't been reached (MAX_PU_DEPTH).
  if (depth < max_depth && depth < MAX_PU_DEPTH) {
    split_cost = 3 * state->global->cur_lambda_cost;

    split_cost += search_intra_trdepth(state, x_px, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(state, x_px + offset, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }

    double tr_split_bit = 0.0;
    double cbf_bits = 0.0;

    // Add bits for split_transform_flag = 1, because transform depth search bypasses
    // the normal recursion in the cost functions.
    if (depth >= 1 && depth <= 3) {
      const cabac_ctx_t *ctx = &(state->cabac.ctx.trans_subdiv_model[5 - (6 - depth)]);
      tr_split_bit += CTX_ENTROPY_FBITS(ctx, 1);
    }

    // Add cost of cbf chroma bits on transform tree.
    // All cbf bits are accumulated to pred_cu.cbf and cbf_is_set returns true
    // if cbf is set at any level >= depth, so cbf chroma is assumed to be 0
    // if this and any previous transform block has no chroma coefficients.
    // When searching the first block we don't actually know the real values,
    // so this will code cbf as 0 and not code the cbf at all for descendants.
    {
      const uint8_t tr_depth = depth - pred_cu->depth;

      const cabac_ctx_t *ctx = &(state->cabac.ctx.qt_cbf_model_chroma[tr_depth]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.u, depth - 1)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.u, depth));
      }
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.v, depth - 1)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.v, depth));
      }
    }

    double bits = tr_split_bit + cbf_bits;
    split_cost += bits * state->global->cur_lambda_cost;
  } else {
    assert(width <= TR_MAX_WIDTH);
  }

  if (depth == 0 || split_cost < nosplit_cost) {
    return split_cost;
  } else {
    lcu_set_trdepth(lcu, x_px, y_px, depth, depth);

    pred_cu->cbf = nosplit_cbf;

    // We only restore the pixel data and not coefficients or cbf data.
    // The only thing we really need are the border pixels.intra_get_dir_luma_predictor
    pixels_blit(nosplit_pixels.y, lcu->rec.y, width, width, width, LCU_WIDTH);
    if (reconstruct_chroma) {
      pixels_blit(nosplit_pixels.u, lcu->rec.u, width_c, width_c, width_c, LCU_WIDTH_C);
      pixels_blit(nosplit_pixels.v, lcu->rec.v, width_c, width_c, width_c, LCU_WIDTH_C);
    }

    return nosplit_cost;
  }
}


static double luma_mode_bits(const encoder_state_t *state, int8_t luma_mode, const int8_t *intra_preds)
{
  double mode_bits;

  bool mode_in_preds = false;
  for (int i = 0; i < 3; ++i) {
    if (luma_mode == intra_preds[i]) {
      mode_in_preds = true;
    }
  }

  const cabac_ctx_t *ctx = &(state->cabac.ctx.intra_mode_model);
  mode_bits = CTX_ENTROPY_FBITS(ctx, mode_in_preds);

  if (mode_in_preds) {
    mode_bits += ((luma_mode == intra_preds[0]) ? 1 : 2);
  } else {
    mode_bits += 5;
  }

  return mode_bits;
}


static double chroma_mode_bits(const encoder_state_t *state, int8_t chroma_mode, int8_t luma_mode)
{
  const cabac_ctx_t *ctx = &(state->cabac.ctx.chroma_pred_model[0]);
  double mode_bits;
  if (chroma_mode == luma_mode) {
    mode_bits = CTX_ENTROPY_FBITS(ctx, 0);
  } else {
    mode_bits = 2.0 + CTX_ENTROPY_FBITS(ctx, 1);
  }

  return mode_bits;
}


static int8_t search_intra_chroma(encoder_state_t * const state,
                                  int x_px, int y_px, int depth,
                                  int8_t intra_mode,
                                  int8_t modes[5], int8_t num_modes,
                                  lcu_t *const lcu)
{
  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);

  if (reconstruct_chroma) {
    const vector2d_t lcu_px = { x_px & 0x3f, y_px & 0x3f };
    cu_info_t *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];

    struct {
      double cost;
      int8_t mode;
    } chroma, best_chroma;

    best_chroma.mode = 0;
    best_chroma.cost = MAX_INT;

    for (int8_t chroma_mode_i = 0; chroma_mode_i < num_modes; ++chroma_mode_i) {
      chroma.mode = modes[chroma_mode_i];

      intra_recon_lcu_chroma(state, x_px, y_px, depth, chroma.mode, NULL, lcu);
      chroma.cost = cu_rd_cost_chroma(state, lcu_px.x, lcu_px.y, depth, tr_cu, lcu);

      double mode_bits = chroma_mode_bits(state, chroma.mode, intra_mode);
      chroma.cost += mode_bits * state->global->cur_lambda_cost;

      if (chroma.cost < best_chroma.cost) {
        best_chroma = chroma;
      }
    }

    return best_chroma.mode;
  }

  return 100;
}

/**
 * \brief Sort modes and costs to ascending order according to costs.
 */
static INLINE void sort_modes(int8_t *__restrict modes, double *__restrict costs, uint8_t length)
{
  // Length is always between 5 and 23, and is either 21, 17, 9 or 8 about
  // 60% of the time, so there should be no need for anything more complex
  // than insertion sort.
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
* \brief Select mode with the smallest cost.
*/
static INLINE uint8_t select_best_mode_index(const int8_t *modes, const double *costs, uint8_t length)
{
  uint8_t best_index = 0;
  double best_cost = costs[0];
  
  for (uint8_t i = 1; i < length; ++i) {
    if (costs[i] < best_cost) {
      best_cost = costs[i];
      best_index = i;
    }
  }

  return best_index;
}

/**
 * \brief Calculate quality of the reconstruction.
 *
 * \param pred  Predicted pixels in continous memory.
 * \param orig_block  Orignal (target) pixels in continous memory.
 * \param satd_func  SATD function for this block size.
 * \param sad_func  SAD function this block size.
 * \param width  Pixel width of the block.
 *
 * \return  Estimated RD cost of the reconstruction and signaling the
 *     coefficients of the residual.
 */
static double get_cost(encoder_state_t * const state, 
                       kvz_pixel *pred, kvz_pixel *orig_block,
                       cost_pixel_nxn_func *satd_func,
                       cost_pixel_nxn_func *sad_func,
                       int width)
{
  double satd_cost = satd_func(pred, orig_block);
  if (TRSKIP_RATIO != 0 && width == 4 && state->encoder_control->trskip_enable) {
    // If the mode looks better with SAD than SATD it might be a good
    // candidate for transform skip. How much better SAD has to be is
    // controlled by TRSKIP_RATIO.

    // Add the offset bit costs of signaling 'luma and chroma use trskip',
    // versus signaling 'luma and chroma don't use trskip' to the SAD cost.
    const cabac_ctx_t *ctx = &state->cabac.ctx.transform_skip_model_luma;
    double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);
    ctx = &state->cabac.ctx.transform_skip_model_chroma;
    trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));

    double sad_cost = TRSKIP_RATIO * sad_func(pred, orig_block) + state->global->cur_lambda_cost_sqrt * trskip_bits;
    if (sad_cost < satd_cost) {
      return sad_cost;
    }
  }
  return satd_cost;
}


static void search_intra_chroma_rough(encoder_state_t * const state,
                                      int x_px, int y_px, int depth,
                                      const kvz_pixel *orig_u, const kvz_pixel *orig_v, int16_t origstride,
                                      const kvz_pixel *rec_u, const kvz_pixel *rec_v, int16_t recstride,
                                      int8_t luma_mode,
                                      int8_t modes[5], double costs[5])
{
  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);
  if (!reconstruct_chroma) return;

  const unsigned width = MAX(LCU_WIDTH_C >> depth, TR_MIN_WIDTH);

  for (int i = 0; i < 5; ++i) {
    costs[i] = 0;
  }

  cost_pixel_nxn_func *const satd_func = pixels_get_satd_func(width);
  //cost_pixel_nxn_func *const sad_func = pixels_get_sad_func(width);

  kvz_pixel _pred[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  kvz_pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);

  kvz_pixel _orig_block[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  kvz_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

  pixels_blit(orig_u, orig_block, width, width, origstride, width);
  for (int i = 0; i < 5; ++i) {
    if (modes[i] == luma_mode) continue;
    intra_get_pred(state->encoder_control, rec_u, NULL, recstride, pred, width, modes[i], 1);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    costs[i] += satd_func(pred, orig_block);
  }

  pixels_blit(orig_v, orig_block, width, width, origstride, width);
  for (int i = 0; i < 5; ++i) {
    if (modes[i] == luma_mode) continue;
    intra_get_pred(state->encoder_control, rec_v, NULL, recstride, pred, width, modes[i], 2);
    //costs[i] += get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    costs[i] += satd_func(pred, orig_block);
  }

  sort_modes(modes, costs, 5);
}


/**
 * \brief  Order the intra prediction modes according to a fast criteria.
 *
 * This function uses SATD to order the intra prediction modes. For 4x4 modes
 * SAD might be used instead, if the cost given by SAD is much better than the
 * one given by SATD, to take into account that 4x4 modes can be coded with
 * transform skip.
 *
 * The modes are searched using halving search and the total number of modes
 * that are tried is dependent on size of the predicted block. More modes
 * are tried for smaller blocks.
 *
 * \param orig  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param orig_stride  Stride of param orig.
 * \param rec  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param rec_stride  Stride of param rec.
 * \param width  Width of the prediction block.
 * \param intra_preds  Array of the 3 predicted intra modes.
 *
 * \param[out] modes  The modes ordered according to their RD costs, from best
 *     to worst. The number of modes and costs output is given by parameter
 *     modes_to_check.
 * \param[out] costs  The RD costs of corresponding modes in param modes.
 *
 * \return  Number of prediction modes in param modes.
 */
static int8_t search_intra_rough(encoder_state_t * const state, 
                                 kvz_pixel *orig, int32_t origstride,
                                 kvz_pixel *rec, int16_t recstride,
                                 int width, int8_t *intra_preds,
                                 int8_t modes[35], double costs[35])
{
  cost_pixel_nxn_func *satd_func = pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = pixels_get_sad_func(width);

  // Temporary block arrays
  kvz_pixel _pred[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  kvz_pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);
  
  kvz_pixel _orig_block[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  kvz_pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);
  
  kvz_pixel rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];

  kvz_pixel *recf = &rec_filtered_temp[recstride + 1];

  assert(width == 4 || width == 8 || width == 16 || width == 32);

  // Store original block for SAD computation
  pixels_blit(orig, orig_block, width, width, origstride, width);

  // Generate filtered reference pixels.
  {
    int16_t x, y;
    for (y = -1; y < recstride; y++) {
      recf[y*recstride - 1] = rec[y*recstride - 1];
    }
    for (x = 0; x < recstride; x++) {
      recf[x - recstride] = rec[x - recstride];
    }
    intra_filter(recf, recstride, width, 0);
  }
  
  int8_t modes_selected = 0;
  unsigned min_cost = UINT_MAX;
  unsigned max_cost = 0;
  
  // Initial offset decides how many modes are tried before moving on to the
  // recursive search.
  int offset;
  if (state->encoder_control->full_intra_search) {
    offset = 1;
  } else if (width == 4) {
    offset = 2;
  } else if (width == 8) {
    offset = 4;
  } else {
    offset = 8;
  }

  // Calculate SAD for evenly spaced modes to select the starting point for 
  // the recursive search.
  for (int mode = 2; mode <= 34; mode += offset) {
    intra_get_pred(state->encoder_control, rec, recf, recstride, pred, width, mode, 0);
    costs[modes_selected] = get_cost(state, pred, orig_block, satd_func, sad_func, width);
    modes[modes_selected] = mode;

    min_cost = MIN(min_cost, costs[modes_selected]);
    max_cost = MAX(max_cost, costs[modes_selected]);

    ++modes_selected;
  }

  int8_t best_mode = modes[select_best_mode_index(modes, costs, modes_selected)];
  double best_cost = min_cost;
  
  // Skip recursive search if all modes have the same cost.
  if (min_cost != max_cost) {
    // Do a recursive search to find the best mode, always centering on the
    // current best mode.
    while (offset > 1) {
      offset >>= 1;

      int8_t center_node = best_mode;
      int8_t mode = center_node - offset;
      if (mode >= 2) {
        intra_get_pred(state->encoder_control, rec, recf, recstride, pred, width, mode, 0);
        costs[modes_selected] = get_cost(state, pred, orig_block, satd_func, sad_func, width);
        modes[modes_selected] = mode;
        if (costs[modes_selected] < best_cost) {
          best_cost = costs[modes_selected];
          best_mode = modes[modes_selected];
        }
        ++modes_selected;
      }

      mode = center_node + offset;
      if (mode <= 34) {
        intra_get_pred(state->encoder_control, rec, recf, recstride, pred, width, mode, 0);
        costs[modes_selected] = get_cost(state, pred, orig_block, satd_func, sad_func, width);
        modes[modes_selected] = mode;
        if (costs[modes_selected] < best_cost) {
          best_cost = costs[modes_selected];
          best_mode = modes[modes_selected];
        }
        ++modes_selected;
      }
    }
  }

  int8_t add_modes[5] = {intra_preds[0], intra_preds[1], intra_preds[2], 0, 1};

  // Add DC, planar and missing predicted modes.
  for (int8_t pred_i = 0; pred_i < 5; ++pred_i) {
    bool has_mode = false;
    int8_t mode = add_modes[pred_i];

    for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
      if (modes[mode_i] == add_modes[pred_i]) {
        has_mode = true;
        break;
      }
    }

    if (!has_mode) {
      intra_get_pred(state->encoder_control, rec, recf, recstride, pred, width, mode, 0);
      costs[modes_selected] = get_cost(state, pred, orig_block, satd_func, sad_func, width);
      modes[modes_selected] = mode;
      ++modes_selected;
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  int lambda_cost = (int)(state->global->cur_lambda_cost_sqrt + 0.5);
  for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
    costs[mode_i] += lambda_cost * luma_mode_bits(state, modes[mode_i], intra_preds);
  }

  return modes_selected;
}


/**
 * \brief  Find best intra mode out of the ones listed in parameter modes.
 *
 * This function perform intra search by doing full quantization,
 * reconstruction and CABAC coding of coefficients. It is very slow
 * but results in better RD quality than using just the rough search.
 *
 * \param x_px  Luma picture coordinate.
 * \param y_px  Luma picture coordinate.
 * \param orig  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param orig_stride  Stride of param orig.
 * \param rec  Pointer to the top-left corner of current CU in the picture
 *     being encoded.
 * \param rec_stride  Stride of param rec.
 * \param intra_preds  Array of the 3 predicted intra modes.
 * \param modes_to_check  How many of the modes in param modes are checked.
 * \param[in] modes  The intra prediction modes that are to be checked.
 * 
 * \param[out] modes  The modes ordered according to their RD costs, from best
 *     to worst. The number of modes and costs output is given by parameter
 *     modes_to_check.
 * \param[out] costs  The RD costs of corresponding modes in param modes.
 * \param[out] lcu  If transform split searching is used, the transform split
 *     information for the best mode is saved in lcu.cu structure.
 */
static int8_t search_intra_rdo(encoder_state_t * const state, 
                             int x_px, int y_px, int depth,
                             kvz_pixel *orig, int32_t origstride,
                             kvz_pixel *rec, int16_t recstride,
                             int8_t *intra_preds,
                             int modes_to_check,
                             int8_t modes[35], double costs[35],
                             lcu_t *lcu)
{
  const int tr_depth = CLIP(1, MAX_PU_DEPTH, depth + state->encoder_control->tr_depth_intra);
  const int width = LCU_WIDTH >> depth;

  kvz_pixel pred[LCU_WIDTH * LCU_WIDTH + 1];
  kvz_pixel orig_block[LCU_WIDTH * LCU_WIDTH + 1];
  int rdo_mode;
  int pred_mode;

  kvz_pixel rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];
  kvz_pixel *recf = &rec_filtered_temp[recstride + 1];

  // Generate filtered reference pixels.
  {
    int x, y;
    for (y = -1; y < recstride; y++) {
      recf[y*recstride - 1] = rec[y*recstride - 1];
    }
    for (x = 0; x < recstride; x++) {
      recf[x - recstride] = rec[x - recstride];
    }
    intra_filter(recf, recstride, width, 0);
  }

  pixels_blit(orig, orig_block, width, width, origstride, width);

  // Check that the predicted modes are in the RDO mode list
  if (modes_to_check < 35) {
    for (pred_mode = 0; pred_mode < 3; pred_mode++) {
      int mode_found = 0;
      for (rdo_mode = 0; rdo_mode < modes_to_check; rdo_mode++) {
        if (intra_preds[pred_mode] == modes[rdo_mode]) {
          mode_found = 1;
          break;
        }
      }
      // Add this prediction mode to RDO checking
      if (!mode_found) {
        modes[modes_to_check] = intra_preds[pred_mode];
        modes_to_check++;
      }
    }
  }

  for(rdo_mode = 0; rdo_mode < modes_to_check; rdo_mode ++) {
    int rdo_bitcost = luma_mode_bits(state, modes[rdo_mode], intra_preds);
    costs[rdo_mode] = rdo_bitcost * (int)(state->global->cur_lambda_cost + 0.5);

    if (0 && width != 4 && tr_depth == depth) {
      // This code path has been disabled for now because it increases bdrate
      // by 1-2 %. Possibly due to not taking chroma into account during luma
      // mode search. Enabling separate chroma search compensates a little,
      // but not enough.

      // The idea for this code path is, that it would do the same thing as
      // the more general search_intra_trdepth, but would only handle cases
      // where transform split or transform skip don't need to be handled.
      intra_get_pred(state->encoder_control, rec, recf, recstride, pred, width, modes[rdo_mode], 0);
      costs[rdo_mode] += rdo_cost_intra(state, pred, orig_block, width, modes[rdo_mode], width == 4 ? 1 : 0);
    } else {
      // Perform transform split search and save mode RD cost for the best one.
      cu_info_t pred_cu;
      pred_cu.depth = depth;
      pred_cu.type = CU_INTRA;
      pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N);
      pred_cu.intra[0].mode = modes[rdo_mode];
      pred_cu.intra[1].mode = modes[rdo_mode];
      pred_cu.intra[2].mode = modes[rdo_mode];
      pred_cu.intra[3].mode = modes[rdo_mode];
      pred_cu.intra[0].mode_chroma = modes[rdo_mode];
      FILL(pred_cu.cbf, 0);

      // Reset transform split data in lcu.cu for this area.
      lcu_set_trdepth(lcu, x_px, y_px, depth, depth);

      double mode_cost = search_intra_trdepth(state, x_px, y_px, depth, tr_depth, modes[rdo_mode], MAX_INT, &pred_cu, lcu);
      costs[rdo_mode] += mode_cost;
    }
  }

  // The best transform split hierarchy is not saved anywhere, so to get the
  // transform split hierarchy the search has to be performed again with the
  // best mode.
  if (tr_depth != depth) {
    cu_info_t pred_cu;
    pred_cu.depth = depth;
    pred_cu.type = CU_INTRA;
    pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N);
    pred_cu.intra[0].mode = modes[0];
    pred_cu.intra[1].mode = modes[0];
    pred_cu.intra[2].mode = modes[0];
    pred_cu.intra[3].mode = modes[0];
    pred_cu.intra[0].mode_chroma = modes[0];
    FILL(pred_cu.cbf, 0);
    search_intra_trdepth(state, x_px, y_px, depth, tr_depth, modes[0], MAX_INT, &pred_cu, lcu);
  }

  return modes_to_check;
}


/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
static double search_cu_intra(encoder_state_t * const state,
                           const int x_px, const int y_px,
                           const int depth, lcu_t *lcu)
{
  const videoframe_t * const frame = state->tile->frame;
  const vector2d_t lcu_px = { x_px & 0x3f, y_px & 0x3f };
  const vector2d_t lcu_cu = { lcu_px.x >> 3, lcu_px.y >> 3 };
  const int8_t cu_width = (LCU_WIDTH >> (depth));
  const int cu_index = LCU_CU_OFFSET + lcu_cu.x + lcu_cu.y * LCU_T_CU_WIDTH;

  cu_info_t *cur_cu = &lcu->cu[cu_index];

  kvz_pixel rec_buffer[(LCU_WIDTH * 2 + 1) * (LCU_WIDTH * 2 + 1)];
  kvz_pixel *cu_in_rec_buffer = &rec_buffer[cu_width * 2 + 8 + 1];

  int8_t candidate_modes[3];

  cu_info_t *left_cu = 0;
  cu_info_t *above_cu = 0;

  // Select left and top CUs if they are available.
  // Top CU is not available across LCU boundary.
  if ((x_px >> 3) > 0) {
    left_cu = &lcu->cu[cu_index - 1];
  }
  if ((y_px >> 3) > 0 && lcu_cu.y != 0) {
    above_cu = &lcu->cu[cu_index - LCU_T_CU_WIDTH];
  }
  intra_get_dir_luma_predictor(x_px, y_px, candidate_modes, cur_cu, left_cu, above_cu);

  if (depth > 0) {
  // Build reconstructed block to use in prediction with extrapolated borders
  intra_build_reference_border(state->encoder_control, x_px, y_px, cu_width * 2 + 8,
                               rec_buffer, cu_width * 2 + 8, 0,
                               frame->width,
                               frame->height,
                               lcu);
  }

  int8_t modes[35];
  double costs[35];

  // Find best intra mode for 2Nx2N.
  kvz_pixel *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];
  unsigned pu_index = PU_INDEX(x_px >> 2, y_px >> 2);

  int8_t number_of_modes;
  bool skip_rough_search = (depth == 0 || state->encoder_control->rdo >= 3);
  if (!skip_rough_search) {
    number_of_modes = search_intra_rough(state,
                                              ref_pixels, LCU_WIDTH,
                                              cu_in_rec_buffer, cu_width * 2 + 8,
                                              cu_width, candidate_modes,
                                              modes, costs);
  } else {
    number_of_modes = 35;
    for (int i = 0; i < number_of_modes; ++i) {
      modes[i] = i;
      costs[i] = MAX_INT;
    }
  }

  // Set transform depth to current depth, meaning no transform splits.
  lcu_set_trdepth(lcu, x_px, y_px, depth, depth);

  // Refine results with slower search or get some results if rough search was skipped.
  if (state->encoder_control->rdo >= 2 || skip_rough_search) {
    int number_of_modes_to_search;
    if (state->encoder_control->rdo == 3) {
      number_of_modes_to_search = 35;
    } else if (state->encoder_control->rdo == 2) {
      number_of_modes_to_search = (cu_width <= 8) ? 8 : 3;
    } else {
      // Check only the predicted modes.
      number_of_modes_to_search = 0;
    }
    int num_modes_to_check = MIN(number_of_modes, number_of_modes_to_search);

    sort_modes(modes, costs, number_of_modes);
    number_of_modes = search_intra_rdo(state,
                      x_px, y_px, depth,
                      ref_pixels, LCU_WIDTH,
                      cu_in_rec_buffer, cu_width * 2 + 8,
                      candidate_modes,
                      num_modes_to_check,
                      modes, costs, lcu);
  }

  uint8_t best_mode_i = select_best_mode_index(modes, costs, number_of_modes);
  cur_cu->intra[pu_index].mode = modes[best_mode_i];

  return costs[best_mode_i];
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
      const cu_info_t *left_cu = ((x > 8) ? &cur_cu[-1] : NULL);
      const cu_info_t *above_cu = ((y > 8) ? &cur_cu[-LCU_T_CU_WIDTH] : NULL);
      intra_get_dir_luma_predictor(x, y, candidate_modes, cur_cu, left_cu, above_cu);
    }

    mode_bits = luma_mode_bits(state, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode, candidate_modes);
    if (PU_INDEX(x >> 2, y >> 2) == 0) {
      mode_bits += chroma_mode_bits(state, cur_cu->intra[0].mode_chroma, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode);
    }
  }

  return mode_bits;
}

static uint8_t get_ctx_cu_split_model(const lcu_t *lcu, int x, int y, int depth)
{
  vector2d_t lcu_cu = { (x & 0x3f) / 8, (y & 0x3f) / 8 };
  const cu_info_t *cu_array = &(lcu)->cu[LCU_CU_OFFSET];
  bool condA = x >= 8 && cu_array[(lcu_cu.x - 1) * lcu_cu.y * LCU_T_CU_WIDTH].depth > depth;
  bool condL = y >= 8 && cu_array[lcu_cu.x * (lcu_cu.y - 1) * LCU_T_CU_WIDTH].depth > depth;
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

  const vector2d_t lcu_px = { x & 0x3f, y & 0x3f };
  lcu_t *const lcu = &work_tree[depth];

  int x_local = (x&0x3f), y_local = (y&0x3f);
#ifdef _DEBUG
  int debug_split = 0;
#endif
  PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_CU);

  // Stop recursion if the CU is completely outside the frame.
  if (x >= frame->width || y >= frame->height) {
    // Return zero cost because this CU does not have to be coded.
    return 0;
  }

  cur_cu = &(&work_tree[depth])->cu[LCU_CU_OFFSET+(x_local>>3) + (y_local>>3)*LCU_T_CU_WIDTH];
  // Assign correct depth
  cur_cu->depth = depth > MAX_DEPTH ? MAX_DEPTH : depth;
  cur_cu->tr_depth = depth > 0 ? depth : 1;
  cur_cu->type = CU_NOTSET;
  cur_cu->part_size = depth > MAX_DEPTH ? SIZE_NxN : SIZE_2Nx2N;
  // If the CU is completely inside the frame at this depth, search for
  // prediction modes at this depth.
  if (x + cu_width <= frame->width &&
      y + cu_width <= frame->height)
  {

    if (state->global->slicetype != SLICE_I &&
        WITHIN(depth, ctrl->pu_depth_inter.min, ctrl->pu_depth_inter.max))
    {
      int mode_cost = search_cu_inter(state, x, y, depth, &work_tree[depth]);
      if (mode_cost < cost) {
        cost = mode_cost;
        cur_cu->type = CU_INTER;
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
      double mode_cost = search_cu_intra(state, x, y, depth, &work_tree[depth]);
      if (mode_cost < cost) {
        cost = mode_cost;
        cur_cu->type = CU_INTRA;
      }
    }

    // Reconstruct best mode because we need the reconstructed pixels for
    // mode search of adjacent CUs.
    if (cur_cu->type == CU_INTRA) {
      int8_t intra_mode = cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode;
      lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                         intra_mode,
                         intra_mode,
                         cur_cu->part_size);
      intra_recon_lcu_luma(state, x, y, depth, intra_mode, NULL, &work_tree[depth]);

      if (PU_INDEX(x >> 2, y >> 2) == 0) {
        int8_t intra_mode_chroma = intra_mode;

        // There is almost no benefit to doing the chroma mode search for
        // rd2. Possibly because the luma mode search already takes chroma
        // into account, so there is less of a chanse of luma mode being
        // really bad for chroma.
        if (state->encoder_control->rdo == 3) {
          const videoframe_t * const frame = state->tile->frame;

          double costs[5];
          int8_t modes[5] = { 0, 26, 10, 1, 34 };
          if (intra_mode != 0 && intra_mode != 26 && intra_mode != 10 && intra_mode != 1) {
            modes[4] = intra_mode;
          }

          // The number of modes to select for slower chroma search. Luma mode
          // is always one of the modes, so 2 means the final decision is made
          // between luma mode and one other mode that looks the best
          // according to search_intra_chroma_rough.
          const int8_t modes_in_depth[5] = { 1, 1, 1, 1, 2 };
          int num_modes = modes_in_depth[depth];

          if (state->encoder_control->rdo == 3) {
            num_modes = 5;
          }

          if (num_modes != 1 && num_modes != 5) {
            kvz_pixel rec_u[(LCU_WIDTH_C * 2 + 8) * (LCU_WIDTH_C * 2 + 8)];
            kvz_pixel rec_v[(LCU_WIDTH_C * 2 + 8) * (LCU_WIDTH_C * 2 + 8)];

            const int16_t width_c = MAX(LCU_WIDTH_C >> depth, TR_MIN_WIDTH);
            const int16_t rec_stride = width_c * 2 + 8;
            const int16_t out_stride = rec_stride;

            intra_build_reference_border(state->encoder_control,
                                         x, y, out_stride,
                                         rec_u, rec_stride, COLOR_U,
                                         frame->width / 2, frame->height / 2,
                                         lcu);
            intra_build_reference_border(state->encoder_control,
                                         x, y, out_stride,
                                         rec_v, rec_stride, COLOR_V,
                                         frame->width / 2, frame->height / 2,
                                         lcu);

            vector2d_t lcu_cpx = { lcu_px.x / 2, lcu_px.y / 2 };
            kvz_pixel *ref_u = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
            kvz_pixel *ref_v = &lcu->ref.v[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];

            search_intra_chroma_rough(state, x, y, depth,
                                      ref_u, ref_v, LCU_WIDTH_C,
                                      &rec_u[rec_stride + 1], &rec_v[rec_stride + 1], rec_stride,
                                      intra_mode, modes, costs);
          }

          if (num_modes > 1) {
            intra_mode_chroma = search_intra_chroma(state, x, y, depth, intra_mode, modes, num_modes, &work_tree[depth]);
            lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                               intra_mode, intra_mode_chroma,
                               cur_cu->part_size);
          }
        }

        intra_recon_lcu_chroma(state, x, y, depth, intra_mode_chroma, NULL, &work_tree[depth]);
      }
    } else if (cur_cu->type == CU_INTER) {
      // Reset transform depth because intra messes with them.
      // This will no longer be necessary if the transform depths are not shared.
      int tr_depth = depth > 0 ? depth : 1;
      lcu_set_trdepth(&work_tree[depth], x, y, depth, tr_depth);

      if (cur_cu->inter.mv_dir == 3) {
        inter_recon_lcu_bipred(state, state->global->ref->images[cur_cu->inter.mv_ref[0]], state->global->ref->images[cur_cu->inter.mv_ref[1]], x, y, LCU_WIDTH >> depth, cur_cu->inter.mv, &work_tree[depth]);
      } else {
        inter_recon_lcu(state, state->global->ref->images[cur_cu->inter.mv_ref[cur_cu->inter.mv_dir - 1]], x, y, LCU_WIDTH >> depth, cur_cu->inter.mv[cur_cu->inter.mv_dir - 1], &work_tree[depth]);
      }

      quantize_lcu_luma_residual(state, x, y, depth, NULL, &work_tree[depth]);
      quantize_lcu_chroma_residual(state, x, y, depth, NULL, &work_tree[depth]);

      int cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);

      if(cur_cu->merged && !cbf) {
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
    cost = cu_rd_cost_luma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    cost += cu_rd_cost_chroma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    double mode_bits = calc_mode_bits(state, cur_cu, x, y);
    cost += mode_bits * state->global->cur_lambda_cost;
  }
  
  // Recursively split all the way to max search depth.
  if (depth < ctrl->pu_depth_intra.max || (depth < ctrl->pu_depth_inter.max && state->global->slicetype != SLICE_I)) {
    int half_cu = cu_width / 2;
    // Using Cost = lambda * 9 to compensate on the price of the split
    double split_cost = state->global->cur_lambda_cost * CU_COST;
    int cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);
        
    if (depth < MAX_DEPTH) {
      uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);

      const cabac_ctx_t *ctx = &(state->cabac.ctx.split_flag_model[split_model]);
      cost += CTX_ENTROPY_FBITS(ctx, 0);
      split_cost += CTX_ENTROPY_FBITS(ctx, 1);
    }

    if (cur_cu->type == CU_INTRA && depth == MAX_DEPTH) {
      const cabac_ctx_t *ctx = &(state->cabac.ctx.part_size_model[0]);
      cost += CTX_ENTROPY_FBITS(ctx, 1);  // 2Nx2N
      split_cost += CTX_ENTROPY_FBITS(ctx, 0);  // NxN
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
      cu_info_t *cu_array_d1 = &(&work_tree[depth + 1])->cu[LCU_CU_OFFSET];
      cu_info_t *cu_d1 = &cu_array_d1[(lcu_cu.x + lcu_cu.y * LCU_T_CU_WIDTH)];

      // If the best CU in depth+1 is intra and the biggest it can be, try it.
      if (cu_d1->type == CU_INTRA && cu_d1->depth == depth + 1) {
        cost = 0;

        cur_cu->intra[0] = cu_d1->intra[0];
        cur_cu->type = CU_INTRA;

        lcu_set_trdepth(&work_tree[depth], x, y, depth, cur_cu->tr_depth);
        lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                           cur_cu->intra[0].mode, cur_cu->intra[0].mode_chroma,
                           cur_cu->part_size);
        intra_recon_lcu_luma(state, x, y, depth, cur_cu->intra[0].mode, NULL, &work_tree[depth]);
        intra_recon_lcu_chroma(state, x, y, depth, cur_cu->intra[0].mode_chroma, NULL, &work_tree[depth]);
        cost += cu_rd_cost_luma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
        cost += cu_rd_cost_chroma(state, x_local, y_local, depth, cur_cu, &work_tree[depth]);

        uint8_t split_model = get_ctx_cu_split_model(lcu, x, y, depth);
        const cabac_ctx_t *ctx = &(state->cabac.ctx.split_flag_model[split_model]);
        cost += CTX_ENTROPY_FBITS(ctx, 0);

        double mode_bits = calc_mode_bits(state, cur_cu, x, y);
        cost += mode_bits * state->global->cur_lambda_cost;
      }
    }

    if (split_cost < cost) {
      // Copy split modes to this depth.
      cost = split_cost;
      work_tree_copy_up(x, y, depth, work_tree);
#if _DEBUG
      debug_split = 1;
#endif
    } else if (depth > 0) {
      // Copy this CU's mode all the way down for use in adjacent CUs mode
      // search.
      work_tree_copy_down(x, y, depth, work_tree);
    }
  }
  
  PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_CU, state->encoder_control->threadqueue, "type=search_cu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d,depth=%d,split=%d,cur_cu_is_intra=%d", state->global->frame, state->tile->id, state->slice->id, 
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
  
  // Copy reference cu_info structs from neighbouring LCUs.
  {
    const int x_cu = x >> MAX_DEPTH;
    const int y_cu = y >> MAX_DEPTH;

    // Use top-left sub-cu of LCU as pointer to lcu->cu array to make things
    // simpler.
    cu_info_t *lcu_cu = &lcu->cu[LCU_CU_OFFSET];

    // Copy top CU row.
    if (y_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info_t *from_cu = videoframe_get_cu_const(frame, x_cu + i, y_cu - 1);
        cu_info_t *to_cu = &lcu_cu[i - LCU_T_CU_WIDTH];
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy left CU column.
    if (x_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info_t *from_cu = videoframe_get_cu_const(frame, x_cu - 1, y_cu + i);
        cu_info_t *to_cu = &lcu_cu[-1 + i * LCU_T_CU_WIDTH];
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy top-left CU.
    if (x_cu > 0 && y_cu > 0) {
      const cu_info_t *from_cu = videoframe_get_cu_const(frame, x_cu - 1, y_cu - 1);
      cu_info_t *to_cu = &lcu_cu[-1 - LCU_T_CU_WIDTH];
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }

    // Copy top-right CU.
    if (y_cu > 0 && x + LCU_WIDTH < frame->width) {
      const cu_info_t *from_cu = videoframe_get_cu_const(frame, x_cu + LCU_CU_WIDTH, y_cu - 1);
      cu_info_t *to_cu = &lcu->cu[LCU_T_CU_WIDTH*LCU_T_CU_WIDTH];
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
      memcpy(&lcu->top_ref.y[x_min_in_lcu], &hor_buf->y[OFFSET_HOR_BUF(x, y, frame, x_min_in_lcu-1)], x_max + (1-x_min_in_lcu));
      memcpy(&lcu->top_ref.u[x_min_in_lcu], &hor_buf->u[OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu-1)], x_max / 2 + (1-x_min_in_lcu));
      memcpy(&lcu->top_ref.v[x_min_in_lcu], &hor_buf->v[OFFSET_HOR_BUF_C(x, y, frame, x_min_in_lcu-1)], x_max / 2 + (1-x_min_in_lcu));
    }
    // Copy left reference pixels.
    if (x > 0) {
      int y_min_in_lcu = (y>0) ? 0 : 1;
      memcpy(&lcu->left_ref.y[y_min_in_lcu], &ver_buf->y[OFFSET_VER_BUF(x, y, frame, y_min_in_lcu-1)], LCU_WIDTH + (1-y_min_in_lcu));
      memcpy(&lcu->left_ref.u[y_min_in_lcu], &ver_buf->u[OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu-1)], LCU_WIDTH / 2 + (1-y_min_in_lcu));
      memcpy(&lcu->left_ref.v[y_min_in_lcu], &ver_buf->v[OFFSET_VER_BUF_C(x, y, frame, y_min_in_lcu-1)], LCU_WIDTH / 2 + (1-y_min_in_lcu));
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

    pixels_blit(&frame->source->y[x + y * frame->source->stride], lcu->ref.y,
                        x_max, y_max, frame->source->stride, LCU_WIDTH);
    pixels_blit(&frame->source->u[x_c + y_c * frame->source->stride/2], lcu->ref.u,
                        x_max_c, y_max_c, frame->source->stride/2, LCU_WIDTH / 2);
    pixels_blit(&frame->source->v[x_c + y_c * frame->source->stride/2], lcu->ref.v,
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

    // Use top-left sub-cu of LCU as pointer to lcu->cu array to make things
    // simpler.
    const cu_info_t *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];

    int x, y;
    for (y = 0; y < LCU_CU_WIDTH; ++y) {
      for (x = 0; x < LCU_CU_WIDTH; ++x) {
        const cu_info_t *from_cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
        cu_info_t *to_cu = videoframe_get_cu(frame, x_cu + x, y_cu + y);
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

    pixels_blit(lcu->rec.y, &pic->rec->y[x_px + y_px * pic->rec->stride],
                        x_max, y_max, LCU_WIDTH, pic->rec->stride);
    coefficients_blit(lcu->coeff.y, &pic->coeff_y[luma_index],
                        x_max, y_max, LCU_WIDTH, pic_width);

    pixels_blit(lcu->rec.u, &pic->rec->u[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    pixels_blit(lcu->rec.v, &pic->rec->v[(x_px / 2) + (y_px / 2) * (pic->rec->stride / 2)],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic->rec->stride / 2);
    coefficients_blit(lcu->coeff.u, &pic->coeff_u[chroma_index],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);
    coefficients_blit(lcu->coeff.v, &pic->coeff_v[chroma_index],
                        x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);
  }
}


/**
 * Search LCU for modes.
 * - Best mode gets copied to current picture.
 */
void search_lcu(encoder_state_t * const state, const int x, const int y, const yuv_t * const hor_buf, const yuv_t * const ver_buf)
{
  lcu_t work_tree[MAX_PU_DEPTH + 1];
  int depth;
  // Initialize work tree.
  for (depth = 0; depth <= MAX_PU_DEPTH; ++depth) {
    FILL(work_tree[depth], 0);
    init_lcu_t(state, x, y, &work_tree[depth], hor_buf, ver_buf);
  }

  // Start search from depth 0.
  search_cu(state, x, y, 0, work_tree);

  copy_lcu_to_cu_data(state, x, y, &work_tree[0]);
}
