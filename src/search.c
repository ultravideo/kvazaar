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

// Temporarily for debugging.
#define SEARCH_MV_FULL_RADIUS 0

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
 * This is used in the hexagon_search to select 3 points to search.
 *
 * The start of the hexagonal pattern has been repeated at the end so that
 * the indices between 1-6 can be used as the start of a 3-point list of new
 * points to search.
 *
 *   6 o-o 1 / 7
 *    /   \
 * 5 o  0  o 2 / 8
 *    \   /
 *   4 o-o 3
 */
const vector2d_t large_hexbs[10] = {
  { 0, 0 },
  { 1, -2 }, { 2, 0 }, { 1, 2 }, { -1, 2 }, { -2, 0 }, { -1, -2 },
  { 1, -2 }, { 2, 0 }
};

/**
 * This is used as the last step of the hexagon search.
 */
const vector2d_t small_hexbs[5] = {
  { 0, 0 },
  { -1, -1 }, { -1, 0 }, { 1, 0 }, { 1, 1 }
};

/*
 *  6 7 8
 *  3 4 5
 *  0 1 2
 */
const vector2d_t square[9] = {
  { -1, 1 },
  { 0, 1 }, { 1, 1 }, { -1, 0 }, { 0, 0 }, { 1, 0 }, { -1, -1 },
  { 0, -1 }, { 1, -1 }
};

static uint32_t get_ep_ex_golomb_bitcost(uint32_t symbol, uint32_t count)
{
  int32_t num_bins = 0;
  while (symbol >= (uint32_t)(1 << count)) {
    ++num_bins;
    symbol -= 1 << count;
    ++count;
  }
  num_bins ++;

  return num_bins;
}

static uint32_t get_mvd_coding_cost(vector2d_t *mvd)
{
  uint32_t bitcost = 0;
  const int32_t mvd_hor = mvd->x;
  const int32_t mvd_ver = mvd->y;
  const int8_t hor_abs_gr0 = mvd_hor != 0;
  const int8_t ver_abs_gr0 = mvd_ver != 0;
  const uint32_t mvd_hor_abs = abs(mvd_hor);
  const uint32_t mvd_ver_abs = abs(mvd_ver);

  // Greater than 0 for x/y
  bitcost += 2;

  if (hor_abs_gr0) {
    if (mvd_hor_abs > 1) {
      bitcost += get_ep_ex_golomb_bitcost(mvd_hor_abs-2, 1) - 2; // TODO: tune the costs
    }
    // Greater than 1 + sign
    bitcost += 2;
  }

  if (ver_abs_gr0) {
    if (mvd_ver_abs > 1) {
      bitcost += get_ep_ex_golomb_bitcost(mvd_ver_abs-2, 1) - 2; // TODO: tune the costs
    }
    // Greater than 1 + sign
    bitcost += 2;
  }

  return bitcost;
}

static int calc_mvd_cost(const encoder_state_t * const state, int x, int y, int mv_shift,
                         int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                         int16_t num_cand,int32_t ref_idx, uint32_t *bitcost)
{
  uint32_t temp_bitcost = 0;
  uint32_t merge_idx;
  int cand1_cost,cand2_cost;
  vector2d_t mvd_temp1, mvd_temp2;
  int8_t merged      = 0;
  int8_t cur_mv_cand = 0;

  x <<= mv_shift;
  y <<= mv_shift;

  // Check every candidate to find a match
  for(merge_idx = 0; merge_idx < (uint32_t)num_cand; merge_idx++) {
    if (merge_cand[merge_idx].dir == 3) continue;
    if (merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][0] == x &&
        merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][1] == y &&
        merge_cand[merge_idx].ref[merge_cand[merge_idx].dir - 1] == ref_idx) {
      temp_bitcost += merge_idx;
      merged = 1;
      break;
    }
  }

  // Check mvd cost only if mv is not merged
  if(!merged) {
    mvd_temp1.x = x - mv_cand[0][0];
    mvd_temp1.y = y - mv_cand[0][1];
    cand1_cost = get_mvd_coding_cost(&mvd_temp1);

    mvd_temp2.x = x - mv_cand[1][0];
    mvd_temp2.y = y - mv_cand[1][1];
    cand2_cost = get_mvd_coding_cost(&mvd_temp2);

    // Select candidate 1 if it has lower cost
    if (cand2_cost < cand1_cost) {
      cur_mv_cand = 1;
    }
    temp_bitcost += cur_mv_cand ? cand2_cost : cand1_cost;
  }
  *bitcost = temp_bitcost;
  return temp_bitcost*(int32_t)(state->global->cur_lambda_cost_sqrt+0.5);
}

unsigned tz_pattern_search(const encoder_state_t * const state, const image_t *pic, const image_t *ref, unsigned pattern_type,
                           const vector2d_t *orig, const int iDist, vector2d_t *mv, unsigned best_cost, int *best_dist,
                           int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS], int16_t num_cand, int32_t ref_idx, uint32_t *best_bitcost,
                           int block_width, int max_lcu_below)
{
  int n_points;
  int best_index = -1;
  int i;
  
  vector2d_t mv_best = { 0, 0 };

  //implemented search patterns
  vector2d_t pattern[4][8] = {
      //diamond (8 points)
      //[ ][ ][ ][ ][1][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][8][ ][ ][ ][5][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[4][ ][ ][ ][o][ ][ ][ ][2]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][7][ ][ ][ ][6][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][3][ ][ ][ ][ ]
      {
        { 0, iDist }, { iDist, 0 }, { 0, -iDist }, { -iDist, 0 },
        { iDist / 2, iDist / 2 }, { iDist / 2, -iDist / 2 }, { -iDist / 2, -iDist / 2 }, { -iDist / 2, iDist / 2 }
      },

      //square (8 points)
      //[8][ ][ ][ ][1][ ][ ][ ][2]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[7][ ][ ][ ][o][ ][ ][ ][3]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[6][ ][ ][ ][5][ ][ ][ ][4]
      {
        { 0, iDist }, { iDist, iDist }, { iDist, 0 }, { iDist, -iDist }, { 0, -iDist },
        { -iDist, -iDist }, { -iDist, 0 }, { -iDist, iDist }
      },

      //octagon (8 points)
      //[ ][ ][5][ ][ ][ ][1][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][2]
      //[4][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][o][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[8][ ][ ][ ][ ][ ][ ][ ][6]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][7][ ][ ][ ][3][ ][ ]
      {
        { iDist / 2, iDist }, { iDist, iDist / 2 }, { iDist / 2, -iDist }, { -iDist, iDist / 2 },
        { -iDist / 2, iDist }, { iDist, -iDist / 2 }, { -iDist / 2, -iDist }, { -iDist, -iDist / 2 }
      },

      //hexagon (6 points)
      //[ ][ ][5][ ][ ][ ][1][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[4][ ][ ][ ][o][ ][ ][ ][2]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][ ][ ][ ][ ][ ][ ][ ]
      //[ ][ ][6][ ][ ][ ][3][ ][ ]
      {
        { iDist / 2, iDist }, { iDist, 0 }, { iDist / 2, -iDist }, { -iDist, 0 },
        { iDist / 2, iDist }, { -iDist / 2, -iDist }, { 0, 0 }, { 0, 0 }
      }

  };

  //make sure parameter pattern_type is within correct range
  if (pattern_type > sizeof pattern - 1)
  {
    pattern_type = sizeof pattern - 1;
  }

  //set the number of points to be checked
  if (iDist == 1)
  {
    switch (pattern_type)
    {
      case 0:
        n_points = 4;
        break;
      case 2:
        n_points = 4;
        break;
      case 3:
        n_points = 4;
        break;
      default:
        n_points = 8;
        break;
    };
  }
  else
  {
    switch (pattern_type)
    {
      case 3:
        n_points = 6;
        break;
      default:
        n_points = 8;
        break;
    };
  }

  //compute SAD values for all chosen points
  for (i = 0; i < n_points; i++)
  {
    vector2d_t *current = &pattern[pattern_type][i];
    unsigned cost;
    uint32_t bitcost;

    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + current->x,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + current->y,
                            block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(state, mv->x + current->x, mv->y + current->y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + current->x,
        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + current->x + block_width,
        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + current->y,
        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + current->y + block_width);
    }

    if (cost < best_cost)
    {
      best_cost = cost;
      *best_bitcost = bitcost;
      best_index = i;
    }

  }

  if (best_index >= 0)
  {
    mv_best = pattern[pattern_type][best_index];
    *best_dist = iDist;
  }
  
  mv->x += mv_best.x;
  mv->y += mv_best.y;

  return best_cost;

}

unsigned tz_raster_search(const encoder_state_t * const state, const image_t *pic, const image_t *ref,
                          const vector2d_t *orig, vector2d_t *mv, unsigned best_cost,
                          int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS], int16_t num_cand, int32_t ref_idx, uint32_t *best_bitcost,
                          int block_width, int iSearchRange, int iRaster, int max_lcu_below)
{
  int i;
  int k;

  vector2d_t mv_best = { 0, 0 };
  
  //compute SAD values for every point in the iRaster downsampled version of the current search area
  for (i = iSearchRange; i >= -iSearchRange; i -= iRaster)
  {
    for (k = -iSearchRange; k <= iSearchRange; k += iRaster)
    {
      vector2d_t current = { k, i };
      unsigned cost;
      uint32_t bitcost;

      {
        PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
        cost = image_calc_sad(pic, ref, orig->x, orig->y,
          (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + k,
          (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + i,
          block_width, block_width, max_lcu_below);
        cost += calc_mvd_cost(state, mv->x + k, mv->y + i, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

        PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
          (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + k,
          (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + k + block_width,
          (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + i,
          (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + i + block_width);
      }

      if (cost < best_cost)
      {
        best_cost = cost;
        *best_bitcost = bitcost;
        mv_best = current;
      }

    }
  }
  
  mv->x += mv_best.x;
  mv->y += mv_best.y;

  return best_cost;

}

static unsigned tz_search(const encoder_state_t * const state, unsigned depth,
                          const image_t *pic, const image_t *ref,
                          const vector2d_t *orig, vector2d_t *mv_in_out,
                          int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                          int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{

  //TZ parameters
  int iSearchRange = 96;  // search range for each stage
  int iRaster = 5;  // search distance limit and downsampling factor for step 3                   
  unsigned step2_type = 0;  // search patterns for steps 2 and 4
  unsigned step4_type = 0;
  bool bRasterRefinementEnable = true;  // enable step 4 mode 1
  bool bStarRefinementEnable = false;   // enable step 4 mode 2 (only one mode will be executed)
  
  int block_width = CU_WIDTH_FROM_DEPTH(depth);

  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };

  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0;
  int iDist;
  int best_dist = 0;
  unsigned best_index = num_cand;
  int max_lcu_below = -1;

  if (state->encoder_control->owf) {
    max_lcu_below = 1;
  }

  //step 1, compare (0,0) vector to predicted vectors
  
  // Check whatever input vector we got, unless its (0, 0) which will be checked later.
  if (mv.x || mv.y) 
  {
    PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);

    best_cost = image_calc_sad(pic, ref, orig->x, orig->y,
                                        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                        block_width, block_width, max_lcu_below);
    best_cost += calc_mvd_cost(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost);

    PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + block_width,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + block_width);
  }

  int i;
  // Select starting point from among merge candidates. These should include
  // both mv_cand vectors and (0, 0).
  for (i = 0; i < num_cand; ++i) 
  {
    if (merge_cand[i].dir == 3) continue;
    mv.x = merge_cand[i].mv[merge_cand[i].dir - 1][0] >> 2;
    mv.y = merge_cand[i].mv[merge_cand[i].dir - 1][1] >> 2;

    PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);

	  uint32_t bitcost;
    unsigned cost = image_calc_sad(pic, ref, orig->x, orig->y,
                                   (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                   (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                   block_width, block_width, max_lcu_below);
    cost += calc_mvd_cost(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

    PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + block_width,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + block_width);

    if (cost < best_cost) {
      best_cost = cost;
      best_index = i;
      best_bitcost = bitcost;
    }
  }
  
  if (best_index < (unsigned)num_cand) {
    mv.x = merge_cand[best_index].mv[merge_cand[best_index].dir - 1][0] >> 2;
    mv.y = merge_cand[best_index].mv[merge_cand[best_index].dir - 1][1] >> 2;
  } else {
    mv.x = mv_in_out->x >> 2;
    mv.y = mv_in_out->y >> 2;
  }

  //step 2, grid search
  for (iDist = 1; iDist <= iSearchRange; iDist *= 2)
  {
    best_cost = tz_pattern_search(state, pic, ref, step2_type, orig, iDist, &mv, best_cost, &best_dist,
                                  mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_lcu_below);
  }

  //step 3, raster scan
  if (best_dist > iRaster)
  {
    best_dist = iRaster;

    best_cost = tz_raster_search(state, pic, ref, orig, &mv, best_cost, mv_cand, merge_cand, 
                                 num_cand, ref_idx, &best_bitcost, block_width, iSearchRange, iRaster, max_lcu_below);
  }

  //step 4

  //raster refinement
  if (bRasterRefinementEnable && best_dist > 0)
  {
    iDist = best_dist >> 1;
    while (iDist > 0)
    {
      best_cost = tz_pattern_search(state, pic, ref, step4_type, orig, iDist, &mv, best_cost, &best_dist,
                                   mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_lcu_below);

      iDist = iDist >> 1;
    }
  }

  //star refinement (repeat step 2 for the current starting point)
  if (bStarRefinementEnable && best_dist > 0)
  {
    for (iDist = 1; iDist <= iSearchRange; iDist *= 2)
    {
      best_cost = tz_pattern_search(state, pic, ref, step4_type, orig, iDist, &mv, best_cost, &best_dist,
                                   mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_lcu_below);
    }
  }

  mv.x = mv.x << 2;
  mv.y = mv.y << 2;

  *mv_in_out = mv;
  *bitcost_out = best_bitcost;

  return best_cost;
}

/**
 * \brief Do motion search using the HEXBS algorithm.
 *
 * \param depth      log2 depth of the search
 * \param pic        Picture motion vector is searched for.
 * \param ref        Picture motion vector is searched from.
 * \param orig       Top left corner of the searched for block.
 * \param mv_in_out  Predicted mv in and best out. Quarter pixel precision.
 *
 * \returns  Cost of the motion vector.
 *
 * Motion vector is searched by first searching iteratively with the large
 * hexagon pattern until the best match is at the center of the hexagon.
 * As a final step a smaller hexagon is used to check the adjacent pixels.
 *
 * If a non 0,0 predicted motion vector predictor is given as mv_in_out,
 * the 0,0 vector is also tried. This is hoped to help in the case where
 * the predicted motion vector is way off. In the future even more additional
 * points like 0,0 might be used, such as vectors from top or left.
 */
static unsigned hexagon_search(const encoder_state_t * const state, unsigned depth,
                               const image_t *pic, const image_t *ref,
                               const vector2d_t *orig, vector2d_t *mv_in_out,
                               int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                               int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{
  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0, bitcost;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.
  int max_lcu_below = -1;
  
  if (state->encoder_control->owf) {
    max_lcu_below = 1;
  }

  // Check mv_in, if it's not in merge candidates.
  bool mv_in_merge_cand = false;
  for (int i = 0; i < num_cand; ++i) {
    if (merge_cand[i].dir == 3) continue;
    if (merge_cand[i].mv[merge_cand[i].dir - 1][0] >> 2 == mv.x &&
        merge_cand[i].mv[merge_cand[i].dir - 1][1] >> 2 == mv.y) {
      mv_in_merge_cand = true;
      break;
    }
  }

  if (!mv_in_merge_cand) {
    PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);

    best_cost = image_calc_sad(pic, ref, orig->x, orig->y,
                                        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                        block_width, block_width, max_lcu_below);
    best_cost += calc_mvd_cost(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);
    best_bitcost = bitcost;
    best_index = num_cand; 

    PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + block_width,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + block_width);
  }

  // Select starting point from among merge candidates. These should include
  // both mv_cand vectors and (0, 0).
  for (i = 0; i < num_cand; ++i) {
    if (merge_cand[i].dir == 3) continue;
    mv.x = merge_cand[i].mv[merge_cand[i].dir - 1][0] >> 2;
    mv.y = merge_cand[i].mv[merge_cand[i].dir - 1][1] >> 2;

    PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);

    unsigned cost = image_calc_sad(pic, ref, orig->x, orig->y,
                                   (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                   (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                   block_width, block_width, max_lcu_below);
    cost += calc_mvd_cost(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

    PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + block_width,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + block_width);

    if (cost < best_cost) {
      best_cost = cost;
      best_index = i;
      best_bitcost = bitcost;
    }
  }
  if (best_index < num_cand) {
    mv.x = merge_cand[best_index].mv[merge_cand[best_index].dir - 1][0] >> 2;
    mv.y = merge_cand[best_index].mv[merge_cand[best_index].dir - 1][1] >> 2;
  } else {
    mv.x = mv_in_out->x >> 2;
    mv.y = mv_in_out->y >> 2;
  }
  
  // Search the initial 7 points of the hexagon.
  best_index = 0;
  for (i = 0; i < 7; ++i) {
    const vector2d_t *pattern = &large_hexbs[i];
    unsigned cost;
    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                             (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x, 
                             (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y,
                             block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(state, mv.x + pattern->x, mv.y + pattern->y, 2, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);

      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
                              (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x, 
                              (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x + block_width, 
                              (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y, 
                              (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y + block_width);
    }

    if (cost < best_cost) {
      best_cost    = cost;
      best_index   = i;
      best_bitcost = bitcost;
    }
  }

  // Iteratively search the 3 new points around the best match, until the best
  // match is in the center.
  while (best_index != 0) {
    unsigned start; // Starting point of the 3 offsets to be searched.
    if (best_index == 1) {
      start = 6;
    } else if (best_index == 8) {
      start = 1;
    } else {
      start = best_index - 1;
    }

    // Move the center to the best match.
    mv.x += large_hexbs[best_index].x;
    mv.y += large_hexbs[best_index].y;
    best_index = 0;

    // Iterate through the next 3 points.
    for (i = 0; i < 3; ++i) {
      const vector2d_t *offset = &large_hexbs[start + i];
      unsigned cost;
      {
        PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
        cost = image_calc_sad(pic, ref, orig->x, orig->y,
                               (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                               (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                               block_width, block_width, max_lcu_below);
        cost += calc_mvd_cost(state, mv.x + offset->x, mv.y + offset->y, 2, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
        PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=large_hexbs_iterative,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
              (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x, 
              (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x + block_width, 
              (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y, 
              (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y + block_width);
      }

      if (cost < best_cost) {
        best_cost    = cost;
        best_index   = start + i;
        best_bitcost = bitcost;
      }
      ++offset;
    }
  }

  // Move the center to the best match.
  mv.x += large_hexbs[best_index].x;
  mv.y += large_hexbs[best_index].y;
  best_index = 0;

  // Do the final step of the search with a small pattern.
  for (i = 1; i < 5; ++i) {
    const vector2d_t *offset = &small_hexbs[i];
    unsigned cost;
    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                             (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                             (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                             block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(state, mv.x + offset->x, mv.y + offset->y, 2, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, state->encoder_control->threadqueue, "type=sad,step=small_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, ref->poc - state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x, 
            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x + block_width, 
            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y, 
            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y + block_width);
    }

    if (cost > 0 && cost < best_cost) {
      best_cost    = cost;
      best_index   = i;
      best_bitcost = bitcost;
    }
  }

  // Adjust the movement vector according to the final best match.
  mv.x += small_hexbs[best_index].x;
  mv.y += small_hexbs[best_index].y;

  // Return final movement vector in quarter-pixel precision.
  mv_in_out->x = mv.x << 2;
  mv_in_out->y = mv.y << 2;

  *bitcost_out = best_bitcost;

  return best_cost;
}


#if SEARCH_MV_FULL_RADIUS
static unsigned search_mv_full(unsigned depth,
                               const picture *pic, const picture *ref,
                               const vector2d *orig, vector2d *mv_in_out,
                               int16_t mv_cand[2][2], int16_t merge_cand[MRG_MAX_NUM_CANDS][3],
                               int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{
  vector2d mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  int x, y;
  uint32_t best_bitcost = 0, bitcost;
  vector2d min_mv, max_mv;

  /*if (abs(mv.x) > SEARCH_MV_FULL_RADIUS || abs(mv.y) > SEARCH_MV_FULL_RADIUS) {
    best_cost = calc_sad(pic, ref, orig->x, orig->y,
                         orig->x, orig->y,
                         block_width, block_width);
    mv.x = 0;
    mv.y = 0;
  }*/

  min_mv.x = mv.x - SEARCH_MV_FULL_RADIUS;
  min_mv.y = mv.y - SEARCH_MV_FULL_RADIUS;
  max_mv.x = mv.x + SEARCH_MV_FULL_RADIUS;
  max_mv.y = mv.y + SEARCH_MV_FULL_RADIUS;

  for (y = min_mv.y; y < max_mv.y; ++y) {
    for (x = min_mv.x; x < max_mv.x; ++x) {
      unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                               orig->x + x,
                               orig->y + y,
                               block_width, block_width);
      cost += calc_mvd_cost(x, y, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
      if (cost < best_cost) {
        best_cost    = cost;
        best_bitcost = bitcost;
        mv.x = x;
        mv.y = y;
      }
    }
  }

  mv_in_out->x = mv.x << 2;
  mv_in_out->y = mv.y << 2;

  *bitcost_out = best_bitcost;

  return best_cost;
}
#endif

/**
 * \brief Do fractional motion estimation
 *
 * \param depth      log2 depth of the search
 * \param pic        Picture motion vector is searched for.
 * \param ref        Picture motion vector is searched from.
 * \param orig       Top left corner of the searched for block.
 * \param mv_in_out  Predicted mv in and best out. Quarter pixel precision.
 *
 * \returns  Cost of the motion vector.
 *
 * Algoritm first searches 1/2-pel positions around integer mv and after best match is found,
 * refines the search by searching best 1/4-pel postion around best 1/2-pel position.
 */
static unsigned search_frac(const encoder_state_t * const state,
                            unsigned depth,
                            const image_t *pic, const image_t *ref,
                            const vector2d_t *orig, vector2d_t *mv_in_out,
                            int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                            int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{

  //Set mv to halfpel precision
  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0, bitcost;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.

  unsigned cost = 0;

  cost_pixel_nxn_func *satd = pixels_get_satd_func(block_width);

  vector2d_t halfpel_offset;

  #define FILTER_SIZE 8
  #define HALF_FILTER (FILTER_SIZE>>1)

  //create buffer for block + extra for filter
  int src_stride = block_width+FILTER_SIZE+1;
  pixel_t src[(LCU_WIDTH+FILTER_SIZE+1) * (LCU_WIDTH+FILTER_SIZE+1)];
  pixel_t* src_off = &src[HALF_FILTER+HALF_FILTER*(block_width+FILTER_SIZE+1)];

  //destination buffer for interpolation
  int dst_stride = (block_width+1)*4;
  pixel_t dst[(LCU_WIDTH+1) * (LCU_WIDTH+1) * 16];
  pixel_t* dst_off = &dst[dst_stride*4+4];

  extend_borders(orig->x, orig->y, mv.x-1, mv.y-1,
                state->tile->lcu_offset_x * LCU_WIDTH,
                state->tile->lcu_offset_y * LCU_WIDTH,
                ref->y, ref->width, ref->height, FILTER_SIZE, block_width+1, block_width+1, src);

  filter_inter_quarterpel_luma(state->encoder_control, src_off, src_stride, block_width+1,
      block_width+1, dst, dst_stride, 1, 1);


  //Set mv to half-pixel precision
  mv.x <<= 1;
  mv.y <<= 1;

  // Search halfpel positions around best integer mv
  for (i = 0; i < 9; ++i) {
    const vector2d_t *pattern = &square[i];

    pixel_t tmp_filtered[LCU_WIDTH*LCU_WIDTH];
    pixel_t tmp_pic[LCU_WIDTH*LCU_WIDTH];

    int y,x;
    for(y = 0; y < block_width; ++y) {
      int dst_y = y*4+pattern->y*2;
      for(x = 0; x < block_width; ++x) {
        int dst_x = x*4+pattern->x*2;
        tmp_filtered[y*block_width+x] = dst_off[dst_y*dst_stride+dst_x];
        tmp_pic[y*block_width+x] = pic->y[orig->x+x + (orig->y+y)*pic->width];
      }
    }

    cost = satd(tmp_pic,tmp_filtered);

    cost += calc_mvd_cost(state, mv.x + pattern->x, mv.y + pattern->y, 1, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);

    if (cost < best_cost) {
      best_cost    = cost;
      best_index   = i;
      best_bitcost = bitcost;

    }
  }

  //Set mv to best match
  mv.x += square[best_index].x;
  mv.y += square[best_index].y;

  halfpel_offset.x = square[best_index].x*2;
  halfpel_offset.y = square[best_index].y*2;

  //Set mv to quarterpel precision
  mv.x <<= 1;
  mv.y <<= 1;

  //Search quarterpel points around best halfpel mv
  for (i = 0; i < 9; ++i) {
    const vector2d_t *pattern = &square[i];

    pixel_t tmp_filtered[LCU_WIDTH*LCU_WIDTH];
    pixel_t tmp_pic[LCU_WIDTH*LCU_WIDTH];

    int y,x;
    for(y = 0; y < block_width; ++y) {
      int dst_y = y*4+halfpel_offset.y+pattern->y;
      for(x = 0; x < block_width; ++x) {
        int dst_x = x*4+halfpel_offset.x+pattern->x;
        tmp_filtered[y*block_width+x] = dst_off[dst_y*dst_stride+dst_x];
        tmp_pic[y*block_width+x] = pic->y[orig->x+x + (orig->y+y)*pic->width];
      }
    }

    cost = satd(tmp_pic,tmp_filtered);

    cost += calc_mvd_cost(state, mv.x + pattern->x, mv.y + pattern->y, 0, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);

    if (cost < best_cost) {
      best_cost    = cost;
      best_index   = i;
      best_bitcost = bitcost;
    }
  }

  //Set mv to best final best match
  mv.x += square[best_index].x;
  mv.y += square[best_index].y;

  mv_in_out->x = mv.x;
  mv_in_out->y = mv.y;

  *bitcost_out = best_bitcost;


  return best_cost;

}

/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
static int search_cu_inter(const encoder_state_t * const state, int x, int y, int depth, lcu_t *lcu)
{
  const videoframe_t * const frame = state->tile->frame;
  uint32_t ref_idx = 0;
  int x_local = (x&0x3f), y_local = (y&0x3f);
  int x_cu = x>>3;
  int y_cu = y>>3;
  int cu_pos = LCU_CU_OFFSET+(x_local>>3) + (y_local>>3)*LCU_T_CU_WIDTH;

  cu_info_t *cur_cu = &lcu->cu[cu_pos];

  int16_t mv_cand[2][2];
  // Search for merge mode candidate
  inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS];
  // Get list of candidates
  int16_t num_cand = inter_get_merge_cand(state, x, y, depth, merge_cand, lcu);

  // Select better candidate
  cur_cu->inter.mv_cand = 0; // Default to candidate 0

  cur_cu->inter.cost = UINT_MAX;

  for (ref_idx = 0; ref_idx < state->global->ref->used_size; ref_idx++) {
    image_t *ref_image = state->global->ref->images[ref_idx];
    uint32_t temp_bitcost = 0;
    uint32_t temp_cost = 0;
    vector2d_t orig, mvd;
    int32_t merged = 0;
    uint8_t cu_mv_cand = 0;
    int8_t merge_idx = 0;
    int8_t ref_list = state->global->refmap[ref_idx].list-1;
    int8_t temp_ref_idx = cur_cu->inter.mv_ref[ref_list];
    int8_t temp_ref_list = cur_cu->inter.mv_dir;
    orig.x = x_cu * CU_MIN_SIZE_PIXELS;
    orig.y = y_cu * CU_MIN_SIZE_PIXELS;
    // Get MV candidates
    cur_cu->inter.mv_ref[ref_list] = ref_idx;
    cur_cu->inter.mv_dir = ref_list+1;
    inter_get_mv_cand(state, x, y, depth, mv_cand, cur_cu, lcu);
    cur_cu->inter.mv_ref[ref_list] = temp_ref_idx;
    cur_cu->inter.mv_dir = temp_ref_list;

    vector2d_t mv = { 0, 0 };
    {
      // Take starting point for MV search from previous frame.
      // When temporal motion vector candidates are added, there is probably
      // no point to this anymore, but for now it helps.
      int mid_x_cu = (x + (LCU_WIDTH >> (depth+1))) / 8;
      int mid_y_cu = (y + (LCU_WIDTH >> (depth+1))) / 8;
      cu_info_t *ref_cu = &state->global->ref->cu_arrays[ref_idx]->data[mid_x_cu + mid_y_cu * (frame->width_in_lcu << MAX_DEPTH)];
      if (ref_cu->type == CU_INTER) {
        if (ref_cu->inter.mv_dir & 1) {
          mv.x = ref_cu->inter.mv[0][0];
          mv.y = ref_cu->inter.mv[0][1];
        } else {
          mv.x = ref_cu->inter.mv[1][0];
          mv.y = ref_cu->inter.mv[1][1];
        }
      }
    }

#if SEARCH_MV_FULL_RADIUS
    temp_cost += search_mv_full(depth, frame, ref_pic, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
#else
    switch (state->encoder_control->cfg->ime_algorithm) {
      case IME_TZ :
        temp_cost += tz_search(state, depth, frame->source, ref_image, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
        break;

      default:
        temp_cost += hexagon_search(state, depth, frame->source, ref_image, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
        break;
      }
#endif
    if (state->encoder_control->cfg->fme_level > 0) {
      temp_cost = search_frac(state, depth, frame->source, ref_image, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
    }

    merged = 0;
    // Check every candidate to find a match
    for(merge_idx = 0; merge_idx < num_cand; merge_idx++) {
      if (merge_cand[merge_idx].dir != 3 &&
          merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][0] == mv.x &&
          merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][1] == mv.y &&          
          (uint32_t)merge_cand[merge_idx].ref == ref_idx) {
        merged = 1;
        break;
      }
    }

    // Only check when candidates are different
    if (!merged && (mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1])) {
      vector2d_t mvd_temp1, mvd_temp2;
      int cand1_cost,cand2_cost;

      mvd_temp1.x = mv.x - mv_cand[0][0];
      mvd_temp1.y = mv.y - mv_cand[0][1];
      cand1_cost = get_mvd_coding_cost(&mvd_temp1);

      mvd_temp2.x = mv.x - mv_cand[1][0];
      mvd_temp2.y = mv.y - mv_cand[1][1];
      cand2_cost = get_mvd_coding_cost(&mvd_temp2);

      // Select candidate 1 if it has lower cost
      if (cand2_cost < cand1_cost) {
        cu_mv_cand = 1;
      }
    }
    mvd.x = mv.x - mv_cand[cu_mv_cand][0];
    mvd.y = mv.y - mv_cand[cu_mv_cand][1];

    if(temp_cost < cur_cu->inter.cost) {

      // Map reference index to L0/L1 pictures
      cur_cu->inter.mv_dir = ref_list+1;
      cur_cu->inter.mv_ref_coded[ref_list] = state->global->refmap[ref_idx].idx;

      cur_cu->merged        = merged;
      cur_cu->merge_idx     = merge_idx;
      cur_cu->inter.mv_ref[ref_list] = ref_idx;
      cur_cu->inter.mv[ref_list][0] = (int16_t)mv.x;
      cur_cu->inter.mv[ref_list][1] = (int16_t)mv.y;
      cur_cu->inter.mvd[ref_list][0] = (int16_t)mvd.x;
      cur_cu->inter.mvd[ref_list][1] = (int16_t)mvd.y;
      cur_cu->inter.cost    = temp_cost;
      cur_cu->inter.bitcost = temp_bitcost + cur_cu->inter.mv_dir - 1 + cur_cu->inter.mv_ref_coded[ref_list];
      cur_cu->inter.mv_cand = cu_mv_cand;
    }
  }

  // Search bi-pred positions
  if (state->global->slicetype == SLICE_B) {
    #define NUM_PRIORITY_LIST 12;
    static const uint8_t priorityList0[] = { 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3 };
    static const uint8_t priorityList1[] = { 1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2 };
    uint8_t cutoff = num_cand;
    for (int32_t idx = 0; idx<cutoff*(cutoff - 1); idx++) {
      uint8_t i = priorityList0[idx];
      uint8_t j = priorityList1[idx];
      if (i >= num_cand || j >= num_cand) break;

      // Find one L0 and L1 candidate according to the priority list
      if ((merge_cand[i].dir & 0x1) && (merge_cand[j].dir & 0x2)) {
        if (merge_cand[i].ref[0] != merge_cand[j].ref[1] ||
          merge_cand[i].mv[0][0] != merge_cand[j].mv[1][0] ||
          merge_cand[i].mv[0][1] != merge_cand[j].mv[1][1]) {          
          int8_t cu_mv_cand = 0;
          // Force L0 and L1 references
          if (state->global->refmap[merge_cand[i].ref[0]].list == 2 || state->global->refmap[merge_cand[j].ref[1]].list == 1) continue;
          cur_cu->inter.mv_dir = 3;
          cur_cu->inter.mv_ref_coded[0] = state->global->refmap[merge_cand[i].ref[0]].idx;
          cur_cu->inter.mv_ref_coded[1] = state->global->refmap[merge_cand[j].ref[1]].idx;

          cur_cu->merged = 0;
          cur_cu->inter.mv_ref[0] = merge_cand[i].ref[0];
          cur_cu->inter.mv_ref[1] = merge_cand[j].ref[1];
          cur_cu->inter.mv[0][0] = merge_cand[i].mv[0][0];
          cur_cu->inter.mv[0][1] = merge_cand[i].mv[0][1];
          cur_cu->inter.mv[1][0] = merge_cand[j].mv[1][0];
          cur_cu->inter.mv[1][1] = merge_cand[j].mv[1][1];

          for (int reflist = 0; reflist < 2; reflist++) {
            if ((mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1])) {
              vector2d_t mvd_temp1, mvd_temp2;
              int cand1_cost, cand2_cost;

              mvd_temp1.x = cur_cu->inter.mv[reflist][0] - mv_cand[0][0];
              mvd_temp1.y = cur_cu->inter.mv[reflist][1] - mv_cand[0][1];
              cand1_cost = get_mvd_coding_cost(&mvd_temp1);

              mvd_temp2.x = cur_cu->inter.mv[reflist][0] - mv_cand[1][0];
              mvd_temp2.y = cur_cu->inter.mv[reflist][1] - mv_cand[1][1];
              cand2_cost = get_mvd_coding_cost(&mvd_temp2);

              // Select candidate 1 if it has lower cost
              if (cand2_cost < cand1_cost) {
                //cu_mv_cand = 1;
              }
            }
            cur_cu->inter.mvd[reflist][0] = cur_cu->inter.mv[reflist][0] - mv_cand[cu_mv_cand][0];
            cur_cu->inter.mvd[reflist][1] = cur_cu->inter.mv[reflist][1] - mv_cand[cu_mv_cand][1];
          }
          cur_cu->inter.cost = 0;
          cur_cu->inter.bitcost = 10 + cur_cu->inter.mv_dir - 1 + cur_cu->inter.mv_ref_coded[0] + cur_cu->inter.mv_ref_coded[1];
          cur_cu->inter.mv_cand = cu_mv_cand;
          break;
        }
      }
    }
  }

  return cur_cu->inter.cost;
}


/**
 * Copy all non-reference CU data from depth+1 to depth.
 */
static void work_tree_copy_up(int x_px, int y_px, int depth, lcu_t work_tree[MAX_PU_DEPTH + 1])
{
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
  const int width = LCU_WIDTH >> depth;
  const int width_c = width > TR_MIN_WIDTH ? width / 2 : width;

  const int offset = width / 2;
  const vector2d_t lcu_px = { x_px & 0x3f, y_px & 0x3f };
  cu_info_t *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];

  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);

  struct {
    pixel_t y[TR_MAX_WIDTH*TR_MAX_WIDTH];
    pixel_t u[TR_MAX_WIDTH*TR_MAX_WIDTH];
    pixel_t v[TR_MAX_WIDTH*TR_MAX_WIDTH];
  } nosplit_pixels;
  cu_cbf_t nosplit_cbf;

  double split_cost = INT32_MAX;
  double nosplit_cost = INT32_MAX;

  assert(width >= TR_MIN_WIDTH);

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
static INLINE int8_t select_best_mode(const int8_t *modes, const double *costs, uint8_t length)
{
  double best_mode = modes[0];
  double best_cost = costs[0];
  
  for (uint8_t i = 1; i < length; ++i) {
    if (costs[i] < best_cost) {
      best_cost = costs[i];
      best_mode = modes[i];
    }
  }

  return best_mode;
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
                       pixel_t *pred, pixel_t *orig_block,
                       cost_pixel_nxn_func *satd_func,
                       cost_pixel_nxn_func *sad_func,
                       int width)
{
  double satd_cost = satd_func(pred, orig_block);
  if (TRSKIP_RATIO != 0 && width == 4) {
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
                                      const pixel_t *orig_u, const pixel_t *orig_v, int16_t origstride,
                                      const pixel_t *rec_u, const pixel_t *rec_v, int16_t recstride,
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

  pixel_t _pred[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel_t *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);

  pixel_t _orig_block[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel_t *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);

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
                                 pixel_t *orig, int32_t origstride,
                                 pixel_t *rec, int16_t recstride,
                                 int width, int8_t *intra_preds,
                                 int8_t modes[35], double costs[35])
{
  cost_pixel_nxn_func *satd_func = pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = pixels_get_sad_func(width);

  // Temporary block arrays
  pixel_t _pred[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel_t *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);
  
  pixel_t _orig_block[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel_t *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);
  
  pixel_t rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];

  pixel_t *recf = &rec_filtered_temp[recstride + 1];

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

  int8_t best_mode = select_best_mode(modes, costs, modes_selected);
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
                             pixel_t *orig, int32_t origstride,
                             pixel_t *rec, int16_t recstride,
                             int8_t *intra_preds,
                             int modes_to_check,
                             int8_t modes[35], double costs[35],
                             lcu_t *lcu)
{
  const int tr_depth = CLIP(1, MAX_PU_DEPTH, depth + state->encoder_control->tr_depth_intra);
  const int width = LCU_WIDTH >> depth;

  pixel_t pred[LCU_WIDTH * LCU_WIDTH + 1];
  pixel_t orig_block[LCU_WIDTH * LCU_WIDTH + 1];
  int rdo_mode;
  int pred_mode;

  pixel_t rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];
  pixel_t *recf = &rec_filtered_temp[recstride + 1];

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
  for(pred_mode = 0; pred_mode < 3; pred_mode++) {
    int mode_found = 0;
    for(rdo_mode = 0; rdo_mode < modes_to_check; rdo_mode ++) {
      if(intra_preds[pred_mode] == modes[rdo_mode]) {
        mode_found = 1;
        break;
      }
    }
    // Add this prediction mode to RDO checking
    if(!mode_found) {
      modes[modes_to_check] = intra_preds[pred_mode];
      modes_to_check++;
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

  pixel_t rec_buffer[(LCU_WIDTH * 2 + 1) * (LCU_WIDTH * 2 + 1)];
  pixel_t *cu_in_rec_buffer = &rec_buffer[cu_width * 2 + 8 + 1];

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
  {
    pixel_t *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];
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

    int8_t best_mode = select_best_mode(modes, costs, number_of_modes);

    cur_cu->intra[pu_index].mode = best_mode;
  }

  return costs[0];
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
static double search_cu(encoder_state_t * const state, int x, int y, int depth, lcu_t work_tree[MAX_PU_DEPTH])
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
            pixel_t rec_u[(LCU_WIDTH_C * 2 + 8) * (LCU_WIDTH_C * 2 + 8)];
            pixel_t rec_v[(LCU_WIDTH_C * 2 + 8) * (LCU_WIDTH_C * 2 + 8)];

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
            pixel_t *ref_u = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];
            pixel_t *ref_v = &lcu->ref.u[lcu_cpx.x + lcu_cpx.y * LCU_WIDTH_C];

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
        pixel_t *temp_lcu_y = MALLOC(pixel_t, 64 * 64);
        pixel_t *temp_lcu_u = MALLOC(pixel_t, 32 * 32);
        pixel_t *temp_lcu_v = MALLOC(pixel_t, 32 * 32);
        int temp_x, temp_y;
        inter_recon_lcu(state, state->global->ref->images[cur_cu->inter.mv_ref[0]], x, y, LCU_WIDTH >> depth, cur_cu->inter.mv[0], &work_tree[depth]);   
        memcpy(temp_lcu_y, work_tree[depth].rec.y, sizeof(pixel_t) * 64 * 64);
        memcpy(temp_lcu_u, work_tree[depth].rec.u, sizeof(pixel_t) * 32 * 32);
        memcpy(temp_lcu_v, work_tree[depth].rec.v, sizeof(pixel_t) * 32 * 32);
        inter_recon_lcu(state, state->global->ref->images[cur_cu->inter.mv_ref[1]], x, y, LCU_WIDTH >> depth, cur_cu->inter.mv[1], &work_tree[depth]);
        for (temp_y = 0; temp_y < LCU_WIDTH >> depth; ++temp_y) {
          int y_in_lcu = ((y + temp_y) & ((LCU_WIDTH)-1));
          for (temp_x = 0; temp_x < LCU_WIDTH >> depth; ++temp_x) {
            int x_in_lcu = ((x + temp_x) & ((LCU_WIDTH)-1));
            lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = (pixel_t)(((int)lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] + 
                                                                    (int)temp_lcu_y[y_in_lcu * LCU_WIDTH + x_in_lcu]) >> 1);
          }
        }
        for (temp_y = 0; temp_y < LCU_WIDTH >> (depth+1); ++temp_y) {
          int y_in_lcu = ((y + temp_y) & ((LCU_WIDTH>>1)-1));
          for (temp_x = 0; temp_x < LCU_WIDTH >> (depth+1); ++temp_x) {
            int x_in_lcu = ((x + temp_x) & ((LCU_WIDTH>>1)-1));
            lcu->rec.u[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu] = (pixel_t)((int)(lcu->rec.u[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu] +
                                                                       (int)temp_lcu_u[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu]) >> 1);
            lcu->rec.v[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu] = (pixel_t)(((int)lcu->rec.v[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu] +
                                                                       (int)temp_lcu_v[y_in_lcu * (LCU_WIDTH >> 1) + x_in_lcu]) >> 1);
          }
        }
        FREE_POINTER(temp_lcu_y);
        FREE_POINTER(temp_lcu_u);
        FREE_POINTER(temp_lcu_v);
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
