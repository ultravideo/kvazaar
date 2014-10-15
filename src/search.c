/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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
const vector2d large_hexbs[10] = {
  { 0, 0 },
  { 1, -2 }, { 2, 0 }, { 1, 2 }, { -1, 2 }, { -2, 0 }, { -1, -2 },
  { 1, -2 }, { 2, 0 }
};

/**
 * This is used as the last step of the hexagon search.
 */
const vector2d small_hexbs[5] = {
  { 0, 0 },
  { -1, -1 }, { -1, 0 }, { 1, 0 }, { 1, 1 }
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

static uint32_t get_mvd_coding_cost(vector2d *mvd)
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

static int calc_mvd_cost(const encoder_state * const encoder_state, int x, int y,
                         int16_t mv_cand[2][2], int16_t merge_cand[MRG_MAX_NUM_CANDS][3],
                         int16_t num_cand,int32_t ref_idx, uint32_t *bitcost)
{
  uint32_t temp_bitcost = 0;
  uint32_t merge_idx;
  int cand1_cost,cand2_cost;
  vector2d mvd_temp1, mvd_temp2;
  int8_t merged      = 0;
  int8_t cur_mv_cand = 0;

  x <<= 2;
  y <<= 2;

  // Check every candidate to find a match
  for(merge_idx = 0; merge_idx < (uint32_t)num_cand; merge_idx++) {
    if (merge_cand[merge_idx][0] == x &&
        merge_cand[merge_idx][1] == y &&
        merge_cand[merge_idx][2] == ref_idx) {
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
  return temp_bitcost*(int32_t)(encoder_state->global->cur_lambda_cost_sqrt+0.5);
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
static unsigned hexagon_search(const encoder_state * const encoder_state, unsigned depth,
                               const image *pic, const image *ref,
                               const vector2d *orig, vector2d *mv_in_out,
                               int16_t mv_cand[2][2], int16_t merge_cand[MRG_MAX_NUM_CANDS][3],
                               int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{
  vector2d mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0, bitcost;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.
  int max_lcu_below = -1;
  
  if (encoder_state->encoder_control->owf) {
    max_lcu_below = 1;
  }
  
  // Search the initial 7 points of the hexagon.
  for (i = 0; i < 7; ++i) {
    const vector2d *pattern = &large_hexbs[i];
    unsigned cost;
    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                             (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x, 
                             (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y,
                             block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(encoder_state, mv.x + pattern->x, mv.y + pattern->y, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);

      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, encoder_state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", encoder_state->global->frame, encoder_state->tile->id, ref->poc - encoder_state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
                              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x, 
                              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x + block_width, 
                              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y, 
                              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y + block_width);
    }

    if (cost < best_cost) {
      best_cost    = cost;
      best_index   = i;
      best_bitcost = bitcost;
    }
  }

  // Try the 0,0 vector.
  if (!(mv.x == 0 && mv.y == 0)) {
    unsigned cost;
    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                             (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x, 
                             (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y,
                             block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(encoder_state, 0, 0, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, encoder_state->encoder_control->threadqueue, "type=sad,step=00vector,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", encoder_state->global->frame, encoder_state->tile->id, ref->poc - encoder_state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
                              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x, 
                              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + block_width, 
                              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y, 
                              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + block_width);
    }

    // If the 0,0 is better, redo the hexagon around that point.
    if (cost < best_cost) {
      best_cost    = cost;
      best_bitcost = bitcost;
      best_index   = 0;
      mv.x = 0;
      mv.y = 0;

      for (i = 1; i < 7; ++i) {
        const vector2d *pattern = &large_hexbs[i];
        unsigned cost;
        {
          PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
          cost = image_calc_sad(pic, ref, orig->x, orig->y,
                                 (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + pattern->x,
                                 (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + pattern->y,
                                 block_width, block_width, max_lcu_below);
          cost += calc_mvd_cost(encoder_state, pattern->x, pattern->y, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
          PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, encoder_state->encoder_control->threadqueue, "type=sad,step=large_hexbs_around00,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", encoder_state->global->frame, encoder_state->tile->id, ref->poc - encoder_state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
                        (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + pattern->x, 
                        (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + pattern->x + block_width, 
                        (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + pattern->y, 
                        (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + pattern->y + block_width);
        }

        if (cost < best_cost) {
          best_cost    = cost;
          best_index   = i;
          best_bitcost = bitcost;
        }
      }
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
      const vector2d *offset = &large_hexbs[start + i];
      unsigned cost;
      {
        PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
        cost = image_calc_sad(pic, ref, orig->x, orig->y,
                               (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                               (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                               block_width, block_width, max_lcu_below);
        cost += calc_mvd_cost(encoder_state, mv.x + offset->x, mv.y + offset->y, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
        PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, encoder_state->encoder_control->threadqueue, "type=sad,step=large_hexbs_iterative,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", encoder_state->global->frame, encoder_state->tile->id, ref->poc - encoder_state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x, 
              (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x + block_width, 
              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y, 
              (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y + block_width);
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
    const vector2d *offset = &small_hexbs[i];
    unsigned cost;
    {
      PERFORMANCE_MEASURE_START(_DEBUG_PERF_SEARCH_PIXELS);
      cost = image_calc_sad(pic, ref, orig->x, orig->y,
                             (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                             (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                             block_width, block_width, max_lcu_below);
      cost += calc_mvd_cost(encoder_state, mv.x + offset->x, mv.y + offset->y, mv_cand,merge_cand,num_cand,ref_idx, &bitcost);
      PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_PIXELS, encoder_state->encoder_control->threadqueue, "type=sad,step=small_hexbs,frame=%d,tile=%d,ref=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", encoder_state->global->frame, encoder_state->tile->id, ref->poc - encoder_state->global->poc, orig->x, orig->x + block_width, orig->y, orig->y + block_width, 
            (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x, 
            (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x + block_width, 
            (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y, 
            (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y + block_width);
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
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
static int search_cu_inter(const encoder_state * const encoder_state, int x, int y, int depth, lcu_t *lcu)
{
  const videoframe * const frame = encoder_state->tile->frame;
  uint32_t ref_idx = 0;
  int x_local = (x&0x3f), y_local = (y&0x3f);
  int x_cu = x>>3;
  int y_cu = y>>3;
  int cu_pos = LCU_CU_OFFSET+(x_local>>3) + (y_local>>3)*LCU_T_CU_WIDTH;

  cu_info *cur_cu = &lcu->cu[cu_pos];

  int16_t mv_cand[2][2];
  // Search for merge mode candidate
  int16_t merge_cand[MRG_MAX_NUM_CANDS][3];
  // Get list of candidates
  int16_t num_cand = inter_get_merge_cand(x, y, depth, merge_cand, lcu);

  // Select better candidate
  cur_cu->inter.mv_cand = 0; // Default to candidate 0

  cur_cu->inter.cost = UINT_MAX;

  for (ref_idx = 0; ref_idx < encoder_state->global->ref->used_size; ref_idx++) {
    image *ref_image = encoder_state->global->ref->images[ref_idx];
    const cu_info *ref_cu = &encoder_state->global->ref->cu_arrays[ref_idx]->data[x_cu + y_cu * (frame->width_in_lcu << MAX_DEPTH)];
    uint32_t temp_bitcost = 0;
    uint32_t temp_cost = 0;
    vector2d orig, mv, mvd;
    int32_t merged = 0;
    uint8_t cu_mv_cand = 0;
    int8_t merge_idx = 0;
    int8_t temp_ref_idx = cur_cu->inter.mv_ref;
    orig.x = x_cu * CU_MIN_SIZE_PIXELS;
    orig.y = y_cu * CU_MIN_SIZE_PIXELS;
    mv.x = 0;
    mv.y = 0;
    if (ref_cu->type == CU_INTER) {
      mv.x = ref_cu->inter.mv[0];
      mv.y = ref_cu->inter.mv[1];
    }
    // Get MV candidates
    cur_cu->inter.mv_ref = ref_idx;
    inter_get_mv_cand(encoder_state, x, y, depth, mv_cand, cur_cu, lcu);
    cur_cu->inter.mv_ref = temp_ref_idx;

#if SEARCH_MV_FULL_RADIUS
    temp_cost += search_mv_full(depth, frame, ref_pic, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
#else
    temp_cost += hexagon_search(encoder_state, depth, frame->source, ref_image, &orig, &mv, mv_cand, merge_cand, num_cand, ref_idx, &temp_bitcost);
#endif

    merged = 0;
    // Check every candidate to find a match
    for(merge_idx = 0; merge_idx < num_cand; merge_idx++) {
      if (merge_cand[merge_idx][0] == mv.x &&
          merge_cand[merge_idx][1] == mv.y &&
          (uint32_t)merge_cand[merge_idx][2] == ref_idx) {
        merged = 1;
        break;
      }
    }

    // Only check when candidates are different
    if (!merged && (mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1])) {
      vector2d mvd_temp1, mvd_temp2;
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
      cur_cu->merged        = merged;
      cur_cu->merge_idx     = merge_idx;
      cur_cu->inter.mv_ref  = ref_idx;
      cur_cu->inter.mv_dir  = 1;
      cur_cu->inter.mv[0]   = (int16_t)mv.x;
      cur_cu->inter.mv[1]   = (int16_t)mv.y;
      cur_cu->inter.mvd[0]  = (int16_t)mvd.x;
      cur_cu->inter.mvd[1]  = (int16_t)mvd.y;
      cur_cu->inter.cost    = temp_cost;
      cur_cu->inter.bitcost = temp_bitcost + ref_idx;
      cur_cu->inter.mv_cand = cu_mv_cand;
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
        const cu_info *from_cu = &work_tree[depth + 1].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
        cu_info *to_cu = &work_tree[depth].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
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
        const cu_info *from_cu = &work_tree[depth].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
        cu_info *to_cu = &work_tree[d].cu[LCU_CU_OFFSET + x + y * LCU_T_CU_WIDTH];
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
  const vector2d lcu_cu = { (x_px & (LCU_WIDTH - 1)) / 8, (y_px & (LCU_WIDTH - 1)) / 8 };
  cu_info *const cur_cu = &lcu->cu[lcu_cu.x + lcu_cu.y * LCU_T_CU_WIDTH + LCU_CU_OFFSET];
  int x, y;

  // Depth 4 doesn't go inside the loop. Set the top-left CU.
  cur_cu->tr_depth = tr_depth;

  for (y = 0; y < width_cu; ++y) {
    for (x = 0; x < width_cu; ++x) {
      cu_info *cu = &cur_cu[x + y * LCU_T_CU_WIDTH];
      cu->tr_depth = tr_depth;
    }
  }
}


static void lcu_set_intra_mode(lcu_t *lcu, int x_px, int y_px, int depth, int pred_mode, int chroma_mode, int part_mode)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;

  // NxN can only be applied to a single CU at a time.
  if (part_mode == SIZE_NxN) {
    cu_info *cu = &lcu_cu[x_cu + y_cu * LCU_T_CU_WIDTH];
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
      cu_info *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
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


static void lcu_set_inter(lcu_t *lcu, int x_px, int y_px, int depth, cu_info *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;
  // Set mode in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
      //Check if this could be moved inside the if
      cu->coded    = 1;
      if (cu != cur_cu) {
        cu->depth    = cur_cu->depth;
        cu->type     = CU_INTER;
        cu->tr_depth = cur_cu->tr_depth;
        cu->merged   = cur_cu->merged;
        cu->skipped  = cur_cu->skipped;
        memcpy(&cu->inter, &cur_cu->inter, sizeof(cu_info_inter));
      }
    }
  }
}


static void lcu_set_coeff(lcu_t *lcu, int x_px, int y_px, int depth, cu_info *cur_cu)
{
  const int width_cu = LCU_CU_WIDTH >> depth;
  const int x_cu = SUB_SCU(x_px) >> MAX_DEPTH;
  const int y_cu = SUB_SCU(y_px) >> MAX_DEPTH;
  cu_info *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];
  int x, y;
  int tr_split = cur_cu->tr_depth-cur_cu->depth;

  // Set coeff flags in every CU covered by part_mode in this depth.
  for (y = y_cu; y < y_cu + width_cu; ++y) {
    for (x = x_cu; x < x_cu + width_cu; ++x) {
      cu_info *cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
      // Use TU top-left CU to propagate coeff flags
      uint32_t mask = ~((width_cu>>tr_split)-1);
      cu_info *cu_from = &lcu_cu[(x & mask) + (y & mask) * LCU_T_CU_WIDTH];
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
static double cu_rd_cost_luma(const encoder_state *const encoder_state,
  const int x_px, const int y_px, const int depth,
  const cu_info *const pred_cu,
  lcu_t *const lcu)
{
  const int rdo = encoder_state->encoder_control->rdo;
  const int width = LCU_WIDTH >> depth;
  const uint8_t pu_index = PU_INDEX(x_px / 4, y_px / 4);

  // cur_cu is used for TU parameters.
  cu_info *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (x_px / 8) + (y_px / 8) * LCU_T_CU_WIDTH];

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
    const cabac_ctx *ctx = &(encoder_state->cabac.ctx.trans_subdiv_model[5 - (6 - depth)]);
    tr_tree_bits += CTX_ENTROPY_FBITS(ctx, tr_depth > 0);
  }

  if (tr_depth > 0) {
    int offset = width / 2;
    double sum = 0;

    sum += cu_rd_cost_luma(encoder_state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(encoder_state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(encoder_state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_luma(encoder_state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * encoder_state->global->cur_lambda_cost;
  }

  // Add transform_tree cbf_luma bit cost.
  if (pred_cu->type == CU_INTRA ||
      tr_depth > 0 ||
      cbf_is_set(tr_cu->cbf.u, depth) ||
      cbf_is_set(tr_cu->cbf.v, depth))
  {
    const cabac_ctx *ctx = &(encoder_state->cabac.ctx.qt_cbf_model_luma[!tr_depth]);
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

  if (rdo >= 1) {
    coefficient coeff_temp[32 * 32];
    int8_t luma_scan_mode = get_scan_order(pred_cu->type, pred_cu->intra[PU_INDEX(x_px / 4, y_px / 4)].mode, depth);

    // Code coeffs using cabac to get a better estimate of real coding costs.
    coefficients_blit(&lcu->coeff.y[(y_px*LCU_WIDTH) + x_px], coeff_temp, width, width, LCU_WIDTH, width);
    coeff_bits += get_coeff_cost(encoder_state, coeff_temp, width, 0, luma_scan_mode);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * LUMA_MULT + bits * encoder_state->global->cur_lambda_cost;
}


static double cu_rd_cost_chroma(const encoder_state *const encoder_state,
  const int x_px, const int y_px, const int depth,
  const cu_info *const pred_cu,
  lcu_t *const lcu)
{
  const vector2d lcu_px = { x_px / 2, y_px / 2 };
  const int rdo = encoder_state->encoder_control->rdo;
  const int width = (depth <= MAX_DEPTH) ? LCU_WIDTH >> (depth + 1) : LCU_WIDTH >> depth;
  cu_info *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x / 4) + (lcu_px.y / 4)*LCU_T_CU_WIDTH];

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
    const cabac_ctx *ctx = &(encoder_state->cabac.ctx.qt_cbf_model_chroma[tr_depth]);
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

    sum += cu_rd_cost_chroma(encoder_state, x_px, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(encoder_state, x_px + offset, y_px, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(encoder_state, x_px, y_px + offset, depth + 1, pred_cu, lcu);
    sum += cu_rd_cost_chroma(encoder_state, x_px + offset, y_px + offset, depth + 1, pred_cu, lcu);

    return sum + tr_tree_bits * encoder_state->global->cur_lambda_cost;
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

  if (rdo >= 1) {
    coefficient coeff_temp[16 * 16];
    int8_t scan_order = get_scan_order(pred_cu->type, pred_cu->intra[0].mode_chroma, depth);
    
    coefficients_blit(&lcu->coeff.u[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += get_coeff_cost(encoder_state, coeff_temp, width, 2, scan_order);

    coefficients_blit(&lcu->coeff.v[(lcu_px.y*(LCU_WIDTH_C)) + lcu_px.x],
                      coeff_temp, width, width, LCU_WIDTH_C, width);
    coeff_bits += get_coeff_cost(encoder_state, coeff_temp, width, 2, scan_order);
  }

  double bits = tr_tree_bits + coeff_bits;
  return (double)ssd * CHROMA_MULT + bits * encoder_state->global->cur_lambda_cost;
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
static double search_intra_trdepth(encoder_state * const encoder_state,
  int x_px, int y_px, int depth, int max_depth,
  int intra_mode, int cost_treshold,
  cu_info *const pred_cu,
  lcu_t *const lcu)
{
  const int width = LCU_WIDTH >> depth;
  const int width_c = width > TR_MIN_WIDTH ? width / 2 : width;

  const int offset = width / 2;
  const vector2d lcu_px = { x_px & 0x3f, y_px & 0x3f };
  cu_info *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];

  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);

  struct {
    pixel y[TR_MAX_WIDTH*TR_MAX_WIDTH];
    pixel u[TR_MAX_WIDTH*TR_MAX_WIDTH];
    pixel v[TR_MAX_WIDTH*TR_MAX_WIDTH];
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

    intra_recon_lcu_luma(encoder_state, x_px, y_px, depth, intra_mode, pred_cu, lcu);
    nosplit_cost += cu_rd_cost_luma(encoder_state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);

    if (reconstruct_chroma) {
      cbf_clear(&pred_cu->cbf.u, depth);
      cbf_clear(&pred_cu->cbf.v, depth);

      intra_recon_lcu_chroma(encoder_state, x_px, y_px, depth, intra_mode, pred_cu, lcu);
      nosplit_cost += cu_rd_cost_chroma(encoder_state, lcu_px.x, lcu_px.y, depth, pred_cu, lcu);
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
    split_cost = 3 * encoder_state->global->cur_lambda_cost;

    split_cost += search_intra_trdepth(encoder_state, x_px, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(encoder_state, x_px + offset, y_px, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(encoder_state, x_px, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }
    if (split_cost < nosplit_cost) {
      split_cost += search_intra_trdepth(encoder_state, x_px + offset, y_px + offset, depth + 1, max_depth, intra_mode, nosplit_cost, pred_cu, lcu);
    }

    double tr_split_bit = 0.0;
    double cbf_bits = 0.0;

    // Add bits for split_transform_flag = 1, because transform depth search bypasses
    // the normal recursion in the cost functions.
    if (depth >= 1 && depth <= 3) {
      const cabac_ctx *ctx = &(encoder_state->cabac.ctx.trans_subdiv_model[5 - (6 - depth)]);
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

      const cabac_ctx *ctx = &(encoder_state->cabac.ctx.qt_cbf_model_chroma[tr_depth]);
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.u, depth - 1)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.u, depth));
      }
      if (tr_depth == 0 || cbf_is_set(pred_cu->cbf.v, depth - 1)) {
        cbf_bits += CTX_ENTROPY_FBITS(ctx, cbf_is_set(pred_cu->cbf.v, depth));
      }
    }

    double bits = tr_split_bit + cbf_bits;
    split_cost += bits * encoder_state->global->cur_lambda_cost;
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


static double luma_mode_bits(const encoder_state *encoder_state, int8_t luma_mode, const int8_t *intra_preds)
{
  double mode_bits;

  bool mode_in_preds = false;
  for (int i = 0; i < 3; ++i) {
    if (luma_mode == intra_preds[i]) {
      mode_in_preds = true;
    }
  }

  const cabac_ctx *ctx = &(encoder_state->cabac.ctx.intra_mode_model);
  mode_bits = CTX_ENTROPY_FBITS(ctx, mode_in_preds);

  if (mode_in_preds) {
    mode_bits += ((luma_mode == intra_preds[0]) ? 1 : 2);
  } else {
    mode_bits += 5;
  }

  return mode_bits;
}


static double chroma_mode_bits(const encoder_state *encoder_state, int8_t chroma_mode, int8_t luma_mode)
{
  const cabac_ctx *ctx = &(encoder_state->cabac.ctx.chroma_pred_model[0]);
  double mode_bits;
  if (chroma_mode == luma_mode) {
    mode_bits = CTX_ENTROPY_FBITS(ctx, 0);
  } else {
    mode_bits = 2.0 + CTX_ENTROPY_FBITS(ctx, 1);
  }

  return mode_bits;
}


static int8_t search_intra_chroma(encoder_state * const encoder_state,
                                int x_px, int y_px, int depth,
                                int8_t intra_mode,
                                lcu_t *const lcu)
{
  const bool reconstruct_chroma = !(x_px & 4 || y_px & 4);

  if (reconstruct_chroma) {
    const vector2d lcu_px = { x_px & 0x3f, y_px & 0x3f };
    cu_info *const tr_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];

    int8_t chroma_modes[5] = { 0, 26, 10, 1, intra_mode };
    const int8_t num_chroma_modes = 5;

    if (intra_mode == 0 || intra_mode == 26 || intra_mode == 10 || intra_mode == 1) {
      chroma_modes[4] = 34;
    }

    struct {
      double cost;
      int8_t mode;
    } chroma, best_chroma;

    best_chroma.mode = 0;
    best_chroma.cost = MAX_INT;

    for (int8_t chroma_mode_i = 0; chroma_mode_i < num_chroma_modes; ++chroma_mode_i) {
      chroma.mode = chroma_modes[chroma_mode_i];

      intra_recon_lcu_chroma(encoder_state, x_px, y_px, depth, chroma.mode, NULL, lcu);
      chroma.cost = cu_rd_cost_chroma(encoder_state, lcu_px.x, lcu_px.y, depth, tr_cu, lcu);

      const cabac_ctx *ctx = &(encoder_state->cabac.ctx.chroma_pred_model[0]);
      double mode_bits = chroma_mode_bits(encoder_state, chroma.mode, intra_mode);
      chroma.cost += mode_bits * encoder_state->global->cur_lambda_cost;

      if (chroma.cost < best_chroma.cost) {
        best_chroma = chroma;
      }
    }

    return best_chroma.mode;
  }

  return 100;
}


static void sort_modes(int8_t *modes, double *costs, int length)
{
  int i, j;
  for (i = 0; i < length; ++i) {
    j = i;
    while (j > 0 && costs[j] < costs[j - 1]) {
      SWAP(costs[j], costs[j - 1], double);
      SWAP(modes[j], modes[j - 1], int8_t);
      --j;
    }
  }
}


static double get_cost(encoder_state * const encoder_state, pixel *pred, pixel *orig_block, cost_pixel_nxn_func *satd_func, cost_pixel_nxn_func *sad_func, int width)
{
  double satd_cost = satd_func(pred, orig_block);
  if (TRSKIP_RATIO != 0 && width == 4) {
    // If the mode looks better with SAD than SATD it might be a good
    // candidate for transform skip. How much better SAD has to be is
    // controlled by TRSKIP_RATIO.
    const cabac_ctx *ctx = &encoder_state->cabac.ctx.transform_skip_model_luma;
    double trskip_bits = CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0);
    ctx = &encoder_state->cabac.ctx.transform_skip_model_chroma;
    trskip_bits += 2.0 * (CTX_ENTROPY_FBITS(ctx, 1) - CTX_ENTROPY_FBITS(ctx, 0));
    double sad_cost = TRSKIP_RATIO * sad_func(pred, orig_block) + encoder_state->global->cur_lambda_cost_sqrt * trskip_bits;
    if (sad_cost < satd_cost) {
      return sad_cost;
    }
  }
  return satd_cost;
}


static int8_t search_intra_rough(encoder_state * const encoder_state, 
                                 pixel *orig, int32_t origstride,
                                 pixel *rec, int16_t recstride,
                                 int width, int8_t *intra_preds,
                                 int8_t modes[35], double costs[35])
{
  cost_pixel_nxn_func *satd_func = pixels_get_satd_func(width);
  cost_pixel_nxn_func *sad_func = pixels_get_sad_func(width);

  // Temporary block arrays
  pixel _pred[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel *pred = ALIGNED_POINTER(_pred, SIMD_ALIGNMENT);
  
  pixel _orig_block[LCU_WIDTH * LCU_WIDTH + 1 + SIMD_ALIGNMENT];
  pixel *orig_block = ALIGNED_POINTER(_orig_block, SIMD_ALIGNMENT);
  
  pixel rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];

  pixel *ref[2] = {rec, &rec_filtered_temp[recstride + 1]};

  assert(width == 4 || width == 8 || width == 16 || width == 32);

  // Store original block for SAD computation
  pixels_blit(orig, orig_block, width, width, origstride, width);

  // Generate filtered reference pixels.
  {
    int16_t x, y;
    for (y = -1; y < recstride; y++) {
      ref[1][y*recstride - 1] = rec[y*recstride - 1];
    }
    for (x = 0; x < recstride; x++) {
      ref[1][x - recstride] = rec[x - recstride];
    }
    intra_filter(ref[1], recstride, width, 0);
  }
  
  int8_t modes_selected = 0;
  unsigned min_cost = UINT_MAX;
  unsigned max_cost = 0;
  
  // Initial offset decides how many modes are tried before moving on to the
  // recursive search.
  int offset;
  if (encoder_state->encoder_control->full_intra_search) {
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
    intra_get_pred(encoder_state->encoder_control, ref, recstride, pred, width, mode, 0);
    costs[modes_selected] = get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
    modes[modes_selected] = mode;

    min_cost = MIN(min_cost, costs[modes_selected]);
    max_cost = MAX(max_cost, costs[modes_selected]);

    ++modes_selected;
  }
  
  // Skip recursive search if all modes have the same cost.
  if (min_cost != max_cost) {
    // Do a recursive search to find the best mode, always centering on the
    // current best mode.
    while (offset > 1) {
      offset >>= 1;
      sort_modes(modes, costs, modes_selected);

      int8_t mode = modes[0] - offset;
      if (mode >= 2) {
        intra_get_pred(encoder_state->encoder_control, ref, recstride, pred, width, mode, 0);
        costs[modes_selected] = get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
        modes[modes_selected] = mode;
        ++modes_selected;
      }

      mode = modes[0] + offset;
      if (mode <= 34) {
        intra_get_pred(encoder_state->encoder_control, ref, recstride, pred, width, mode, 0);
        costs[modes_selected] = get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
        modes[modes_selected] = mode;
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
      intra_get_pred(encoder_state->encoder_control, ref, recstride, pred, width, mode, 0);
      costs[modes_selected] = get_cost(encoder_state, pred, orig_block, satd_func, sad_func, width);
      modes[modes_selected] = mode;
      ++modes_selected;
    }
  }

  // Add prediction mode coding cost as the last thing. We don't want this
  // affecting the halving search.
  int lambda_cost = (int)(encoder_state->global->cur_lambda_cost_sqrt + 0.5);
  for (int mode_i = 0; mode_i < modes_selected; ++mode_i) {
    costs[mode_i] += lambda_cost * luma_mode_bits(encoder_state, modes[mode_i], intra_preds);
  }

  sort_modes(modes, costs, modes_selected);
  return modes_selected;
}


static void search_intra_rdo(encoder_state * const encoder_state, 
                             int x_px, int y_px, int depth,
                             pixel *orig, int32_t origstride,
                             pixel *rec, int16_t recstride,
                             int8_t *intra_preds,
                             int modes_to_check,
                             int8_t modes[35], double costs[35],
                             lcu_t *lcu)
{
  const int tr_depth = CLIP(1, MAX_PU_DEPTH, depth + encoder_state->encoder_control->tr_depth_intra);
  const int width = LCU_WIDTH >> depth;

  pixel pred[LCU_WIDTH * LCU_WIDTH + 1];
  pixel orig_block[LCU_WIDTH * LCU_WIDTH + 1];
  int rdo_mode;
  int pred_mode;

  pixel rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];
  pixel *ref[2] = {rec, &rec_filtered_temp[recstride + 1]};

  // Generate filtered reference pixels.
  {
    int x, y;
    for (y = -1; y < recstride; y++) {
      ref[1][y*recstride - 1] = rec[y*recstride - 1];
    }
    for (x = 0; x < recstride; x++) {
      ref[1][x - recstride] = rec[x - recstride];
    }
    intra_filter(ref[1], recstride, width, 0);
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
    int rdo_bitcost = luma_mode_bits(encoder_state, modes[rdo_mode], intra_preds);
    costs[rdo_mode] = rdo_bitcost * (int)(encoder_state->global->cur_lambda_cost + 0.5);

    if (0 && tr_depth == depth) {
      // The reconstruction is calculated again here, it could be saved from before..
      intra_get_pred(encoder_state->encoder_control, ref, recstride, pred, width, modes[rdo_mode], 0);
      costs[rdo_mode] += rdo_cost_intra(encoder_state, pred, orig_block, width, modes[rdo_mode], width == 4 ? 1 : 0);
    } else {
      // Perform transform split search and save mode RD cost for the best one.
      cu_info pred_cu;
      pred_cu.depth = depth;
      pred_cu.type = CU_INTRA;
      pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N);
      pred_cu.intra[0].mode = modes[rdo_mode];
      pred_cu.intra[1].mode = modes[rdo_mode];
      pred_cu.intra[2].mode = modes[rdo_mode];
      pred_cu.intra[3].mode = modes[rdo_mode];
      pred_cu.intra[0].mode_chroma = modes[rdo_mode];
      memset(&pred_cu.cbf, 0, sizeof(pred_cu.cbf));

      // Reset transform split data in lcu.cu for this area.
      lcu_set_trdepth(lcu, x_px, y_px, depth, depth);

      double mode_cost = search_intra_trdepth(encoder_state, x_px, y_px, depth, tr_depth, modes[rdo_mode], MAX_INT, &pred_cu, lcu);
      costs[rdo_mode] += mode_cost;
    }
  }

  sort_modes(modes, costs, modes_to_check);

  if (tr_depth != depth) {
    cu_info pred_cu;
    pred_cu.depth = depth;
    pred_cu.type = CU_INTRA;
    pred_cu.part_size = ((depth == MAX_PU_DEPTH) ? SIZE_NxN : SIZE_2Nx2N);
    pred_cu.intra[0].mode = modes[0];
    pred_cu.intra[1].mode = modes[0];
    pred_cu.intra[2].mode = modes[0];
    pred_cu.intra[3].mode = modes[0];
    pred_cu.intra[0].mode_chroma = modes[0];
    memset(&pred_cu.cbf, 0, sizeof(pred_cu.cbf));
    search_intra_trdepth(encoder_state, x_px, y_px, depth, tr_depth, modes[0], MAX_INT, &pred_cu, lcu);
  }
}


/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
static double search_cu_intra(encoder_state * const encoder_state,
                           const int x_px, const int y_px,
                           const int depth, lcu_t *lcu)
{
  const videoframe * const frame = encoder_state->tile->frame;
  const vector2d lcu_px = { x_px & 0x3f, y_px & 0x3f };
  const vector2d lcu_cu = { lcu_px.x >> 3, lcu_px.y >> 3 };
  const int8_t cu_width = (LCU_WIDTH >> (depth));
  const int cu_index = LCU_CU_OFFSET + lcu_cu.x + lcu_cu.y * LCU_T_CU_WIDTH;

  cu_info *cur_cu = &lcu->cu[cu_index];

  pixel rec_buffer[(LCU_WIDTH * 2 + 1) * (LCU_WIDTH * 2 + 1)];
  pixel *cu_in_rec_buffer = &rec_buffer[cu_width * 2 + 8 + 1];

  int8_t candidate_modes[3];

  cu_info *left_cu = 0;
  cu_info *above_cu = 0;

  if ((x_px >> 3) > 0) {
    left_cu = &lcu->cu[cu_index - 1];
  }
  // Don't take the above CU across the LCU boundary.
  if ((y_px >> 3) > 0 && lcu_cu.y != 0) {
    above_cu = &lcu->cu[cu_index - LCU_T_CU_WIDTH];
  }

  // Get intra predictors
  intra_get_dir_luma_predictor(x_px, y_px, candidate_modes, cur_cu, left_cu, above_cu);

  if (depth > 0) {
  // Build reconstructed block to use in prediction with extrapolated borders
  intra_build_reference_border(encoder_state->encoder_control, x_px, y_px, cu_width * 2 + 8,
                               rec_buffer, cu_width * 2 + 8, 0,
                               frame->width,
                               frame->height,
                               lcu);
  }

  int8_t modes[35];
  double costs[35];

  // Find best intra mode for 2Nx2N.
  {
    pixel *ref_pixels = &lcu->ref.y[lcu_px.x + lcu_px.y * LCU_WIDTH];
    unsigned pu_index = PU_INDEX(x_px >> 2, y_px >> 2);

    int8_t number_of_modes;
    bool skip_rough_search = (depth == 0 || encoder_state->encoder_control->rdo >= 3);
    if (!skip_rough_search) {
      number_of_modes = search_intra_rough(encoder_state,
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

    if (encoder_state->encoder_control->rdo >= 2) {
      int number_of_modes_to_search = (cu_width <= 8) ? 8 : 3;
      if (encoder_state->encoder_control->rdo == 3) {
        number_of_modes_to_search = 35;
      }
      int num_modes_to_check = MIN(number_of_modes, number_of_modes_to_search);
      search_intra_rdo(encoder_state, 
                       x_px, y_px, depth,
                       ref_pixels, LCU_WIDTH,
                       cu_in_rec_buffer, cu_width * 2 + 8,
                       candidate_modes,
                       num_modes_to_check,
                       modes, costs, lcu);
    }

    cur_cu->intra[pu_index].mode = modes[0];
  }

  return costs[0];
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
static double search_cu(encoder_state * const encoder_state, int x, int y, int depth, lcu_t work_tree[MAX_PU_DEPTH])
{
  const videoframe * const frame = encoder_state->tile->frame;
  int cu_width = LCU_WIDTH >> depth;
  double cost = MAX_INT;
  cu_info *cur_cu;
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

    if (encoder_state->global->slicetype != SLICE_I &&
        depth >= MIN_INTER_SEARCH_DEPTH &&
        depth <= MAX_INTER_SEARCH_DEPTH)
    {
      int mode_cost = search_cu_inter(encoder_state, x, y, depth, &work_tree[depth]);
      if (mode_cost < cost) {
        cost = mode_cost;
        cur_cu->type = CU_INTER;
      }
    }

    if (depth >= MIN_INTRA_SEARCH_DEPTH &&
        depth <= MAX_INTRA_SEARCH_DEPTH)
    {
      double mode_cost = search_cu_intra(encoder_state, x, y, depth, &work_tree[depth]);
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
                         100,
                         cur_cu->part_size);
      intra_recon_lcu_luma(encoder_state, x, y, depth, intra_mode, NULL, &work_tree[depth]);

      if (PU_INDEX(x >> 2, y >> 2) == 0) {
        int8_t intra_mode_chroma = intra_mode;
        if (encoder_state->encoder_control->rdo >= 2) {
          intra_mode_chroma = search_intra_chroma(encoder_state, x, y, depth, intra_mode, &work_tree[depth]);
        }
        lcu_set_intra_mode(&work_tree[depth], x, y, depth,
                           intra_mode, intra_mode_chroma,
                           cur_cu->part_size);
        intra_recon_lcu_chroma(encoder_state, x, y, depth, intra_mode_chroma, NULL, &work_tree[depth]);
      }
    } else if (cur_cu->type == CU_INTER) {
      int cbf;
      inter_recon_lcu(encoder_state, encoder_state->global->ref->images[cur_cu->inter.mv_ref], x, y, LCU_WIDTH>>depth, cur_cu->inter.mv, &work_tree[depth]);
      quantize_lcu_luma_residual(encoder_state, x, y, depth, NULL, &work_tree[depth]);
      quantize_lcu_chroma_residual(encoder_state, x, y, depth, NULL, &work_tree[depth]);

      cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);

      if(cur_cu->merged && !cbf) {
        cur_cu->merged = 0;
        cur_cu->skipped = 1;
        // Selecting skip reduces bits needed to code the CU
        cur_cu->inter.bitcost--;
      }
      lcu_set_inter(&work_tree[depth], x, y, depth, cur_cu);
      lcu_set_coeff(&work_tree[depth], x, y, depth, cur_cu);
    }
  }
  if (cur_cu->type == CU_INTRA || cur_cu->type == CU_INTER) {
    cost = cu_rd_cost_luma(encoder_state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    cost += cu_rd_cost_chroma(encoder_state, x_local, y_local, depth, cur_cu, &work_tree[depth]);
    
    double mode_bits;
    // Bitcost
    if (cur_cu->type == CU_INTER) {
      mode_bits = cur_cu->inter.bitcost;
    } else {
      int8_t candidate_modes[3];
      {
        lcu_t *lcu = &work_tree[depth];
        const vector2d lcu_px = { x & 0x3f, y & 0x3f };
        const vector2d lcu_cu = { lcu_px.x >> 3, lcu_px.y >> 3 };
        const cu_info *left_cu = ((x >> 3) ? &cur_cu[-1] : NULL);
        const cu_info *above_cu = ((lcu_cu.y) ? &cur_cu[-LCU_T_CU_WIDTH] : NULL);
        intra_get_dir_luma_predictor(x, y, candidate_modes, cur_cu, left_cu, above_cu);
      }

      mode_bits = luma_mode_bits(encoder_state, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode, candidate_modes);
      if (PU_INDEX(x >> 2, y >> 2) == 0) {
        mode_bits += chroma_mode_bits(encoder_state, cur_cu->intra[0].mode_chroma, cur_cu->intra[PU_INDEX(x >> 2, y >> 2)].mode);
      }
    }
    cost += mode_bits * encoder_state->global->cur_lambda_cost;
  }
  
  // Recursively split all the way to max search depth.
  if (depth < MAX_INTRA_SEARCH_DEPTH || (depth < MAX_INTER_SEARCH_DEPTH && encoder_state->global->slicetype != SLICE_I)) {
    int half_cu = cu_width / 2;
    // Using Cost = lambda * 9 to compensate on the price of the split
    double split_cost = encoder_state->global->cur_lambda_cost * CU_COST;
    int cbf = cbf_is_set(cur_cu->cbf.y, depth) || cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth);
        
    if (depth < MAX_DEPTH) {
      vector2d lcu_cu = { x_local / 8, y_local / 8 };
      cu_info *cu_array = &(&work_tree[depth])->cu[LCU_CU_OFFSET];
      bool condA = x >= 8 && cu_array[(lcu_cu.x - 1) * lcu_cu.y * LCU_T_CU_WIDTH].depth > depth;
      bool condL = y >= 8 && cu_array[lcu_cu.x * (lcu_cu.y - 1) * LCU_T_CU_WIDTH].depth > depth;
      uint8_t split_model = condA + condL;

      const cabac_ctx *ctx = &(encoder_state->cabac.ctx.split_flag_model[split_model]);
      cost += CTX_ENTROPY_FBITS(ctx, 0);
      split_cost += CTX_ENTROPY_FBITS(ctx, 1);
    }

    if (cur_cu->type == CU_INTRA && depth == MAX_DEPTH) {
      const cabac_ctx *ctx = &(encoder_state->cabac.ctx.part_size_model[0]);
      cost += CTX_ENTROPY_FBITS(ctx, 1);  // 2Nx2N
      split_cost += CTX_ENTROPY_FBITS(ctx, 0);  // NxN
    }

    // If skip mode was selected for the block, skip further search.
    // Skip mode means there's no coefficients in the block, so splitting
    // might not give any better results but takes more time to do.
    if (cur_cu->type == CU_NOTSET || cbf || FULL_CU_SPLIT_SEARCH) {
      split_cost += search_cu(encoder_state, x,           y,           depth + 1, work_tree);
      split_cost += search_cu(encoder_state, x + half_cu, y,           depth + 1, work_tree);
      split_cost += search_cu(encoder_state, x,           y + half_cu, depth + 1, work_tree);
      split_cost += search_cu(encoder_state, x + half_cu, y + half_cu, depth + 1, work_tree);
    } else {
      split_cost = INT_MAX;
    }
    if (split_cost < cost) {
      // Copy split modes to this depth.
      cost = split_cost;
      work_tree_copy_up(x, y, depth, work_tree);
#if _DEBUG
      debug_split = 1;
#endif
    } else {
      // Copy this CU's mode all the way down for use in adjacent CUs mode
      // search.
      work_tree_copy_down(x, y, depth, work_tree);
    }
  }
  
  PERFORMANCE_MEASURE_END(_DEBUG_PERF_SEARCH_CU, encoder_state->encoder_control->threadqueue, "type=search_cu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d,depth=%d,split=%d,cur_cu_is_intra=%d", encoder_state->global->frame, encoder_state->tile->id, encoder_state->slice->id, 
                          (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + x,
                          (encoder_state->tile->lcu_offset_x * LCU_WIDTH) + x + (LCU_WIDTH >> depth), 
                          (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + y,
                          (encoder_state->tile->lcu_offset_y * LCU_WIDTH) + y + (LCU_WIDTH >> depth), 
                          depth, debug_split, (cur_cu->type==CU_INTRA)?1:0);

  return cost;
}


/**
 * Initialize lcu_t for search.
 * - Copy reference CUs.
 * - Copy reference pixels from neighbouring LCUs.
 * - Copy reference pixels from this LCU.
 */
static void init_lcu_t(const encoder_state * const encoder_state, const int x, const int y, lcu_t *lcu, const yuv_t *hor_buf, const yuv_t *ver_buf)
{
  const videoframe * const frame = encoder_state->tile->frame;
  
  // Copy reference cu_info structs from neighbouring LCUs.
  {
    const int x_cu = x >> MAX_DEPTH;
    const int y_cu = y >> MAX_DEPTH;

    // Use top-left sub-cu of LCU as pointer to lcu->cu array to make things
    // simpler.
    cu_info *lcu_cu = &lcu->cu[LCU_CU_OFFSET];

    // Copy top CU row.
    if (y_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info *from_cu = videoframe_get_cu_const(frame, x_cu + i, y_cu - 1);
        cu_info *to_cu = &lcu_cu[i - LCU_T_CU_WIDTH];
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy left CU column.
    if (x_cu > 0) {
      int i;
      for (i = 0; i < LCU_CU_WIDTH; ++i) {
        const cu_info *from_cu = videoframe_get_cu_const(frame, x_cu - 1, y_cu + i);
        cu_info *to_cu = &lcu_cu[-1 + i * LCU_T_CU_WIDTH];
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
    // Copy top-left CU.
    if (x_cu > 0 && y_cu > 0) {
      const cu_info *from_cu = videoframe_get_cu_const(frame, x_cu - 1, y_cu - 1);
      cu_info *to_cu = &lcu_cu[-1 - LCU_T_CU_WIDTH];
      memcpy(to_cu, from_cu, sizeof(*to_cu));
    }

    // Copy top-right CU.
    if (y_cu > 0 && x + LCU_WIDTH < frame->width) {
      const cu_info *from_cu = videoframe_get_cu_const(frame, x_cu + LCU_CU_WIDTH, y_cu - 1);
      cu_info *to_cu = &lcu->cu[LCU_T_CU_WIDTH*LCU_T_CU_WIDTH];
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
    const videoframe * const frame = encoder_state->tile->frame;
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
static void copy_lcu_to_cu_data(const encoder_state * const encoder_state, int x_px, int y_px, const lcu_t *lcu)
{
  // Copy non-reference CUs to picture.
  {
    const int x_cu = x_px >> MAX_DEPTH;
    const int y_cu = y_px >> MAX_DEPTH;
    videoframe * const frame = encoder_state->tile->frame;

    // Use top-left sub-cu of LCU as pointer to lcu->cu array to make things
    // simpler.
    const cu_info *const lcu_cu = &lcu->cu[LCU_CU_OFFSET];

    int x, y;
    for (y = 0; y < LCU_CU_WIDTH; ++y) {
      for (x = 0; x < LCU_CU_WIDTH; ++x) {
        const cu_info *from_cu = &lcu_cu[x + y * LCU_T_CU_WIDTH];
        cu_info *to_cu = videoframe_get_cu(frame, x_cu + x, y_cu + y);
        memcpy(to_cu, from_cu, sizeof(*to_cu));
      }
    }
  }

  // Copy pixels to picture.
  {
    videoframe * const pic = encoder_state->tile->frame;
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
void search_lcu(encoder_state * const encoder_state, const int x, const int y, const yuv_t * const hor_buf, const yuv_t * const ver_buf)
{
  lcu_t work_tree[MAX_PU_DEPTH + 1];
  int depth;
  // Initialize work tree.
  for (depth = 0; depth <= MAX_PU_DEPTH; ++depth) {
    memset(&work_tree[depth], 0, sizeof(work_tree[depth]));
    init_lcu_t(encoder_state, x, y, &work_tree[depth], hor_buf, ver_buf);
  }

  // Start search from depth 0.
  search_cu(encoder_state, x, y, 0, work_tree);

  copy_lcu_to_cu_data(encoder_state, x, y, &work_tree[0]);
}
