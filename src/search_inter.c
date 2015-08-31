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

#include "search_inter.h"

#include <stdlib.h>

#include "inter.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-ipol.h"
#include "rdo.h"


// Temporarily for debugging.
#define SEARCH_MV_FULL_RADIUS 0


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


static uint32_t get_mvd_coding_cost(vector2d_t *mvd, cabac_data_t* cabac)
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
    cand1_cost = get_mvd_coding_cost(&mvd_temp1, NULL);

    mvd_temp2.x = x - mv_cand[1][0];
    mvd_temp2.y = y - mv_cand[1][1];
    cand2_cost = get_mvd_coding_cost(&mvd_temp2, NULL);

    // Select candidate 1 if it has lower cost
    if (cand2_cost < cand1_cost) {
      cur_mv_cand = 1;
    }
    temp_bitcost += cur_mv_cand ? cand2_cost : cand1_cost;
  }
  *bitcost = temp_bitcost;
  return temp_bitcost*(int32_t)(state->global->cur_lambda_cost_sqrt+0.5);
}


unsigned kvz_tz_pattern_search(const encoder_state_t * const state, const kvz_picture *pic, const kvz_picture *ref, unsigned pattern_type,
                           const vector2d_t *orig, const int iDist, vector2d_t *mv, unsigned best_cost, int *best_dist,
                           int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS], int16_t num_cand, int32_t ref_idx, uint32_t *best_bitcost,
                           int block_width, int max_px_below_lcu)
{
  int n_points;
  int best_index = -1;
  int i;
  
  vector2d_t mv_best = { 0, 0 };


  int(*calc_mvd)(const encoder_state_t * const, int, int, int,
    int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
    int16_t, int32_t, uint32_t *) = calc_mvd_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    calc_mvd = kvz_calc_mvd_cost_cabac;
  }

  assert(pattern_type < 4);

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
      PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);
      cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                            (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + current->x,
                            (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + current->y,
                            block_width, block_width, max_px_below_lcu);
      cost += calc_mvd(state, mv->x + current->x, mv->y + current->y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

      PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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


unsigned kvz_tz_raster_search(const encoder_state_t * const state, const kvz_picture *pic, const kvz_picture *ref,
                          const vector2d_t *orig, vector2d_t *mv, unsigned best_cost,
                          int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS], int16_t num_cand, int32_t ref_idx, uint32_t *best_bitcost,
                          int block_width, int iSearchRange, int iRaster, int max_px_below_lcu)
{
  int i;
  int k;

  vector2d_t mv_best = { 0, 0 };

  int(*calc_mvd)(const encoder_state_t * const, int, int, int,
    int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
    int16_t, int32_t, uint32_t *) = calc_mvd_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    calc_mvd = kvz_calc_mvd_cost_cabac;
  }
  
  //compute SAD values for every point in the iRaster downsampled version of the current search area
  for (i = iSearchRange; i >= -iSearchRange; i -= iRaster)
  {
    for (k = -iSearchRange; k <= iSearchRange; k += iRaster)
    {
      vector2d_t current = { k, i };
      unsigned cost;
      uint32_t bitcost;

      {
        PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);
        cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
          (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv->x + k,
          (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv->y + i,
          block_width, block_width, max_px_below_lcu);
        cost += calc_mvd(state, mv->x + k, mv->y + i, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

        PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
                          const kvz_picture *pic, const kvz_picture *ref,
                          const vector2d_t *orig, vector2d_t *mv_in_out,
                          int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                          int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{

  //TZ parameters
  const int iSearchRange = 96;  // search range for each stage
  const int iRaster = 5;  // search distance limit and downsampling factor for step 3                   
  const unsigned step2_type = 0;  // search patterns for steps 2 and 4
  const unsigned step4_type = 0;
  const bool bRasterRefinementEnable = true;  // enable step 4 mode 1
  const bool bStarRefinementEnable = false;   // enable step 4 mode 2 (only one mode will be executed)
  
  const int block_width = CU_WIDTH_FROM_DEPTH(depth);

  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };

  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0;
  int iDist;
  int best_dist = 0;
  unsigned best_index = num_cand;
  int max_px_below_lcu = -1;

  int(*calc_mvd)(const encoder_state_t * const, int, int, int,
    int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
    int16_t, int32_t, uint32_t *) = calc_mvd_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    calc_mvd = kvz_calc_mvd_cost_cabac;
  }

  if (state->encoder_control->owf) {
    max_px_below_lcu = LCU_WIDTH;
    if (state->encoder_control->fme_level > 0) {
      // Fractional motion estimation can change the mv by at most 1 pixel.
      max_px_below_lcu -= 1;
    }
    if (state->encoder_control->deblock_enable) {
      // Strong deblock filter modifies 3 pixels.
      max_px_below_lcu -= 3;
    }
  }

  //step 1, compare (0,0) vector to predicted vectors
  
  // Check whatever input vector we got, unless its (0, 0) which will be checked later.
  if (mv.x || mv.y) 
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);

    best_cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                                        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                        block_width, block_width, max_px_below_lcu);
    best_cost += calc_mvd(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost);

    PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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

    PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);

	  uint32_t bitcost;
    unsigned cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                                   (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                   (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                   block_width, block_width, max_px_below_lcu);
    cost += calc_mvd(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

    PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
    best_cost = kvz_tz_pattern_search(state, pic, ref, step2_type, orig, iDist, &mv, best_cost, &best_dist,
                                  mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_px_below_lcu);
  }

  //step 3, raster scan
  if (best_dist > iRaster)
  {
    best_dist = iRaster;

    best_cost = kvz_tz_raster_search(state, pic, ref, orig, &mv, best_cost, mv_cand, merge_cand, 
                                 num_cand, ref_idx, &best_bitcost, block_width, iSearchRange, iRaster, max_px_below_lcu);
  }

  //step 4

  //raster refinement
  if (bRasterRefinementEnable && best_dist > 0)
  {
    iDist = best_dist >> 1;
    while (iDist > 0)
    {
      best_cost = kvz_tz_pattern_search(state, pic, ref, step4_type, orig, iDist, &mv, best_cost, &best_dist,
                                   mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_px_below_lcu);

      iDist = iDist >> 1;
    }
  }

  //star refinement (repeat step 2 for the current starting point)
  if (bStarRefinementEnable && best_dist > 0)
  {
    for (iDist = 1; iDist <= iSearchRange; iDist *= 2)
    {
      best_cost = kvz_tz_pattern_search(state, pic, ref, step4_type, orig, iDist, &mv, best_cost, &best_dist,
                                   mv_cand, merge_cand, num_cand, ref_idx, &best_bitcost, block_width, max_px_below_lcu);
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
                               const kvz_picture *pic, const kvz_picture *ref,
                               const vector2d_t *orig, vector2d_t *mv_in_out,
                               int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                               int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{
  // The start of the hexagonal pattern has been repeated at the end so that
  // the indices between 1-6 can be used as the start of a 3-point list of new
  // points to search.
  //   6--1,7
  //  /     \    =)
  // 5   0  2,8
  //  \     /
  //   4---3
  static const vector2d_t large_hexbs[10] = {
      { 0, 0 },
      { 1, -2 }, { 2, 0 }, { 1, 2 }, { -1, 2 }, { -2, 0 }, { -1, -2 },
      { 1, -2 }, { 2, 0 }
  };
  // This is used as the last step of the hexagon search.
  //   1
  // 2 0 3
  //   4
  static const vector2d_t small_hexbs[5] = {
      { 0, 0 },
      { 0, -1 }, { -1, 0 }, { 1, 0 }, { 0, 1 }
  };

  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0, bitcost;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.
  int max_px_below_lcu = -1;

  int (*calc_mvd)(const encoder_state_t * const, int, int, int,
    int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
    int16_t, int32_t, uint32_t *) = calc_mvd_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    calc_mvd = kvz_calc_mvd_cost_cabac;
  }

  
  if (state->encoder_control->owf) {
    max_px_below_lcu = LCU_WIDTH;
    if (state->encoder_control->fme_level > 0) {
      // Fractional motion estimation can change the mv by at most 1 pixel.
      max_px_below_lcu -= 1;
    }
    if (state->encoder_control->deblock_enable) {
      // Strong deblock filter modifies 3 pixels.
      max_px_below_lcu -= 3;
    }
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
    PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);

    best_cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                                        (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                        (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                        block_width, block_width, max_px_below_lcu);
    best_cost += calc_mvd(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);
    best_bitcost = bitcost;
    best_index = num_cand; 

    PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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

    PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);

    unsigned cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                                   (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x,
                                   (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y,
                                   block_width, block_width, max_px_below_lcu);
    cost += calc_mvd(state, mv.x, mv.y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

    PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
      PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);
      cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                             (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + pattern->x, 
                             (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + pattern->y,
                             block_width, block_width, max_px_below_lcu);
      cost += calc_mvd(state, mv.x + pattern->x, mv.y + pattern->y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

      PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
        PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);
        cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                               (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                               (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                               block_width, block_width, max_px_below_lcu);
        cost += calc_mvd(state, mv.x + offset->x, mv.y + offset->y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);
        PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=large_hexbs_iterative,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
      PERFORMANCE_MEASURE_START(KVZ_PERF_SEARCHPX);
      cost = kvz_image_calc_sad(pic, ref, orig->x, orig->y,
                             (state->tile->lcu_offset_x * LCU_WIDTH) + orig->x + mv.x + offset->x,
                             (state->tile->lcu_offset_y * LCU_WIDTH) + orig->y + mv.y + offset->y,
                             block_width, block_width, max_px_below_lcu);
      cost += calc_mvd(state, mv.x + offset->x, mv.y + offset->y, 2, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);
      PERFORMANCE_MEASURE_END(KVZ_PERF_SEARCHPX, state->encoder_control->threadqueue, "type=sad,step=small_hexbs,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d,ref_px_x=%d-%d,ref_px_y=%d-%d", state->global->frame, state->tile->id, orig->x, orig->x + block_width, orig->y, orig->y + block_width,
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
                            const kvz_picture *pic, const kvz_picture *ref,
                            const vector2d_t *orig, vector2d_t *mv_in_out,
                            int16_t mv_cand[2][2], inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                            int16_t num_cand, int32_t ref_idx, uint32_t *bitcost_out)
{
  // Map indexes to relative coordinates in the following way:
  // 6 7 8
  // 3 4 5
  // 0 1 2
  static const vector2d_t square[9] = {
      { -1, 1 },  { 0, 1 },  { 1, 1 },
      { -1, 0 },  { 0, 0 },  { 1, 0 },
      { -1, -1 }, { 0, -1 }, { 1, -1 }
  };

  //Set mv to halfpel precision
  vector2d_t mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = UINT32_MAX;
  uint32_t best_bitcost = 0, bitcost;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.

  unsigned cost = 0;

  cost_pixel_nxn_func *satd = kvz_pixels_get_satd_func(block_width);

  vector2d_t halfpel_offset;

  #define FILTER_SIZE 8
  #define HALF_FILTER (FILTER_SIZE>>1)

  kvz_extended_block src = { 0, 0, 0 };

  //destination buffer for interpolation
  int dst_stride = (block_width+1)*4;
  kvz_pixel dst[(LCU_WIDTH+1) * (LCU_WIDTH+1) * 16];
  kvz_pixel* dst_off = &dst[dst_stride*4+4];

  int(*calc_mvd)(const encoder_state_t * const, int, int, int,
    int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
    int16_t, int32_t, uint32_t *) = calc_mvd_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    calc_mvd = kvz_calc_mvd_cost_cabac;
  }

  kvz_get_extended_block(orig->x, orig->y, mv.x-1, mv.y-1,
                state->tile->lcu_offset_x * LCU_WIDTH,
                state->tile->lcu_offset_y * LCU_WIDTH,
                ref->y, ref->width, ref->height, FILTER_SIZE, block_width+1, block_width+1, &src);

  kvz_filter_inter_quarterpel_luma(state->encoder_control, src.orig_topleft, src.stride, block_width+1,
      block_width+1, dst, dst_stride, 1, 1);

  if (src.malloc_used) free(src.buffer);

  //Set mv to half-pixel precision
  mv.x <<= 1;
  mv.y <<= 1;

  kvz_pixel tmp_filtered[LCU_WIDTH*LCU_WIDTH];
  kvz_pixel tmp_pic[LCU_WIDTH*LCU_WIDTH];
  kvz_pixels_blit(pic->y + orig->y*pic->width + orig->x, tmp_pic, block_width, block_width, pic->stride, block_width);

  // Search halfpel positions around best integer mv
  for (i = 0; i < 9; ++i) {
    const vector2d_t *pattern = &square[i];

    int y,x;
    for(y = 0; y < block_width; ++y) {
      int dst_y = y*4+pattern->y*2;
      for(x = 0; x < block_width; ++x) {
        int dst_x = x*4+pattern->x*2;
        tmp_filtered[y*block_width+x] = dst_off[dst_y*dst_stride+dst_x];
      }
    }

    cost = satd(tmp_pic,tmp_filtered);

    cost += calc_mvd(state, mv.x + pattern->x, mv.y + pattern->y, 1, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

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

    int y,x;
    for(y = 0; y < block_width; ++y) {
      int dst_y = y*4+halfpel_offset.y+pattern->y;
      for(x = 0; x < block_width; ++x) {
        int dst_x = x*4+halfpel_offset.x+pattern->x;
        tmp_filtered[y*block_width+x] = dst_off[dst_y*dst_stride+dst_x];
      }
    }

    cost = satd(tmp_pic,tmp_filtered);

    cost += calc_mvd(state, mv.x + pattern->x, mv.y + pattern->y, 0, mv_cand, merge_cand, num_cand, ref_idx, &bitcost);

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
 * \brief Perform inter search for a single reference frame.
 */
static void search_cu_inter_ref(const encoder_state_t * const state,
                                int x, int y, int depth,
                                lcu_t *lcu, cu_info_t *cur_cu,
                                int16_t mv_cand[2][2],
                                inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                                int16_t num_cand,
                                unsigned ref_idx,
                                uint32_t(*get_mvd_cost)(vector2d_t *, cabac_data_t*))
{
  const int x_cu = x >> 3;
  const int y_cu = y >> 3;
  const videoframe_t * const frame = state->tile->frame;
  kvz_picture *ref_image = state->global->ref->images[ref_idx];
  uint32_t temp_bitcost = 0;
  uint32_t temp_cost = 0;
  vector2d_t orig, mvd;
  int32_t merged = 0;
  uint8_t cu_mv_cand = 0;
  int8_t merge_idx = 0;
  int8_t ref_list = state->global->refmap[ref_idx].list-1;
  int8_t temp_ref_idx = cur_cu->inter.mv_ref[ref_list];
  orig.x = x_cu * CU_MIN_SIZE_PIXELS;
  orig.y = y_cu * CU_MIN_SIZE_PIXELS;
  // Get MV candidates
  cur_cu->inter.mv_ref[ref_list] = ref_idx;
  kvz_inter_get_mv_cand(state, x, y, LCU_WIDTH >> depth, LCU_WIDTH >> depth, mv_cand, cur_cu, lcu, ref_list);
  cur_cu->inter.mv_ref[ref_list] = temp_ref_idx;


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
    case KVZ_IME_TZ:
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
        (uint32_t)merge_cand[merge_idx].ref[merge_cand[merge_idx].dir - 1] == ref_idx) {
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
    cand1_cost = get_mvd_cost(&mvd_temp1, (cabac_data_t*)&state->cabac);

    mvd_temp2.x = mv.x - mv_cand[1][0];
    mvd_temp2.y = mv.y - mv_cand[1][1];
    cand2_cost = get_mvd_cost(&mvd_temp2, (cabac_data_t*)&state->cabac);

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
    cur_cu->inter.mv_cand[ref_list] = cu_mv_cand;
  }
}


/**
 * Update lcu to have best modes at this depth.
 * \return Cost of best mode.
 */
int kvz_search_cu_inter(const encoder_state_t * const state, int x, int y, int depth, lcu_t *lcu)
{
  const videoframe_t * const frame = state->tile->frame;
  uint32_t ref_idx = 0;
  int x_local = SUB_SCU(x);
  int y_local = SUB_SCU(y);
  cu_info_t *cur_cu = LCU_GET_CU(lcu, x_local >> 3, y_local >> 3);

  int16_t mv_cand[2][2];
  // Search for merge mode candidate
  inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS];
  // Get list of candidates
  int16_t num_cand = kvz_inter_get_merge_cand(state, x, y, depth, merge_cand, lcu);

  uint32_t(*get_mvd_cost)(vector2d_t *, cabac_data_t*) = get_mvd_coding_cost;
  if (state->encoder_control->cfg->mv_rdo) {
    get_mvd_cost = kvz_get_mvd_coding_cost_cabac;
  }

  int max_px_below_lcu = -1;
  if (state->encoder_control->owf) {
    max_px_below_lcu = LCU_WIDTH;
    if (state->encoder_control->fme_level > 0) {
      // Fractional motion estimation can change the mv by at most 1 pixel.
      max_px_below_lcu -= 1;
    }
    if (state->encoder_control->deblock_enable) {
      // Strong deblock filter modifies 3 pixels.
      max_px_below_lcu -= 3;
    }
  }

  // Default to candidate 0
  cur_cu->inter.mv_cand[0] = 0;
  cur_cu->inter.mv_cand[1] = 0;

  cur_cu->inter.cost = UINT_MAX;

  for (ref_idx = 0; ref_idx < state->global->ref->used_size; ref_idx++) {
    search_cu_inter_ref(state,
                        x, y, depth,
                        lcu, cur_cu,
                        mv_cand, merge_cand, num_cand,
                        ref_idx,
                        get_mvd_cost);
  }

  // Search bi-pred positions
  if (state->global->slicetype == KVZ_SLICE_B && state->encoder_control->cfg->bipred) {
    lcu_t *templcu = MALLOC(lcu_t, 1);
    cost_pixel_nxn_func *satd = kvz_pixels_get_satd_func(LCU_WIDTH >> depth);
    #define NUM_PRIORITY_LIST 12;
    static const uint8_t priorityList0[] = { 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3 };
    static const uint8_t priorityList1[] = { 1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2 };
    uint8_t cutoff = num_cand;


    int(*calc_mvd)(const encoder_state_t * const, int, int, int,
      int16_t[2][2], inter_merge_cand_t[MRG_MAX_NUM_CANDS],
      int16_t, int32_t, uint32_t *) = calc_mvd_cost;
    if (state->encoder_control->cfg->mv_rdo) {
      calc_mvd = kvz_calc_mvd_cost_cabac;
    }

    for (int32_t idx = 0; idx<cutoff*(cutoff - 1); idx++) {
      uint8_t i = priorityList0[idx];
      uint8_t j = priorityList1[idx];
      if (i >= num_cand || j >= num_cand) break;

      // Find one L0 and L1 candidate according to the priority list
      if ((merge_cand[i].dir & 0x1) && (merge_cand[j].dir & 0x2)) {
        if (merge_cand[i].ref[0] != merge_cand[j].ref[1] ||
          merge_cand[i].mv[0][0] != merge_cand[j].mv[1][0] ||
          merge_cand[i].mv[0][1] != merge_cand[j].mv[1][1]) {
          uint32_t bitcost[2];
          uint32_t cost = 0;
          int8_t cu_mv_cand = 0;
          int16_t mv[2][2];
          kvz_pixel tmp_block[64 * 64];
          kvz_pixel tmp_pic[64 * 64];
          // Force L0 and L1 references
          if (state->global->refmap[merge_cand[i].ref[0]].list == 2 || state->global->refmap[merge_cand[j].ref[1]].list == 1) continue;

          mv[0][0] = merge_cand[i].mv[0][0];
          mv[0][1] = merge_cand[i].mv[0][1];
          mv[1][0] = merge_cand[j].mv[1][0];
          mv[1][1] = merge_cand[j].mv[1][1];

          // Check boundaries when using owf to process multiple frames at the same time
          if (max_px_below_lcu >= 0) {
            // The following has been modified to fit testing two mv's, but it's
            // equivalent to: y + mv + block_height - next_lcu_row_px > max_px_below_lcu
            int next_lcu_row_px = ((y >> LOG2_LCU_WIDTH) + 1) << LOG2_LCU_WIDTH;
            int block_height = LCU_WIDTH >> depth;
            int max_mv = next_lcu_row_px - y - block_height + max_px_below_lcu;
            int ceil_mv_l0 = ((mv[0][1] + 3) >> 2);
            int ceil_mv_l1 = ((mv[1][1] + 3) >> 2);

            if (ceil_mv_l0 < max_mv || ceil_mv_l1 < max_mv) continue;
          }

          kvz_inter_recon_lcu_bipred(state, state->global->ref->images[merge_cand[i].ref[0]], state->global->ref->images[merge_cand[j].ref[1]], x, y, LCU_WIDTH >> depth, mv, templcu);

          for (int ypos = 0; ypos < LCU_WIDTH >> depth; ++ypos) {
            int dst_y = ypos*(LCU_WIDTH >> depth);
            for (int xpos = 0; xpos < (LCU_WIDTH >> depth); ++xpos) {
              tmp_block[dst_y + xpos] = templcu->rec.y[
                SUB_SCU(y + ypos) * LCU_WIDTH + SUB_SCU(x + xpos)];
              tmp_pic[dst_y + xpos] = frame->source->y[x + xpos + (y + ypos)*frame->source->width];
            }
          }

          cost = satd(tmp_pic, tmp_block);

          cost += calc_mvd(state, merge_cand[i].mv[0][0], merge_cand[i].mv[0][1], 0, mv_cand, merge_cand, 0, ref_idx, &bitcost[0]);
          cost += calc_mvd(state, merge_cand[i].mv[1][0], merge_cand[i].mv[1][1], 0, mv_cand, merge_cand, 0, ref_idx, &bitcost[1]);

          if (cost < cur_cu->inter.cost) {

            cur_cu->inter.mv_dir = 3;
            cur_cu->inter.mv_ref_coded[0] = state->global->refmap[merge_cand[i].ref[0]].idx;
            cur_cu->inter.mv_ref_coded[1] = state->global->refmap[merge_cand[j].ref[1]].idx;



            cur_cu->inter.mv_ref[0] = merge_cand[i].ref[0];
            cur_cu->inter.mv_ref[1] = merge_cand[j].ref[1];

            cur_cu->inter.mv[0][0] = merge_cand[i].mv[0][0];
            cur_cu->inter.mv[0][1] = merge_cand[i].mv[0][1];
            cur_cu->inter.mv[1][0] = merge_cand[j].mv[1][0];
            cur_cu->inter.mv[1][1] = merge_cand[j].mv[1][1];
            cur_cu->merged = 0;
                        
            // Check every candidate to find a match
            for(int merge_idx = 0; merge_idx < num_cand; merge_idx++) {
              if (
                  merge_cand[merge_idx].mv[0][0] == cur_cu->inter.mv[0][0] &&
                  merge_cand[merge_idx].mv[0][1] == cur_cu->inter.mv[0][1] &&     
                  merge_cand[merge_idx].mv[1][0] == cur_cu->inter.mv[1][0] &&
                  merge_cand[merge_idx].mv[1][1] == cur_cu->inter.mv[1][1] &&    
                  merge_cand[merge_idx].ref[0] == cur_cu->inter.mv_ref[0] && 
                  merge_cand[merge_idx].ref[1] == cur_cu->inter.mv_ref[1]) {
                cur_cu->merged = 1;
                cur_cu->merge_idx = merge_idx;
                break;
              }
            }

            // Each motion vector has its own candidate
            for (int reflist = 0; reflist < 2; reflist++) {
              cu_mv_cand = 0;
              kvz_inter_get_mv_cand(state, x, y, LCU_WIDTH >> depth, LCU_WIDTH >> depth, mv_cand, cur_cu, lcu, reflist);
              if ((mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1])) {
                vector2d_t mvd_temp1, mvd_temp2;
                int cand1_cost, cand2_cost;

                mvd_temp1.x = cur_cu->inter.mv[reflist][0] - mv_cand[0][0];
                mvd_temp1.y = cur_cu->inter.mv[reflist][1] - mv_cand[0][1];
                cand1_cost = get_mvd_cost(&mvd_temp1, (cabac_data_t*)&state->cabac);

                mvd_temp2.x = cur_cu->inter.mv[reflist][0] - mv_cand[1][0];
                mvd_temp2.y = cur_cu->inter.mv[reflist][1] - mv_cand[1][1];
                cand2_cost = get_mvd_cost(&mvd_temp2, (cabac_data_t*)&state->cabac);

                // Select candidate 1 if it has lower cost
                if (cand2_cost < cand1_cost) {
                  cu_mv_cand = 1;                  
                }
              }
              cur_cu->inter.mvd[reflist][0] = cur_cu->inter.mv[reflist][0] - mv_cand[cu_mv_cand][0];
              cur_cu->inter.mvd[reflist][1] = cur_cu->inter.mv[reflist][1] - mv_cand[cu_mv_cand][1];
              cur_cu->inter.mv_cand[reflist] = cu_mv_cand;
            }
            cur_cu->inter.cost = cost;
            cur_cu->inter.bitcost = bitcost[0] + bitcost[1] + cur_cu->inter.mv_dir - 1 + cur_cu->inter.mv_ref_coded[0] + cur_cu->inter.mv_ref_coded[1];
          }
        }
      }
    }
    FREE_POINTER(templcu);
  }

  return cur_cu->inter.cost;
}
