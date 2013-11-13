/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "search.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "debug.h"


// Temporarily for debugging.
#define USE_INTRA_IN_P 1
//#define RENDER_CU encoder->frame==2
#define RENDER_CU 0
#define SEARCH_MV_FULL_RADIUS 0

#define IN_FRAME(x, y, width, height, block_width, block_height) \
  ((x) >= 0 && (y) >= 0 \
  && (x) + (block_width) <= (width) \
  && (y) + (block_height) <= (height))

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

int calc_mvd_cost(int x, int y, const vector2d *pred)
{
  int cost = 0;

  // Get the absolute difference vector and count the bits.
  x = abs(abs(x) - abs(pred->x));
  y = abs(abs(y) - abs(pred->y));
  while (x >>= 1) {
    ++cost;
  }
  while (y >>= 1) {
    ++cost;
  }

  // I don't know what is a good cost function for this. It probably doesn't
  // have to aproximate the actual cost of encoding the vector, but it's a
  // place to start.

  // Add two for quarter pixel resolution and multiply by two for Exp-Golomb.
  return (cost ? (cost + 2) << 1 : 0);
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
unsigned hexagon_search(unsigned depth, 
                        const picture *pic, const picture *ref,
                        const vector2d *orig, vector2d *mv_in_out)
{
  vector2d mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = -1;
  unsigned i;
  unsigned best_index = 0; // Index of large_hexbs or finally small_hexbs.

  // Search the initial 7 points of the hexagon.
  for (i = 0; i < 7; ++i) {
    const vector2d *pattern = &large_hexbs[i];
    unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                             orig->x + mv.x + pattern->x, orig->y + mv.y + pattern->y,
                             block_width, block_width);
    cost += calc_mvd_cost(mv.x + pattern->x, orig->y + mv.y + pattern->y, mv_in_out);

    if (cost < best_cost) {
      best_cost = cost;
      best_index = i;
    }
  }

  // Try the 0,0 vector.
  if (!(mv.x == 0 && mv.y == 0)) {
    unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                             orig->x, orig->y,
                             block_width, block_width);
    cost += calc_mvd_cost(0, 0, mv_in_out);
    
    // If the 0,0 is better, redo the hexagon around that point.
    if (cost < best_cost) {
      best_cost = cost;
      best_index = 0;
      mv.x = 0;
      mv.y = 0;

      for (i = 1; i < 7; ++i) {
        const vector2d *pattern = &large_hexbs[i];
        unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                                 orig->x + pattern->x, 
                                 orig->y + pattern->y,
                                 block_width, block_width);
        cost += calc_mvd_cost(pattern->x, pattern->y, mv_in_out);

        if (cost < best_cost) {
          best_cost = cost;
          best_index = i;
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
      unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                               orig->x + mv.x + offset->x, 
                               orig->y + mv.y + offset->y,
                               block_width, block_width);
      cost += calc_mvd_cost(mv.x + offset->x, mv.y + offset->y, mv_in_out);

      if (cost < best_cost) {
        best_cost = cost;
        best_index = start + i;
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
    unsigned cost = calc_sad(pic, ref, orig->x, orig->y,
                             orig->x + mv.x + offset->x,
                             orig->y + mv.y + offset->y,
                             block_width, block_width);
    cost += calc_mvd_cost(mv.x + offset->x, mv.y + offset->y, mv_in_out);

    if (cost > 0 && cost < best_cost) {
      best_cost = cost;
      best_index = i;
    }
  }

  // Adjust the movement vector according to the final best match.
  mv.x += small_hexbs[best_index].x;
  mv.y += small_hexbs[best_index].y;

  // Return final movement vector in quarter-pixel precision.
  mv_in_out->x = mv.x << 2;
  mv_in_out->y = mv.y << 2;

  return best_cost;
}

unsigned search_mv_full(unsigned depth, 
                        const picture *pic, const picture *ref,
                        const vector2d *orig, vector2d *mv_in_out)
{
  vector2d mv = { mv_in_out->x >> 2, mv_in_out->y >> 2 };
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned best_cost = -1;
  int x, y;
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
      cost += calc_mvd_cost(x, y, mv_in_out);
      if (cost < best_cost) {
        best_cost = cost;
        mv.x = x;
        mv.y = y;
      }
    }
  }

  mv_in_out->x = mv.x << 2;
  mv_in_out->y = mv.y << 2;

  return best_cost;
}

/**
 * \brief
 */
void search_buildReferenceBorder(picture *pic, int32_t x_ctb, int32_t y_ctb,
                                 int16_t outwidth, pixel *dst, 
                                 int32_t dststride, int8_t chroma)
{
  int32_t left_col; // left column iterator
  pixel val;      // variable to store extrapolated value
  int32_t i;        // index iterator
  pixel dc_val = 1 << (g_bitdepth - 1); // default predictor value
  int32_t top_row;  // top row iterator
  int32_t src_width = (pic->width >> (chroma ? 1 : 0));   // source picture width
  int32_t src_height = (pic->height >> (chroma ? 1 : 0)); // source picture height
  pixel *src_pic = (!chroma) ? pic->y_data : ((chroma == 1) ? pic->u_data : pic->v_data); // input picture pointer
  int16_t scu_width = LCU_WIDTH >> (MAX_DEPTH + (chroma ? 1 : 0)); // Smallest Coding Unit width
  pixel *src_shifted = &src_pic[x_ctb * scu_width + (y_ctb * scu_width) * src_width]; // input picture pointer shifted to start from the left-top corner of the current block
  int32_t width_in_scu = pic->width_in_lcu << MAX_DEPTH; // picture width in SCU

  // Fill left column
  if (x_ctb) {
    // Loop SCU's
    for (left_col = 1; left_col < outwidth / scu_width; left_col++) {
      // If over the picture height or block not yet searched, stop
      if ((y_ctb + left_col) * scu_width >= src_height
          || pic->cu_array[MAX_DEPTH][x_ctb - 1 + (y_ctb + left_col) * width_in_scu].type == CU_NOTSET) {
        break;
      }
    }

    // Copy the pixels to output
    for (i = 0; i < left_col * scu_width - 1; i++) {
      dst[(i + 1) * dststride] = src_shifted[i * src_width - 1];
    }

    // if the loop was not completed, extrapolate the last pixel pushed to output
    if (left_col != outwidth / scu_width) {
      val = src_shifted[(left_col * scu_width - 1) * src_width - 1];
      for (i = (left_col * scu_width); i < outwidth; i++) {
        dst[i * dststride] = val;
      }
    }
  } else { // If left column not available, copy from toprow or use the default predictor
    val = y_ctb ? src_shifted[-src_width] : dc_val;
    for (i = 0; i < outwidth; i++) {
      dst[i * dststride] = val;
    }
  }

  if (y_ctb) {
    // Loop top SCU's
    for (top_row = 1; top_row < outwidth / scu_width; top_row++) {
      if ((x_ctb + top_row) * scu_width >= src_width
          || pic->cu_array[MAX_DEPTH][x_ctb + top_row + (y_ctb - 1) * width_in_scu].type
              == CU_NOTSET) {
        break;
      }
    }

    for (i = 0; i < top_row * scu_width - 1; i++) {
      dst[i + 1] = src_shifted[i - src_width];
    }

    if (top_row != outwidth / scu_width) {
      val = src_shifted[(top_row * scu_width) - src_width - 1];
      for (i = (top_row * scu_width); i < outwidth; i++) {
        dst[i] = val;
      }
    }
  } else {
    val = x_ctb ? src_shifted[-1] : dc_val;
    for (i = 1; i < outwidth; i++) {
      dst[i] = val;
    }
  }
  // Topleft corner
  dst[0] = (x_ctb && y_ctb) ? src_shifted[-src_width - 1] : dst[dststride];

}

/**
 * \brief
 */
void search_tree(encoder_control *encoder, 
                 uint16_t x_ctb, uint16_t y_ctb, uint8_t depth)
{
  uint8_t border_x = ((encoder->in.width) < (x_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_y = ((encoder->in.height) < (y_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_split_x = ((encoder->in.width) < ((x_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border_split_y = ((encoder->in.height) < ((y_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border = border_x | border_y; // are we in any border CU

  picture *cur_pic = encoder->in.cur_pic;
  cu_info *cur_cu = &cur_pic->cu_array[depth][x_ctb + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)];

  cur_cu->intra.cost = 0xffffffff;
  cur_cu->inter.cost = 0xffffffff;

  // Force split on border
  if (depth != MAX_DEPTH) {
    if (border) {
      uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
      search_tree(encoder, x_ctb, y_ctb, depth + 1);
      if (!border_x || border_split_x) {
        search_tree(encoder, x_ctb + change, y_ctb, depth + 1);
      }
      if (!border_y || border_split_y) {
        search_tree(encoder, x_ctb, y_ctb + change, depth + 1);
      }
      if (!border || (border_split_x && border_split_y)) {
        search_tree(encoder, x_ctb + change, y_ctb + change, depth + 1);
      }
      return;
    }
  }

  // INTER SEARCH
  if (cur_pic->slicetype != SLICE_I
      && depth >= MIN_INTER_SEARCH_DEPTH && depth <= MAX_INTER_SEARCH_DEPTH) {
    
    picture *ref_pic = encoder->ref->pics[0];
    unsigned width_in_scu = NO_SCU_IN_LCU(ref_pic->width_in_lcu);
    cu_info *ref_cu = &ref_pic->cu_array[MAX_DEPTH][y_ctb * width_in_scu + x_ctb];

    vector2d orig, mv;
    orig.x = x_ctb * CU_MIN_SIZE_PIXELS;
    orig.y = y_ctb * CU_MIN_SIZE_PIXELS;
    mv.x = 0;
    mv.y = 0;
    if (ref_cu->type == CU_INTER) {
      mv.x = ref_cu->inter.mv[0];
      mv.y = ref_cu->inter.mv[1];
    }

#if SEARCH_MV_FULL_RADIUS
    cur_cu->inter.cost = search_mv_full(depth, cur_pic, ref_pic, &orig, &mv);
#else
    cur_cu->inter.cost = hexagon_search(depth, cur_pic, ref_pic, &orig, &mv);
#endif
    
    cur_cu->inter.mv_dir = 1;
    cur_cu->inter.mv[0] = mv.x;
    cur_cu->inter.mv[1] = mv.y;
  }

  // INTRA SEARCH
  if (depth >= MIN_INTRA_SEARCH_DEPTH && depth <= MAX_INTRA_SEARCH_DEPTH
      && (encoder->in.cur_pic->slicetype == SLICE_I || USE_INTRA_IN_P)) {
    int x = 0, y = 0;
    pixel *base = &encoder->in.cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH))) * encoder->in.width];
    uint32_t width = LCU_WIDTH >> depth;

    // INTRAPREDICTION
    pixel pred[LCU_WIDTH * LCU_WIDTH + 1];
    pixel rec[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8)];
    pixel *recShift = &rec[(LCU_WIDTH >> (depth)) * 2 + 8 + 1];

    // Build reconstructed block to use in prediction with extrapolated borders
    search_buildReferenceBorder(encoder->in.cur_pic, x_ctb, y_ctb,
        (LCU_WIDTH >> (depth)) * 2 + 8, rec, (LCU_WIDTH >> (depth)) * 2 + 8, 0);
    cur_cu->intra.mode = (uint8_t) intra_prediction(encoder->in.cur_pic->y_data,
        encoder->in.width, recShift, (LCU_WIDTH >> (depth)) * 2 + 8,
        x_ctb * (LCU_WIDTH >> (MAX_DEPTH)), y_ctb * (LCU_WIDTH >> (MAX_DEPTH)),
        width, pred, width, &cur_cu->intra.cost);
    //free(pred);
    //free(rec);
  }

  // Split and search to max_depth
  if (depth < MAX_INTRA_SEARCH_DEPTH && depth < MAX_INTER_SEARCH_DEPTH) {
    // Split blocks and remember to change x and y block positions
    uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
    search_tree(encoder, x_ctb,          y_ctb,          depth + 1);
    search_tree(encoder, x_ctb + change, y_ctb,          depth + 1);
    search_tree(encoder, x_ctb,          y_ctb + change, depth + 1);
    search_tree(encoder, x_ctb + change, y_ctb + change, depth + 1);
  }
}

/**
 * \brief
 */
uint32_t search_best_mode(encoder_control *encoder, 
                          uint16_t x_ctb, uint16_t y_ctb, uint8_t depth)
{
  cu_info *cur_cu = &encoder->in.cur_pic->cu_array[depth]
                     [x_ctb + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)];
  uint32_t best_intra_cost = cur_cu->intra.cost;
  uint32_t best_inter_cost = cur_cu->inter.cost;
  uint32_t lambda_cost = (4 * g_lambda_cost[encoder->QP]) << 4; //<<5; //TODO: Correct cost calculation

  if (depth < MAX_INTRA_SEARCH_DEPTH && depth < MAX_INTER_SEARCH_DEPTH) {
    uint32_t cost = lambda_cost;
    uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
    cost += search_best_mode(encoder, x_ctb,          y_ctb,          depth + 1);
    cost += search_best_mode(encoder, x_ctb + change, y_ctb,          depth + 1);
    cost += search_best_mode(encoder, x_ctb,          y_ctb + change, depth + 1);
    cost += search_best_mode(encoder, x_ctb + change, y_ctb + change, depth + 1);

    if (cost < best_intra_cost && cost < best_inter_cost)
    {
      // Better value was found at a lower level.
      return cost;
    }
  } 

  // If search hasn't been peformed at all for this block, the cost will be
  // max value, so it is safe to just compare costs. It just has to be made
  // sure that no value overflows.
  if (best_inter_cost <= best_intra_cost) {
    inter_set_block(encoder->in.cur_pic, x_ctb, y_ctb, depth, cur_cu);
    return best_inter_cost;
  } else {
    intra_set_block_mode(encoder->in.cur_pic, x_ctb, y_ctb, depth,
        cur_cu->intra.mode);
    return best_intra_cost;
  }
}



/**
 * \brief
 */
void search_slice_data(encoder_control *encoder)
{
  int16_t x_lcu, y_lcu;
  FILE *fp = 0, *fp2 = 0;

  if (RENDER_CU) {
    fp = open_cu_file("cu_search.html");
    fp2 = open_cu_file("cu_best.html");
  }

  // Loop through every LCU in the slice
  for (y_lcu = 0; y_lcu < encoder->in.height_in_lcu; y_lcu++) {
    for (x_lcu = 0; x_lcu < encoder->in.width_in_lcu; x_lcu++) {
      uint8_t depth = 0;
      // Recursive function for looping through all the sub-blocks
      search_tree(encoder, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, depth);
      if (RENDER_CU) {
        render_cu_file(encoder, encoder->in.cur_pic, depth, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, fp);
      }

      // Decide actual coding modes
      search_best_mode(encoder, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, depth);
      if (RENDER_CU) {
        render_cu_file(encoder, encoder->in.cur_pic, depth, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, fp2);
      }

      encode_block_residual(encoder, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, depth);

    }
  }



  if (RENDER_CU && fp) {
    close_cu_file(fp);
    fp = 0;
  }
  if (RENDER_CU && fp2) {
    close_cu_file(fp2);
    fp2 = 0;
  }
}
