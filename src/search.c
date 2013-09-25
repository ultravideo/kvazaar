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
#define USE_INTRA_IN_P 0
#define RENDER_CU 1


/**
 * \brief Search motions vectors for a block and all it's sub-blocks.
 *
 * \param pic
 * \param pic_data  picture color data starting from the block MV is being searched for.
 * \param ref_data  picture color data starting from the beginning of reference pic.
 * \param cur_cu
 */
void search_motion_vector(picture *pic, uint8_t *pic_data, uint8_t *ref_data,
                          cu_info *cur_cu, unsigned step, 
                          int orig_x, int orig_y, int x, int y, unsigned depth)
{
  // TODO: Inter: Handle non-square blocks.
  int block_width = CU_WIDTH_FROM_DEPTH(depth);
  int block_height = block_width;
  unsigned cost;

  // TODO: Inter: Calculating error outside picture borders.
  // This prevents choosing vectors that need interpolating of borders to work.
  if (orig_x + x < 0 || orig_y + y < 0 || orig_x + x > pic->width - block_width
      || orig_y + y > pic->height - block_height) return;

  cost = sad(pic_data, &ref_data[(orig_y + y) * pic->width + (orig_x + x)],
      block_width, block_height, pic->width) + 1;
  if (cost < cur_cu->inter.cost) {
    cur_cu->inter.cost = cost;
    cur_cu->inter.mv[0] = x << 2;
    cur_cu->inter.mv[1] = y << 2;
  }

  step /= 2;
  if (step > 0) {
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y,
                         x, y - step, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y,
                         x - step, y, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y,
                         x + step, y, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y,
                         x, y + step, depth);
  }
}

/**
 * \brief
 */
void search_buildReferenceBorder(picture *pic, int32_t x_ctb, int32_t y_ctb,
                                 int16_t outwidth, int16_t *dst, 
                                 int32_t dststride, int8_t chroma)
{
  int32_t left_col; // left column iterator
  int16_t val;      // variable to store extrapolated value
  int32_t i;        // index iterator
  int16_t dc_val = 1 << (g_bitdepth - 1); // default predictor value
  int32_t top_row;  // top row iterator
  int32_t src_width = (pic->width >> (chroma ? 1 : 0));   // source picture width
  int32_t src_height = (pic->height >> (chroma ? 1 : 0)); // source picture height
  uint8_t *src_pic = (!chroma) ? pic->y_data : ((chroma == 1) ? pic->u_data : pic->v_data); // input picture pointer
  int16_t scu_width = LCU_WIDTH >> (MAX_DEPTH + (chroma ? 1 : 0)); // Smallest Coding Unit width
  uint8_t *src_shifted = &src_pic[x_ctb * scu_width + (y_ctb * scu_width) * src_width]; // input picture pointer shifted to start from the left-top corner of the current block
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
  cu_info *cur_cu = &encoder->in.cur_pic->cu_array[depth][x_ctb + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)];

  cur_cu->intra.cost = 0xffffffff;
  cur_cu->inter.cost = 0xffffffff;

  // Force split on border
  if (depth != MAX_DEPTH) {
    if (border) {
      // Split blocks and remember to change x and y block positions
      uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
      SET_SPLITDATA(cur_cu, 1);
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
      // We don't need to do anything else here
      return;
    }
  }

  // INTER SEARCH
  if (depth >= MIN_INTER_SEARCH_DEPTH && depth <= MAX_INTER_SEARCH_DEPTH
      && encoder->in.cur_pic->slicetype != SLICE_I) {
    // Motion estimation on P-frame
    if (encoder->in.cur_pic->slicetype != SLICE_B) {

    }

    {
      unsigned mv[2] = { 0, 0 }; // TODO: Take initial MV from adjacent blocks.
      picture *cur_pic = encoder->in.cur_pic;

      picture *ref_pic = encoder->ref->pics[0];

      int x = x_ctb * CU_MIN_SIZE_PIXELS;
      int y = y_ctb * CU_MIN_SIZE_PIXELS;
      uint8_t *cur_data = &cur_pic->y_data[(y * cur_pic->width) + x];
      search_motion_vector(cur_pic, cur_data, ref_pic->y_data, cur_cu, 8, x, y,
                           0, 0, depth);
    }

    cur_cu->type = CU_INTER;
    cur_cu->inter.mv_dir = 1;
    inter_set_block(encoder->in.cur_pic, x_ctb, y_ctb, depth, cur_cu);
  }

  // INTRA SEARCH
  if (depth >= MIN_INTRA_SEARCH_DEPTH && depth <= MAX_INTRA_SEARCH_DEPTH
      && (encoder->in.cur_pic->slicetype == SLICE_I || USE_INTRA_IN_P)) {
    int x = 0, y = 0;
    uint8_t *base = &encoder->in.cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH))) * encoder->in.width];
    uint32_t width = LCU_WIDTH >> depth;

    // INTRAPREDICTION
    int16_t pred[LCU_WIDTH * LCU_WIDTH + 1];
    int16_t rec[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8)];
    int16_t *recShift = &rec[(LCU_WIDTH >> (depth)) * 2 + 8 + 1];

    //int16_t *pred = (int16_t*)malloc(LCU_WIDTH*LCU_WIDTH*sizeof(int16_t));
    //int16_t *rec = (int16_t*)malloc((LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)*sizeof(int16_t));

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
  cu_info *cur_cu = &encoder->in.cur_pic->cu_array[depth][x_ctb
      + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)];
  uint32_t best_intra_cost = cur_cu->intra.cost;
  uint32_t best_inter_cost = cur_cu->inter.cost;
  uint32_t best_cost = 0;
  uint32_t cost = 0;
  uint32_t lambdaCost = (4 * g_lambda_cost[encoder->QP]) << 4; //<<5; //TODO: Correct cost calculation

  // Split and search to max_depth
  if (depth != MAX_INTRA_SEARCH_DEPTH) {
    // Split blocks and remember to change x and y block positions
    uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
    cost =  search_best_mode(encoder, x_ctb,          y_ctb,          depth + 1);
    cost += search_best_mode(encoder, x_ctb + change, y_ctb,          depth + 1);
    cost += search_best_mode(encoder, x_ctb,          y_ctb + change, depth + 1);
    cost += search_best_mode(encoder, x_ctb + change, y_ctb + change, depth + 1);

    // We split if the cost is better (0 cost -> not checked)
    if (cost != 0 
        && (best_intra_cost != 0 && cost + lambdaCost < best_intra_cost)
        && (best_inter_cost != 0
            && cost + lambdaCost < best_inter_cost
            && encoder->in.cur_pic->slicetype != SLICE_I))
    {
      // Set split to 1
      picture_set_block_split(encoder->in.cur_pic, x_ctb, y_ctb, depth, 1);
      best_cost = cost + lambdaCost;
    } else if (best_inter_cost != 0 // Else, check if inter cost is smaller or the same as intra 
        && (best_inter_cost <= best_intra_cost || best_intra_cost == 0)
        && encoder->in.cur_pic->slicetype != SLICE_I)
    {
      // Set split to 0 and mode to inter.mode
      picture_set_block_split(encoder->in.cur_pic, x_ctb, y_ctb, depth, 0);
      inter_set_block(encoder->in.cur_pic, x_ctb, y_ctb, depth, cur_cu);
      best_cost = best_inter_cost;
    } else { // Else, dont split and recursively set block mode
      // Set split to 0 and mode to intra.mode
      picture_set_block_split(encoder->in.cur_pic, x_ctb, y_ctb, depth, 0);
      intra_set_block_mode(encoder->in.cur_pic, x_ctb, y_ctb, depth,
          cur_cu->intra.mode);
      best_cost = best_intra_cost;
    }
  } else if (best_inter_cost != 0
             && (best_inter_cost <= best_intra_cost || best_intra_cost == 0)
             && encoder->in.cur_pic->slicetype != SLICE_I)
  {
    // Set split to 0 and mode to inter.mode
    picture_set_block_split(encoder->in.cur_pic, x_ctb, y_ctb, depth, 0);
    inter_set_block(encoder->in.cur_pic, x_ctb, y_ctb, depth, cur_cu);
    best_cost = best_inter_cost;
  } else {
    // Set split to 0 and mode to intra.mode
    picture_set_block_split(encoder->in.cur_pic, x_ctb, y_ctb, depth, 0);
    intra_set_block_mode(encoder->in.cur_pic, x_ctb, y_ctb, depth,
        cur_cu->intra.mode);
    best_cost = best_intra_cost;
  }

  return best_cost;
}

/**
 * \brief
 */
void search_slice_data(encoder_control *encoder)
{
  int16_t x_lcu, y_lcu;
  FILE *fp = 0, *fp2 = 0;

  if (RENDER_CU && encoder->frame == 1) {
    fp = open_cu_file("cu_search.html");
    fp2 = open_cu_file("cu_best.html");
  }

  // Loop through every LCU in the slice
  for (y_lcu = 0; y_lcu < encoder->in.height_in_lcu; y_lcu++) {
    for (x_lcu = 0; x_lcu < encoder->in.width_in_lcu; x_lcu++) {
      uint8_t depth = 0;
      // Recursive function for looping through all the sub-blocks
      search_tree(encoder, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, depth);
      if (RENDER_CU && encoder->frame == 1) {
        render_cu_file(encoder, depth, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, fp);
      }

      // Decide actual coding modes
      search_best_mode(encoder, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, depth);
      if (RENDER_CU && encoder->frame == 1) {
        render_cu_file(encoder, depth, x_lcu << MAX_DEPTH, y_lcu << MAX_DEPTH, fp2);
      }
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
