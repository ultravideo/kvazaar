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

/**
 * \file
 * \brief Functions for handling intra frames.
 */

#include "intra.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "encoder.h"
#include "transform.h"
#include "rdo.h"


const uint8_t kvz_intra_hor_ver_dist_thres[5] = {0,7,1,0,0};


/**
 * \brief Set intrablock mode (and init typedata)
 * \param pic picture to use
 * \param xCtb x CU position (smallest CU)
 * \param yCtb y CU position (smallest CU)
 * \param depth current CU depth
 * \param mode mode to set
 * \returns Void
 */
void kvz_intra_set_block_mode(videoframe_t *frame,uint32_t x_cu, uint32_t y_cu, uint8_t depth, uint8_t mode, uint8_t part_mode)
{
  uint32_t x, y;
  int block_scu_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);

  if (part_mode == SIZE_NxN) {
    cu_info_t *cur_cu = kvz_videoframe_get_cu(frame, x_cu, y_cu);
    // Modes are already set.
    cur_cu->depth = depth;
    cur_cu->type = CU_INTRA;
    cur_cu->tr_depth = depth + 1;
    return;
  }

  // Loop through all the blocks in the area of cur_cu
  for (y = y_cu; y < y_cu + block_scu_width; y++) {
    for (x = x_cu; x < x_cu + block_scu_width; x++) {
      cu_info_t *cur_cu = kvz_videoframe_get_cu(frame, x_cu, y_cu);
      cur_cu->depth = depth;
      cur_cu->type = CU_INTRA;
      cur_cu->intra[0].mode = mode;
      cur_cu->intra[1].mode = mode;
      cur_cu->intra[2].mode = mode;
      cur_cu->intra[3].mode = mode;
      cur_cu->part_size = part_mode;
      cur_cu->tr_depth = depth;
    }
  }
}

/**
 * \brief get intrablock mode
 * \param pic picture data to use
 * \param picwidth width of the picture data
 * \param xpos x-position
 * \param ypos y-position
 * \param width block width
 * \returns DC prediction
*/
kvz_pixel kvz_intra_get_dc_pred(const kvz_pixel *pic, uint16_t picwidth, uint8_t width)
{
  int32_t i, sum = 0;

  // pixels on top and left
  for (i = -picwidth; i < width - picwidth; i++) {
    sum += pic[i];
  }
  for (i = -1; i < width * picwidth - 1; i += picwidth) {
    sum += pic[i];
  }

  // return the average
  return (kvz_pixel)((sum + width) / (width + width));
}

/**
 * \brief Function for deriving intra luma predictions
 * \param pic picture to use
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param preds output buffer for 3 predictions
 * \returns (predictions are found)?1:0
 */
int8_t kvz_intra_get_dir_luma_predictor(const uint32_t x, const uint32_t y, int8_t* preds,
                                    const cu_info_t * const cur_cu, const cu_info_t * const left_cu, const cu_info_t * const above_cu)
{
  int y_cu = y>>3;

  // The default mode if block is not coded yet is INTRA_DC.
  int8_t left_intra_dir  = 1;
  int8_t above_intra_dir = 1;

  if (x & 4) {
    // If current CU is NxN and PU is on the right half, take mode from the
    // left half of the same CU.
    left_intra_dir = cur_cu->intra[PU_INDEX(0, y >> 2)].mode;
  } else if (left_cu && left_cu->type == CU_INTRA) {
    // Otherwise take the mode from the right side of the CU on the left.
    left_intra_dir = left_cu->intra[PU_INDEX(1, y >> 2)].mode;
  }

  if (y & 4) {
    // If current CU is NxN and PU is on the bottom half, take mode from the
    // top half of the same CU.
    above_intra_dir = cur_cu->intra[PU_INDEX(x >> 2, 0)].mode;
  } else if (above_cu && above_cu->type == CU_INTRA &&
             (y_cu * (LCU_WIDTH>>MAX_DEPTH)) % LCU_WIDTH != 0)
  {
    // Otherwise take the mode from the bottom half of the CU above.
    above_intra_dir = above_cu->intra[PU_INDEX(x >> 2, 1)].mode;
  }

  // If the predictions are the same, add new predictions
  if (left_intra_dir == above_intra_dir) {
    if (left_intra_dir > 1) { // angular modes
      preds[0] = left_intra_dir;
      preds[1] = ((left_intra_dir + 29) % 32) + 2;
      preds[2] = ((left_intra_dir - 1 ) % 32) + 2;
    } else { //non-angular
      preds[0] = 0;//PLANAR_IDX;
      preds[1] = 1;//DC_IDX;
      preds[2] = 26;//VER_IDX;
    }
  } else { // If we have two distinct predictions
    preds[0] = left_intra_dir;
    preds[1] = above_intra_dir;

    // add planar mode if it's not yet present
    if (left_intra_dir && above_intra_dir ) {
      preds[2] = 0; // PLANAR_IDX;
    } else {  // Add DC mode if it's not present, otherwise 26.
      preds[2] =  (left_intra_dir+above_intra_dir)<2? 26 : 1;
    }
  }

  return 1;
}

/**
 * \brief Intra filtering of the border samples
 * \param ref reference picture data
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param depth current CU depth
 * \param preds output buffer for 3 predictions
 * \returns (predictions are found)?1:0
 */
void kvz_intra_filter(kvz_pixel *ref, int32_t stride,int32_t width, int8_t mode)
{
  #define FWIDTH (LCU_WIDTH*2+1)
  kvz_pixel filtered[FWIDTH * FWIDTH]; //!< temporary buffer for filtered samples
  kvz_pixel *filteredShift = &filtered[FWIDTH+1]; //!< pointer to temporary buffer with offset (1,1)
  int x,y;

  if (!mode) {
    // pF[ -1 ][ -1 ] = ( p[ -1 ][ 0 ] + 2*p[ -1 ][ -1 ] + p[ 0 ][ -1 ] + 2 )  >>  2	(8 35)
    filteredShift[-FWIDTH-1] = (ref[-1] + 2*ref[-(int32_t)stride-1] + ref[-(int32_t)stride] + 2) >> 2;

    // pF[ -1 ][ y ] = ( p[ -1 ][ y + 1 ] + 2*p[ -1 ][ y ] + p[ -1 ][ y - 1 ] + 2 )  >>  2 for y = 0..nTbS * 2 - 2	(8 36)
    for (y = 0; y < (int32_t)width * 2 - 1; y++) {
      filteredShift[y*FWIDTH-1] = (ref[(y + 1) * stride - 1] + 2*ref[y * stride - 1] + ref[(y - 1) * stride - 1] + 2) >> 2;
    }

    // pF[ -1 ][ nTbS * 2 - 1 ] = p[ -1 ][ nTbS * 2 - 1 ]		(8 37)
    filteredShift[(width * 2 - 1) * FWIDTH - 1] = ref[(width * 2 - 1) * stride - 1];

    // pF[ x ][ -1 ] = ( p[ x - 1 ][ -1 ] + 2*p[ x ][ -1 ] + p[ x + 1 ][ -1 ] + 2 )  >>  2 for x = 0..nTbS * 2 - 2	(8 38)
    for(x = 0; x < (int32_t)width*2-1; x++) {
      filteredShift[x - FWIDTH] = (ref[x - 1 - stride] + 2*ref[x - stride] + ref[x + 1 - stride] + 2) >> 2;
    }

    // pF[ nTbS * 2 - 1 ][ -1 ] = p[ nTbS * 2 - 1 ][ -1 ]
    filteredShift[(width * 2 - 1) - FWIDTH] = ref[(width * 2 - 1) - stride];

    // Copy filtered samples to the input array
    for (x = -1; x < (int32_t)width * 2; x++) {
      ref[x - stride] = filtered[x + 1];
    }
    for(y = 0; y < (int32_t)width * 2; y++)  {
      ref[y * stride - 1] = filtered[(y + 1) * FWIDTH];
    }
  } else  {
    printf("UNHANDLED: %s: %d\r\n", __FILE__, __LINE__);
    exit(1);
  }
  #undef FWIDTH
}


static void intra_filter_reference(int_fast8_t log2_width, kvz_intra_references *refs)
{
  if (refs->filtered_initialized) {
    return;
  } else {
    refs->filtered_initialized = true;
  }

  const int_fast8_t ref_width = 2 * (1 << log2_width) + 1;
  kvz_intra_ref *ref = &refs->ref;
  kvz_intra_ref *filtered_ref = &refs->filtered_ref;

  filtered_ref->left[0] = (ref->left[1] + 2 * ref->left[0] + ref->top[1] + 2) / 4;
  filtered_ref->top[0] = filtered_ref->left[0];

  for (int_fast8_t y = 1; y < ref_width - 1; ++y) {
    kvz_pixel *p = &ref->left[y];
    filtered_ref->left[y] = (p[-1] + 2 * p[0] + p[1] + 2) / 4;
  }
  filtered_ref->left[ref_width - 1] = ref->left[ref_width - 1];

  for (int_fast8_t x = 1; x < ref_width - 1; ++x) {
    kvz_pixel *p = &ref->top[x];
    filtered_ref->top[x] = (p[-1] + 2 * p[0] + p[1] + 2) / 4;
  }
  filtered_ref->top[ref_width - 1] = ref->top[ref_width - 1];
}


static void post_process_intra_angular(
  unsigned width,
  unsigned stride,
  const kvz_pixel *ref,
  kvz_pixel *block)
{
  kvz_pixel ref2 = ref[0];
  for (unsigned i = 0; i < width; i++) {
    kvz_pixel val = block[i * stride];
    kvz_pixel ref1 = ref[i + 1];
    block[i * stride] = CLIP_TO_PIXEL(val + ((ref1 - ref2) >> 1));
  }
}


/**
 * \brief Generage angular predictions.
 * \param log2_width    Log2 of width, range 2..5.
 * \param intra_mode    Angular mode in range 2..34.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_intra_pred_angular(
  const int_fast8_t log2_width,
  const int_fast8_t intra_mode,
  const kvz_pixel *const in_ref_above,
  const kvz_pixel *const in_ref_left,
  kvz_pixel *const dst)
{
  assert(log2_width >= 2 && log2_width <= 5);
  assert(intra_mode >= 2 && intra_mode <= 34);

  static const int8_t modedisp2sampledisp[9] = {0, 2, 5, 9, 13, 17, 21, 26, 32};
  static const int16_t modedisp2invsampledisp[9] = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / sampledisp

  // Temporary buffer for modes 11-25.
  // It only needs to be big enough to hold indices from -width to width-1.
  kvz_pixel tmp_ref[2 * 32];
  const int_fast8_t width = 1 << log2_width;

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 18;
  // Modes distance to horizontal or vertical mode.
  const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
  // Sample displacement per column in fractions of 32.
  const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

  // Pointer for the reference we are interpolating from.
  const kvz_pixel *ref_main;
  // Pointer for the other reference.
  const kvz_pixel *ref_side;

  // Set ref_main and ref_side such that, when indexed with 0, they point to
  // index 0 in block coordinates.
  if (sample_disp < 0) {
    // Negative sample_disp means, we need to use both references.

    ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

    // Move the reference pixels to start from the middle to the later half of
    // the tmp_ref, so there is room for negative indices.
    for (int_fast8_t x = -1; x < width; ++x) {
      tmp_ref[x + width] = ref_main[x];
    }
    // Get a pointer to block index 0 in tmp_ref.
    ref_main = &tmp_ref[width];

    // Extend the side reference to the negative indices of main reference.
    int_fast32_t col_sample_disp = 128; // rounding for the ">> 8"
    int_fast16_t inv_abs_sample_disp = modedisp2invsampledisp[abs(mode_disp)];
    int_fast8_t most_negative_index = (width * sample_disp) >> 5;
    for (int_fast8_t x = -2; x >= most_negative_index; --x) {
      col_sample_disp += inv_abs_sample_disp;
      int_fast8_t side_index = col_sample_disp >> 8;
      tmp_ref[x + width] = ref_side[side_index - 1];
    }
  } else {
    // sample_disp >= 0 means we don't need to refer to negative indices,
    // which means we can just use the references as is.
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
    ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
  }

  if (sample_disp != 0) {
    // The mode is not horizontal or vertical, we have to do interpolation.

    int_fast16_t delta_pos = 0;
    for (int_fast8_t y = 0; y < width; ++y) {
      delta_pos += sample_disp;
      int_fast8_t delta_int = delta_pos >> 5;
      int_fast8_t delta_fract = delta_pos & (32 - 1);

      if (delta_fract) {
        // Do linear filtering
        for (int_fast8_t x = 0; x < width; ++x) {
          kvz_pixel ref1 = ref_main[x + delta_int];
          kvz_pixel ref2 = ref_main[x + delta_int + 1];
          dst[y * width + x] = ((32 - delta_fract) * ref1 + delta_fract * ref2 + 16) >> 5;
        }
      } else {
        // Just copy the integer samples
        for (int_fast8_t x = 0; x < width; x++) {
          dst[y * width + x] = ref_main[x + delta_int];
        }
      }
    }
  } else {
    // Mode is horizontal or vertical, just copy the pixels.

    for (int_fast8_t y = 0; y < width; ++y) {
      for (int_fast8_t x = 0; x < width; ++x) {
        dst[y * width + x] = ref_main[x];
      }
    }
  }

  // Flip the block if this is was a horizontal mode.
  if (!vertical_mode) {
    for (int_fast8_t y = 0; y < width - 1; ++y) {
      for (int_fast8_t x = y + 1; x < width; ++x) {
        SWAP(dst[y * width + x], dst[x * width + y], kvz_pixel);
      }
    }
  }
}


/**
 * \brief Generage planar prediction.
 * \param log2_width    Log2 of width, range 2..5.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_intra_pred_planar(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const dst)
{
  assert(log2_width >= 2 && log2_width <= 5);

  const int_fast8_t width = 1 << log2_width;
  const kvz_pixel top_right = ref_top[width + 1];
  const kvz_pixel bottom_left = ref_left[width + 1];

#if 0
  // Unoptimized version for reference.
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int_fast16_t hor = (width - 1 - x) * ref_left[y + 1] + (x + 1) * top_right;
      int_fast16_t ver = (width - 1 - y) * ref_top[x + 1] + (y + 1) * bottom_left;
      dst[y * width + x] = (ver + hor + width) >> (log2_width + 1);
    }
  }
#else
  int_fast16_t top[32];
  for (int i = 0; i < width; ++i) {
    top[i] = ref_top[i + 1] << log2_width;
  }

  for (int y = 0; y < width; ++y) {
    int_fast16_t hor = (ref_left[y + 1] << log2_width) + width;
    for (int x = 0; x < width; ++x) {
      hor += top_right - ref_left[y + 1];
      top[x] += bottom_left - ref_top[x + 1];
      dst[y * width + x] = (hor + top[x]) >> (log2_width + 1);
    }
  }
#endif
}


/**
* \brief Generage planar prediction.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
*/
static void kvz_intra_pred_dc(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const out_block)
{
  int_fast8_t width = 1 << log2_width;

  int_fast16_t sum = 0;
  for (int_fast8_t i = 0; i < width; ++i) {
    sum += ref_top[i + 1];
    sum += ref_left[i + 1];
  }

  const kvz_pixel dc_val = (sum + width) >> (log2_width + 1);
  const int_fast16_t block_size = 1 << (log2_width * 2);

  for (int_fast16_t i = 0; i < block_size; ++i) {
    out_block[i] = dc_val;
  }
}


/**
* \brief Generage intra DC prediction with post filtering applied.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
*/
static void kvz_intra_pred_filtered_dc(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const out_block)
{
  assert(log2_width >= 2 && log2_width <= 5);

  const int_fast8_t width = 1 << log2_width;

  int_fast16_t sum = 0;
  for (int_fast8_t i = 0; i < width; ++i) {
    sum += ref_top[i + 1];
    sum += ref_left[i + 1];
  }

  const kvz_pixel dc_val = (sum + width) >> (log2_width + 1);

  // Filter top-left with ([1 2 1] / 4)
  out_block[0] = (ref_left[1] + 2 * dc_val + ref_top[1] + 2) / 4;

  // Filter rest of the boundary with ([1 3] / 4)
  for (int_fast8_t x = 1; x < width; ++x) {
    out_block[x] = (ref_top[x + 1] + 3 * dc_val + 2) / 4;
  }
  for (int_fast8_t y = 1; y < width; ++y) {
    out_block[y * width] = (ref_left[y + 1] + 3 * dc_val + 2) / 4;
    for (int_fast8_t x = 1; x < width; ++x) {
      out_block[y * width + x] = dc_val;
    }
  }
}


void kvz_intra_get_pred(const encoder_control_t * const encoder, const kvz_pixel *rec, const kvz_pixel *rec_filtered, int recstride, kvz_pixel *dst, int width, int mode, int is_chroma)
{
  const kvz_pixel *ref_pixels = rec;
  if (is_chroma || mode == 1 || width == 4) {
    // For chroma, DC and 4x4 blocks, always use unfiltered reference.
  } else if (mode == 0) {
    // Otherwise, use filtered for planar.
    ref_pixels = rec_filtered;
  } else {
    // Angular modes use smoothed reference pixels, unless the mode is close
    // to being either vertical or horizontal.
    int filter_threshold = kvz_intra_hor_ver_dist_thres[g_to_bits[width]];
    int dist_from_vert_or_hor = MIN(abs(mode - 26), abs(mode - 10));
    if (dist_from_vert_or_hor > filter_threshold) {
      ref_pixels = rec_filtered;
    }
  }

  if (mode == 0) {
    kvz_intra_get_planar_pred(ref_pixels, recstride, width, dst, width);
  } else if (mode == 1) {
    int i;
    kvz_pixel val = kvz_intra_get_dc_pred(ref_pixels, recstride, width);
    for (i = 0; i < width * width; i++) {
      dst[i] = val;
    }
    // Do extra post filtering for edge pixels of luma DC mode.
    if (!is_chroma && width < 32) {
      kvz_intra_dc_pred_filtering(ref_pixels, recstride, dst, width, width, width);
    }
  } else {
    int filter = !is_chroma && width < 32;
    kvz_intra_get_angular_pred(encoder, ref_pixels, recstride, dst, width, width, mode, filter);
  }
}


void kvz_intra_get_pred_new(
  kvz_intra_references *refs,
  int_fast8_t log2_width,
  int_fast8_t mode,
  color_t color,
  kvz_pixel *dst)
{
  const int_fast8_t width = 1 << log2_width;

  const kvz_intra_ref *used_ref = &refs->ref;
  if (color != COLOR_Y || mode == 1 || width == 4) {
    // For chroma, DC and 4x4 blocks, always use unfiltered reference.
  } else if (mode == 0) {
    // Otherwise, use filtered for planar.
    used_ref = &refs->filtered_ref;
  } else {
    // Angular modes use smoothed reference pixels, unless the mode is close
    // to being either vertical or horizontal.
    int filter_threshold = kvz_intra_hor_ver_dist_thres[g_to_bits[width]];
    int dist_from_vert_or_hor = MIN(abs(mode - 26), abs(mode - 10));
    if (dist_from_vert_or_hor > filter_threshold) {
      used_ref = &refs->filtered_ref;
    }
  }

  if (used_ref == &refs->filtered_ref && !refs->filtered_initialized) {
    intra_filter_reference(log2_width, refs);
  }

  if (mode == 0) {
    kvz_intra_pred_planar(log2_width, used_ref->top, used_ref->left, dst);
  } else if (mode == 1) {
    // Do extra post filtering for edge pixels of luma DC mode.
    if (color == COLOR_Y && width < 32) {
      kvz_intra_pred_filtered_dc(log2_width, used_ref->top, used_ref->left, dst);
    } else {
      kvz_intra_pred_dc(log2_width, used_ref->top, used_ref->left, dst);
    }
  } else {
    kvz_intra_pred_angular(log2_width, mode, used_ref->top, used_ref->left, dst);
    if (color == COLOR_Y && width < 32) {
      if (mode == 10) {
        post_process_intra_angular(width, 1, used_ref->top, dst);
      } else if (mode == 26) {
        post_process_intra_angular(width, width, used_ref->left, dst);
      }
    }
  }
}


/**
 * \brief Reconstruct intra block according to prediction
 * \param rec reconstructed picture data
 * \param recstride reconstructed picture stride
 * \param width block size to predict
 * \param dst destination buffer for best prediction
 * \param dststride destination width
 * \param mode intra mode to use
 * \param chroma chroma-block flag
*/
void kvz_intra_recon(const encoder_control_t * const encoder, kvz_pixel* rec, int32_t recstride, uint32_t width, kvz_pixel* dst, int32_t dststride, int8_t mode, int8_t chroma)
{
  kvz_pixel pred[LCU_WIDTH * LCU_WIDTH];
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
    kvz_intra_filter(recf, recstride, width, 0);
  }

  kvz_intra_get_pred(encoder, rec, recf, recstride, pred, width, mode, chroma);

  kvz_pixels_blit(pred, dst, width, width, width, dststride);
}

void kvz_intra_recon_new(
  kvz_intra_references *refs, 
  uint32_t log2_width, 
  kvz_pixel* dst, 
  int32_t dst_stride, 
  int8_t mode, 
  color_t color)
{
  kvz_pixel pred[32 * 32];
  const int_fast8_t width = 1 << log2_width;
  
  kvz_intra_get_pred_new(refs, log2_width, mode, color, pred);

  kvz_pixels_blit(pred, dst, width, width, width, dst_stride);
}

/**
 * \brief Build top and left borders for a reference block.
 * \param pic picture to use as a source
 * \param outwidth width of the prediction block
 * \param chroma signaling if chroma is used, 0 = luma, 1 = U and 2 = V
 *
 * The end result is 2*width+8 x 2*width+8 array, with only the top and left
 * edge pixels filled with the reconstructed pixels.
 */
void kvz_intra_build_reference_border(const encoder_control_t * const encoder, int32_t x_luma, int32_t y_luma, int16_t out_width,
                                      kvz_pixel *dst, int32_t dst_stride, int8_t chroma,
                                      int32_t pic_width, int32_t pic_height,
                                      lcu_t *lcu)
{
  // Some other function might make use of the arrays num_ref_pixels_top and
  // num_ref_pixels_left in the future, but until that happens lets leave
  // them here.

  /**
   * \brief Table for looking up the number of intra reference pixels based on
   *        prediction units coordinate within an LCU.
   *
   * This table was generated by "tools/generate_ref_pixel_tables.py".
   */
  static const uint8_t num_ref_pixels_top[16][16] = {
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 }
  };

  /**
   * \brief Table for looking up the number of intra reference pixels based on
   *        prediction units coordinate within an LCU.
   *
   * This table was generated by "tools/generate_ref_pixel_tables.py".
   */
  static const uint8_t num_ref_pixels_left[16][16] = {
    { 64,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 }
  };

  const kvz_pixel dc_val = 1 << (encoder->bitdepth - 1);
  const int is_chroma = chroma ? 1 : 0;

  // input picture pointer
  //const pixel * const src = (!chroma) ? pic->y_recdata : ((chroma == 1) ? pic->u_recdata : pic->v_recdata);

  // Convert luma coordinates to chroma coordinates for chroma.
  const int x = chroma ? x_luma / 2 : x_luma;
  const int y = chroma ? y_luma / 2 : y_luma;

  const int y_in_lcu = y_luma % LCU_WIDTH;
  const int x_in_lcu = x_luma % LCU_WIDTH;

  int x_local = (x_luma&0x3f)>>is_chroma, y_local = (y_luma&0x3f)>>is_chroma;

  kvz_pixel *left_ref = !chroma ? &lcu->left_ref.y[1] : (chroma == 1) ? &lcu->left_ref.u[1] : &lcu->left_ref.v[1];
  kvz_pixel *top_ref  = !chroma ? &lcu->top_ref.y[1]  : (chroma == 1) ? &lcu->top_ref.u[1]  : &lcu->top_ref.v[1];
  kvz_pixel *rec_ref  = !chroma ? lcu->rec.y : (chroma == 1) ? lcu->rec.u : lcu->rec.v;

  kvz_pixel *left_border = &left_ref[y_local];
  kvz_pixel *top_border = &top_ref[x_local];
  uint32_t left_stride = 1;

  if(x_local) {
    left_border = &rec_ref[x_local - 1 + y_local * (LCU_WIDTH>>is_chroma)];
    left_stride = LCU_WIDTH>>is_chroma;
  }

  if(y_local) {
    top_border = &rec_ref[x_local + (y_local - 1) * (LCU_WIDTH>>is_chroma)];
  }

  // Copy pixels for left edge.
  if (x > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int num_ref_pixels = num_ref_pixels_left[y_in_lcu / 4][x_in_lcu / 4] >> is_chroma;
    int i;
    kvz_pixel nearest_pixel;

    // Max pixel we can copy from src is yy + outwidth - 1 because the dst
    // extends one pixel to the left.
    num_ref_pixels = MIN(num_ref_pixels, out_width - 1);
    // There are no coded pixels below the frame.
    num_ref_pixels = MIN(num_ref_pixels, pic_height - y);
    // There are no coded pixels below the bottom of the LCU due to raster
    // scan order.
    num_ref_pixels = MIN(num_ref_pixels, (LCU_WIDTH - y_in_lcu) >> is_chroma);

    // Copy pixels from coded CUs.
    for (i = 0; i < num_ref_pixels; ++i) {
      dst[(i + 1) * dst_stride] = left_border[i*left_stride];
    }
    // Extend the last pixel for the rest of the reference values.
    nearest_pixel = dst[i * dst_stride];
    for (i = num_ref_pixels; i < out_width - 1; ++i) {
      dst[i * dst_stride] = nearest_pixel;
    }
  } else {
    // If we are on the left edge, extend the first pixel of the top row.
    kvz_pixel nearest_pixel = y > 0 ? top_border[0] : dc_val;
    int i;
    for (i = 1; i < out_width - 1; i++) {
      dst[i * dst_stride] = nearest_pixel;
    }
  }

  // Copy pixels for top edge.
  if (y > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int num_ref_pixels = num_ref_pixels_top[y_in_lcu / 4][x_in_lcu / 4] >> is_chroma;
    int i;
    kvz_pixel nearest_pixel;

    // Max pixel we can copy from src is yy + outwidth - 1 because the dst
    // extends one pixel to the left.
    num_ref_pixels = MIN(num_ref_pixels, out_width - 1);
    // All LCUs in the row above have been coded.
    num_ref_pixels = MIN(num_ref_pixels, pic_width - x);

    // Copy pixels from coded CUs.
    for (i = 0; i < num_ref_pixels; ++i) {
      dst[i + 1] = top_border[i];
    }
    // Extend the last pixel for the rest of the reference values.
    nearest_pixel = top_border[num_ref_pixels - 1];
    for (; i < out_width - 1; ++i) {
      dst[i + 1] = nearest_pixel;
    }
  } else {
    // Extend nearest pixel.
    kvz_pixel nearest_pixel = x > 0 ? left_border[0] : dc_val;
    int i;
    for(i = 1; i < out_width; i++)
    {
      dst[i] = nearest_pixel;
    }
  }

  // If top-left corner sample doesn't exist, use the sample from below.
  // Unavailable samples on the left boundary are copied from below if
  // available. This is the only place they are available because we don't
  // support constrained intra prediction.
  if (x > 0 && y > 0) {
    // Make sure we always take the top-left pixel from the LCU reference
    // pixel arrays if they are available.
    if (x_local == 0) {
      dst[0] = left_border[-1];
    } else {
      dst[0] = top_border[-1];
    }
  } else {
    dst[0] = dst[dst_stride];
  }
}


void kvz_intra_build_reference(
  const int_fast8_t log2_width,
  const color_t color,
  const vector2d_t *const luma_px,
  const vector2d_t *const pic_px,
  const lcu_t *const lcu,
  kvz_intra_references *const refs)
{
  assert(log2_width >= 2 && log2_width <= 5);

  // Tables for looking up the number of intra reference pixels based on
  // prediction units coordinate within an LCU.
  // generated by "tools/generate_ref_pixel_tables.py".
  static const uint8_t num_ref_pixels_top[16][16] = {
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 }
  };
  static const uint8_t num_ref_pixels_left[16][16] = {
    { 64,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 60,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 56,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 52,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 48,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 44,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 40,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 36,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 },
    { 32,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 28,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 24,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 20,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 }
  };

  refs->filtered_initialized = false;
  kvz_pixel *out_left_ref = &refs->ref.left[0];
  kvz_pixel *out_top_ref = &refs->ref.top[0];

  const kvz_pixel dc_val = 1 << (KVZ_BIT_DEPTH - 1);
  const int is_chroma = color != COLOR_Y ? 1 : 0;
  const int_fast8_t width = 1 << log2_width;

  // Convert luma coordinates to chroma coordinates for chroma.
  const vector2d_t lcu_px = {
    luma_px->x % LCU_WIDTH,
    luma_px->y % LCU_WIDTH
  };
  const vector2d_t px = {
    lcu_px.x >> is_chroma,
    lcu_px.y >> is_chroma,
  };

  // Init pointers to LCUs reconstruction buffers, such that index 0 refers to block coordinate 0.
  const kvz_pixel *left_ref = !color ? &lcu->left_ref.y[1] : (color == 1) ? &lcu->left_ref.u[1] : &lcu->left_ref.v[1];
  const kvz_pixel *top_ref = !color ? &lcu->top_ref.y[1] : (color == 1) ? &lcu->top_ref.u[1] : &lcu->top_ref.v[1];
  const kvz_pixel *rec_ref = !color ? lcu->rec.y : (color == 1) ? lcu->rec.u : lcu->rec.v;

  // Init top borders pointer to point to the correct place in the correct reference array.
  const kvz_pixel *top_border;
  if (px.y) {
    top_border = &rec_ref[px.x + (px.y - 1) * (LCU_WIDTH >> is_chroma)];
  } else {
    top_border = &top_ref[px.x];
  }

  // Init left borders pointer to point to the correct place in the correct reference array.
  const kvz_pixel *left_border;
  int left_stride; // Distance between reference samples.
  if (px.x) {
    left_border = &rec_ref[px.x - 1 + px.y * (LCU_WIDTH >> is_chroma)];
    left_stride = LCU_WIDTH >> is_chroma;
  } else {
    left_border = &left_ref[px.y];
    left_stride = 1;
  }

  // Generate left reference.
  if (luma_px->x > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int px_available_left = num_ref_pixels_left[lcu_px.y / 4][lcu_px.x / 4] >> is_chroma;

    // Limit the number of available pixels based on block size and dimensions
    // of the picture.
    px_available_left = MIN(px_available_left, width * 2);
    px_available_left = MIN(px_available_left, (pic_px->y - luma_px->y) >> is_chroma);

    // Copy pixels from coded CUs.
    for (int i = 0; i < px_available_left; ++i) {
      out_left_ref[i + 1] = left_border[i * left_stride];
    }
    // Extend the last pixel for the rest of the reference values.
    kvz_pixel nearest_pixel = out_left_ref[px_available_left];
    for (int i = px_available_left; i < width * 2; ++i) {
      out_left_ref[i + 1] = nearest_pixel;
    }
  } else {
    // If we are on the left edge, extend the first pixel of the top row.
    kvz_pixel nearest_pixel = luma_px->y > 0 ? top_border[0] : dc_val;
    for (int i = 0; i < width * 2; i++) {
      out_left_ref[i + 1] = nearest_pixel;
    }
  }

  // Generate top-left reference.
  if (luma_px->x > 0 && luma_px->y > 0) {
    // If the block is at an LCU border, the top-left must be copied from
    // the border that points to the LCUs 1D reference buffer.
    if (px.x == 0) {
      out_left_ref[0] = left_border[-1 * left_stride];
      out_top_ref[0] = left_border[-1 * left_stride];
    } else {
      out_left_ref[0] = top_border[-1];
      out_top_ref[0] = top_border[-1];
    }
  } else {
    // Copy reference clockwise.
    out_left_ref[0] = out_left_ref[1];
    out_top_ref[0] = out_left_ref[1];
  }

  // Generate top reference.
  if (luma_px->y > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int px_available_top = num_ref_pixels_top[lcu_px.y / 4][lcu_px.x / 4] >> is_chroma;

    // Limit the number of available pixels based on block size and dimensions
    // of the picture.
    px_available_top = MIN(px_available_top, width * 2);
    px_available_top = MIN(px_available_top, (pic_px->x - luma_px->x) >> is_chroma);

    // Copy all the pixels we can.
    for (int i = 0; i < px_available_top; ++i) {
      out_top_ref[i + 1] = top_border[i];
    }
    // Extend the last pixel for the rest of the reference values.
    kvz_pixel nearest_pixel = top_border[px_available_top - 1];
    for (int i = px_available_top; i < width * 2; ++i) {
      out_top_ref[i + 1] = nearest_pixel;
    }
  } else {
    // Extend nearest pixel.
    kvz_pixel nearest_pixel = luma_px->x > 0 ? left_border[0] : dc_val;
    for (int i = 0; i < width * 2; i++) {
      out_top_ref[i + 1] = nearest_pixel;
    }
  }
}


/**
 * \brief this functions constructs the angular intra prediction from border samples
 *
 */
void kvz_intra_get_angular_pred(const encoder_control_t * const encoder, const kvz_pixel* src, int32_t src_stride, kvz_pixel* dst, int32_t dst_stride, int32_t width, int32_t dir_mode, int8_t filter)
{
  static const int32_t kvz_ang_table[9] = { 0, 2, 5, 9, 13, 17, 21, 26, 32 };
  static const int32_t kvz_inv_ang_table[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / Angle

  int32_t k,l;
  int32_t blk_size        = width;

  // Map the mode index to main prediction direction and angle
  bool mode_ver       = dir_mode >= 18;
  int32_t intra_pred_angle = mode_ver ? dir_mode - 26 : 10 - dir_mode;
  int32_t abs_ang       = abs(intra_pred_angle);
  int32_t sign_ang      = intra_pred_angle < 0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32_t inv_angle       = kvz_inv_ang_table[abs_ang];

  // Do angular predictions
  kvz_pixel *ref_main;
  kvz_pixel *ref_side;
  kvz_pixel  ref_above[2 * LCU_WIDTH + 1];
  kvz_pixel  ref_left[2 * LCU_WIDTH + 1];

  // Tell clang-analyzer that everything is ok.
  assert(width == 4 || width == 8 || width == 16 || width == 32);

  abs_ang           = kvz_ang_table[abs_ang];
  intra_pred_angle  = sign_ang * abs_ang;

  // Initialise the Main and Left reference array.
  if (intra_pred_angle < 0) {
    int32_t invAngleSum = 128; // rounding for (shift by 8)
    for (k = 0; k < blk_size + 1; k++) {
      ref_above[k + blk_size - 1] = src[k - src_stride - 1];
      ref_left[k + blk_size - 1]  = src[(k - 1) * src_stride - 1];
    }

    ref_main = (mode_ver ? ref_above : ref_left) + (blk_size - 1);
    ref_side = (mode_ver ? ref_left : ref_above) + (blk_size - 1);

    // Extend the Main reference to the left.
    for (k = -1; k > blk_size * intra_pred_angle>>5; k--) {
      invAngleSum += inv_angle;
      ref_main[k] = ref_side[invAngleSum>>8];
    }
  } else {
    for (k = 0; k < 2 * blk_size + 1; k++) {
      ref_above[k] = src[k - src_stride - 1];
      ref_left[k]  = src[(k - 1) * src_stride - 1];
    }
    ref_main = mode_ver ? ref_above : ref_left;
    ref_side = mode_ver ? ref_left  : ref_above;
  }

  if (intra_pred_angle == 0) {
    for (k = 0; k < blk_size; k++) {
      for (l = 0; l < blk_size; l++) {
        dst[k * dst_stride + l] = ref_main[l + 1];
      }
    }

    if (filter) {
      for (k=0;k<blk_size;k++) {
        dst[k * dst_stride] = CLIP(0, (1<<encoder->bitdepth) - 1, dst[k * dst_stride] + (( ref_side[k + 1] - ref_side[0]) >> 1));
      }
    }
  } else {
    int32_t delta_pos=0;
    int32_t delta_int;
    int32_t delta_fract;
    int32_t minus_delta_fract;
    int32_t ref_main_index;
    for (k = 0; k < blk_size; k++) {
      delta_pos += intra_pred_angle;
      delta_int   = delta_pos >> 5;
      delta_fract = delta_pos & (32 - 1);


      if (delta_fract) {
        minus_delta_fract = (32 - delta_fract);
        // Do linear filtering
        for (l = 0; l < blk_size; l++) {
          ref_main_index        = l + delta_int + 1;
          dst[k * dst_stride + l] = (kvz_pixel) ( (minus_delta_fract * ref_main[ref_main_index]
                                                 + delta_fract * ref_main[ref_main_index + 1] + 16) >> 5);
        }
      } else {
        // Just copy the integer samples
        for (l = 0; l < blk_size; l++) {
          dst[k * dst_stride + l] = ref_main[l + delta_int + 1];
        }
      }
    }
  }

  // Flip the block if this is the horizontal mode
  if (!mode_ver) {
    kvz_pixel tmp;
    for (k=0;k<blk_size-1;k++) {
      for (l=k+1;l<blk_size;l++) {
        tmp                 = dst[k * dst_stride + l];
        dst[k * dst_stride + l] = dst[l * dst_stride + k];
        dst[l * dst_stride + k] = tmp;
      }
    }
  }
}




void kvz_intra_dc_pred_filtering(const kvz_pixel *src, int32_t src_stride, kvz_pixel *dst, int32_t dst_stride, int32_t width, int32_t height )
{
  int32_t x, y, dst_stride2, src_stride2;

  // boundary pixels processing
  dst[0] = ((src[-src_stride] + src[-1] + 2 * dst[0] + 2) >> 2);

  for (x = 1; x < width; x++) {
    dst[x] = ((src[x - src_stride] +  3 * dst[x] + 2) >> 2);
  }
  for ( y = 1, dst_stride2 = dst_stride, src_stride2 = src_stride-1;
        y < height; y++, dst_stride2+=dst_stride, src_stride2+=src_stride ) {
    dst[dst_stride2] = ((src[src_stride2] + 3 * dst[dst_stride2] + 2) >> 2);
  }
  return;
}

/**
 * \brief Function for deriving planar intra prediction.
 * \param src source pixel array
 * \param srcstride source width
 * \param width block size to predict
 * \param dst destination buffer for prediction
 * \param dststride destination width

  This function derives the prediction samples for planar mode (intra coding).
*/
void kvz_intra_get_planar_pred(const kvz_pixel* src, int32_t srcstride, uint32_t width, kvz_pixel* dst, int32_t dststride)
{
  int32_t k, l, bottom_left, top_right;
  int32_t hor_pred;
  int32_t left_column[LCU_WIDTH+1], top_row[LCU_WIDTH+1], bottom_row[LCU_WIDTH+1], right_column[LCU_WIDTH+1];
  uint32_t blk_size = width;
  uint32_t offset_2d = width;
  uint32_t shift_1d = kvz_g_convert_to_bit[ width ] + 2;
  uint32_t shift_2d = shift_1d + 1;

  // Get left and above reference column and row
  for (k = 0; k < (int32_t)blk_size + 1; k++) {
    top_row[k] = src[k - srcstride];
    left_column[k] = src[k * srcstride - 1];
  }

  // Prepare intermediate variables used in interpolation
  bottom_left = left_column[blk_size];
  top_right   = top_row[blk_size];
  for (k = 0; k < (int32_t)blk_size; k++) {
    bottom_row[k]   = bottom_left - top_row[k];
    right_column[k] = top_right   - left_column[k];
    top_row[k]      <<= shift_1d;
    left_column[k]  <<= shift_1d;
  }

  // Generate prediction signal
  for (k = 0; k < (int32_t)blk_size; k++) {
    hor_pred = left_column[k] + offset_2d;
    for (l = 0; l < (int32_t)blk_size; l++) {
      hor_pred += right_column[k];
      top_row[l] += bottom_row[l];
      dst[k * dststride + l] = (kvz_pixel)((hor_pred + top_row[l]) >> shift_2d);
    }
  }
}

void kvz_intra_recon_lcu_luma(encoder_state_t * const state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const vector2d_t lcu_px = { x & 0x3f, y & 0x3f };
  if (cur_cu == NULL) {
    cur_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];
  }
  const int8_t width = LCU_WIDTH >> depth;

  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;

    kvz_intra_recon_lcu_luma(state, x,          y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x + offset, y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x,          y + offset, depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x + offset, y + offset, depth+1, intra_mode, NULL, lcu);

    if (depth < MAX_DEPTH) {
      cu_info_t *cu_a = &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset) >> 3) + (lcu_px.y >> 3)        *LCU_T_CU_WIDTH];
      cu_info_t *cu_b = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + ((lcu_px.y + offset) >> 3)*LCU_T_CU_WIDTH];
      cu_info_t *cu_c = &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset) >> 3) + ((lcu_px.y + offset) >> 3)*LCU_T_CU_WIDTH];
      if (cbf_is_set(cu_a->cbf.y, depth+1) || cbf_is_set(cu_b->cbf.y, depth+1) || cbf_is_set(cu_c->cbf.y, depth+1)) {
        cbf_set(&cur_cu->cbf.y, depth);
      }
    }

    return;
  }
  {
    const uint32_t pic_width = state->tile->frame->width;
    const uint32_t pic_height = state->tile->frame->height;

    // Pointers to reconstruction arrays
    kvz_pixel *recbase_y = &lcu->rec.y[lcu_px.x + lcu_px.y * LCU_WIDTH];

    kvz_pixel rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    kvz_pixel *rec_shift  = &rec[width * 2 + 8 + 1];

    int32_t rec_stride = LCU_WIDTH;

    kvz_intra_build_reference_border(encoder, x, y,(int16_t)width * 2 + 8, rec, (int16_t)width * 2 + 8, 0,
                                 pic_width, pic_height, lcu);
    kvz_intra_recon(encoder, rec_shift, width * 2 + 8,
                width, recbase_y, rec_stride, intra_mode, 0);

    kvz_quantize_lcu_luma_residual(state, x, y, depth, cur_cu, lcu);
  }
}

void kvz_intra_recon_lcu_chroma(encoder_state_t * const state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const vector2d_t lcu_px = { x & 0x3f, y & 0x3f };
  const int8_t width = LCU_WIDTH >> depth;
  const int8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  if (cur_cu == NULL) {
    cur_cu = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + (lcu_px.y >> 3)*LCU_T_CU_WIDTH];
  }

  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;

    kvz_intra_recon_lcu_chroma(state, x,          y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x + offset, y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x,          y + offset, depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x + offset, y + offset, depth+1, intra_mode, NULL, lcu);

    if (depth < MAX_DEPTH) {
      cu_info_t *cu_a = &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset) >> 3) + (lcu_px.y >> 3)        *LCU_T_CU_WIDTH];
      cu_info_t *cu_b = &lcu->cu[LCU_CU_OFFSET + (lcu_px.x >> 3) + ((lcu_px.y + offset) >> 3)*LCU_T_CU_WIDTH];
      cu_info_t *cu_c = &lcu->cu[LCU_CU_OFFSET + ((lcu_px.x + offset) >> 3) + ((lcu_px.y + offset) >> 3)*LCU_T_CU_WIDTH];
      if (cbf_is_set(cu_a->cbf.u, depth+1) || cbf_is_set(cu_b->cbf.u, depth+1) || cbf_is_set(cu_c->cbf.u, depth+1)) {
        cbf_set(&cur_cu->cbf.u, depth);
      }
      if (cbf_is_set(cu_a->cbf.v, depth+1) || cbf_is_set(cu_b->cbf.v, depth+1) || cbf_is_set(cu_c->cbf.v, depth+1)) {
        cbf_set(&cur_cu->cbf.v, depth);
      }
    }

    return;
  }

  {
    const uint32_t pic_width = state->tile->frame->width;
    const uint32_t pic_height = state->tile->frame->height;

    // Pointers to reconstruction arrays
    kvz_pixel *recbase_u = &lcu->rec.u[lcu_px.x/2 + (lcu_px.y * LCU_WIDTH)/4];
    kvz_pixel *recbase_v = &lcu->rec.v[lcu_px.x/2 + (lcu_px.y * LCU_WIDTH)/4];

    kvz_pixel rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];

    int32_t rec_stride = LCU_WIDTH;

    // Reconstruct chroma.
    if (!(x & 4 || y & 4)) {
      kvz_pixel *rec_shift_c  = &rec[width_c * 2 + 8 + 1];
      kvz_intra_build_reference_border(encoder, x, y,(int16_t)width_c * 2 + 8, rec, (int16_t)width_c * 2 + 8, 1,
                                   pic_width/2, pic_height/2, lcu);
      kvz_intra_recon(encoder,
                  rec_shift_c,
                  width_c * 2 + 8,
                  width_c,
                  recbase_u,
                  rec_stride >> 1,
                  intra_mode,
                  1);

      kvz_intra_build_reference_border(encoder, x, y,(int16_t)width_c * 2 + 8, rec, (int16_t)width_c * 2 + 8, 2,
                                   pic_width/2, pic_height/2, lcu);
      kvz_intra_recon(encoder,
                  rec_shift_c,
                  width_c * 2 + 8,
                  width_c,
                  recbase_v,
                  rec_stride >> 1,
                  intra_mode,
                  2);

      kvz_quantize_lcu_chroma_residual(state, x, y, depth, cur_cu, lcu);
    }
  }
}
