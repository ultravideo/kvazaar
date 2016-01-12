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

#include "intra.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "encoder.h"
#include "transform.h"
#include "strategies/strategies-intra.h"
#include "strategies/strategies-picture.h"


int8_t kvz_intra_get_dir_luma_predictor(
  const uint32_t x,
  const uint32_t y,
  int8_t *preds,
  const cu_info_t *const cur_cu,
  const cu_info_t *const left_cu,
  const cu_info_t *const above_cu)
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


static void intra_filter_reference(
  int_fast8_t log2_width,
  kvz_intra_references *refs)
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


static void intra_post_process_angular(
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
* \brief Generage planar prediction.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
*/
static void intra_pred_dc(
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
static void intra_pred_filtered_dc(
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


void kvz_intra_predict(
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
    static const int kvz_intra_hor_ver_dist_thres[5] = { 0, 7, 1, 0, 0 };
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
      intra_pred_filtered_dc(log2_width, used_ref->top, used_ref->left, dst);
    } else {
      intra_pred_dc(log2_width, used_ref->top, used_ref->left, dst);
    }
  } else {
    kvz_angular_pred(log2_width, mode, used_ref->top, used_ref->left, dst);
    if (color == COLOR_Y && width < 32) {
      if (mode == 10) {
        intra_post_process_angular(width, 1, used_ref->top, dst);
      } else if (mode == 26) {
        intra_post_process_angular(width, width, used_ref->left, dst);
      }
    }
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


void kvz_intra_recon_lcu_luma(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  int8_t intra_mode,
  cu_info_t *cur_cu,
  lcu_t *lcu)
{
  const vector2d_t lcu_px = { SUB_SCU(x), SUB_SCU(y) };
  if (cur_cu == NULL) {
    cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  }
  const int8_t width = LCU_WIDTH >> depth;

  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;

    kvz_intra_recon_lcu_luma(state, x,          y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x + offset, y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x,          y + offset, depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_luma(state, x + offset, y + offset, depth+1, intra_mode, NULL, lcu);

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

  // Perform intra prediction and put the result in correct place lcu.
  vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };
  vector2d_t luma_px = { x, y };
  kvz_intra_references refs;
  const int_fast8_t log2_width = kvz_g_convert_to_bit[width] + 2;
  kvz_intra_build_reference(log2_width, COLOR_Y, &luma_px, &pic_px, lcu, &refs);

  kvz_pixel pred[32 * 32];
  kvz_intra_predict(&refs, log2_width, intra_mode, COLOR_Y, pred);
  
  kvz_pixel *block_in_lcu = &lcu->rec.y[lcu_px.x + lcu_px.y * LCU_WIDTH];
  kvz_pixels_blit(pred, block_in_lcu, width, width, width, LCU_WIDTH);

  kvz_quantize_lcu_luma_residual(state, x, y, depth, cur_cu, lcu);
}


void kvz_intra_recon_lcu_chroma(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  int8_t intra_mode,
  cu_info_t *cur_cu,
  lcu_t *lcu)
{
  const vector2d_t lcu_px = { SUB_SCU(x), SUB_SCU(y) };
  const int8_t width = LCU_WIDTH >> depth;
  const int8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  if (cur_cu == NULL) {
    cur_cu = LCU_GET_CU_AT_PX(lcu, lcu_px.x, lcu_px.y);
  }

  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width / 2;

    kvz_intra_recon_lcu_chroma(state, x,          y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x + offset, y,          depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x,          y + offset, depth+1, intra_mode, NULL, lcu);
    kvz_intra_recon_lcu_chroma(state, x + offset, y + offset, depth+1, intra_mode, NULL, lcu);

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

  if (!(x & 4 || y & 4)) {
    const int_fast8_t log2_width_c = kvz_g_convert_to_bit[width_c] + 2;
    const vector2d_t luma_px = { x, y };
    const vector2d_t pic_px = { state->tile->frame->width, state->tile->frame->height };

    // Intra predict U-plane and put the result in lcu buffer.
    {
      kvz_intra_references refs;
      kvz_intra_build_reference(log2_width_c, COLOR_U, &luma_px, &pic_px, lcu, &refs);

      kvz_pixel pred[32 * 32];
      kvz_intra_predict(&refs, log2_width_c, intra_mode, COLOR_U, pred);

      kvz_pixel *pu_in_lcu = &lcu->rec.u[lcu_px.x / 2 + (lcu_px.y * LCU_WIDTH) / 4];
      kvz_pixels_blit(pred, pu_in_lcu, width_c, width_c, width_c, LCU_WIDTH_C);
    }

    // Intra predict V-plane and put the result in lcu buffer.
    {
      kvz_intra_references refs;
      kvz_intra_build_reference(log2_width_c, COLOR_V, &luma_px, &pic_px, lcu, &refs);
      
      kvz_pixel pred[32 * 32];
      kvz_intra_predict(&refs, log2_width_c, intra_mode, COLOR_V, pred);

      kvz_pixel *pu_in_lcu = &lcu->rec.v[lcu_px.x / 2 + (lcu_px.y * LCU_WIDTH) / 4];
      kvz_pixels_blit(pred, pu_in_lcu, width_c, width_c, width_c, LCU_WIDTH_C);
    }

    kvz_quantize_lcu_chroma_residual(state, x, y, depth, cur_cu, lcu);
  }
}
