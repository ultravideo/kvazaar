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

#include "strategies/generic/intra-generic.h"

#include <stdlib.h>

#include "kvazaar.h"
#include "strategyselector.h"


 /**
 * \brief Generage angular predictions.
 * \param log2_width    Log2 of width, range 2..5.
 * \param intra_mode    Angular mode in range 2..34.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_angular_pred_generic(
  const int_fast8_t log2_width,
  const int_fast8_t intra_mode,
  const kvz_pixel *const in_ref_above,
  const kvz_pixel *const in_ref_left,
  kvz_pixel *const dst)
{
  assert(log2_width >= 2 && log2_width <= 5);
  assert(intra_mode >= 2 && intra_mode <= 34);

  static const int8_t modedisp2sampledisp[9] = { 0, 2, 5, 9, 13, 17, 21, 26, 32 };
  static const int16_t modedisp2invsampledisp[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / sampledisp

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
  }
  else {
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
      }
      else {
        // Just copy the integer samples
        for (int_fast8_t x = 0; x < width; x++) {
          dst[y * width + x] = ref_main[x + delta_int];
        }
      }
    }
  }
  else {
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
 * \brief Generate planar prediction.
 * \param log2_width    Log2 of width, range 2..5.
 * \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
 * \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
 * \param dst           Buffer of size width*width.
 */
static void kvz_intra_pred_planar_generic(
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
* \brief Generage intra DC prediction with post filtering applied.
* \param log2_width    Log2 of width, range 2..5.
* \param in_ref_above  Pointer to -1 index of above reference, length=width*2+1.
* \param in_ref_left   Pointer to -1 index of left reference, length=width*2+1.
* \param dst           Buffer of size width*width.
*/
static void kvz_intra_pred_filtered_dc_generic(
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


int kvz_strategy_register_intra_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "angular_pred", "generic", 0, &kvz_angular_pred_generic);
  success &= kvz_strategyselector_register(opaque, "intra_pred_planar", "generic", 0, &kvz_intra_pred_planar_generic);
  success &= kvz_strategyselector_register(opaque, "intra_pred_filtered_dc", "generic", 0, &kvz_intra_pred_filtered_dc_generic);

  return success;
}
