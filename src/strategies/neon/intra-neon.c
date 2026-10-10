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

/**
 * \ingroup Optimization
 * \file
 * AArch64 NEON implementation of intra prediction (angular, planar, filtered DC).
 *
 * Bit-exact with the generic implementation for 8-bit video.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "strategies/neon/intra-neon.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-intra.h"
#include "strategyselector.h"

#include <stdlib.h>
#include <string.h>

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

static const int8_t modedisp2sampledisp[9] = { 0, 2, 5, 9, 13, 17, 21, 26, 32 };
static const int16_t modedisp2invsampledisp[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / sampledisp

/*
 * One row of an angular prediction with a non-zero fractional offset:
 *   out[x] = ((32 - f) * r[x] + f * r[x + 1] + 16) >> 5
 * The 16-bit sum is at most (32 * 255) + 16 which fits in uint16, and
 * vrshrn by 5 is exactly (sum + 16) >> 5.
 */
static inline void angular_filter_row(const kvz_pixel *r, int width, int f, kvz_pixel *out)
{
  const uint16_t wa = (uint16_t)(32 - f);
  const uint16_t wb = (uint16_t)f;
  int x = 0;
  for (; x + 8 <= width; x += 8) {
    uint16x8_t a = vmovl_u8(vld1_u8(r + x));
    uint16x8_t b = vmovl_u8(vld1_u8(r + x + 1));
    uint16x8_t s = vmulq_n_u16(a, wa);
    s = vmlaq_n_u16(s, b, wb);
    vst1_u8(out + x, vrshrn_n_u16(s, 5));
  }
  for (; x < width; ++x) {
    out[x] = (kvz_pixel)(((32 - f) * r[x] + f * r[x + 1] + 16) >> 5);
  }
}

static void angular_pred_neon(const int_fast8_t log2_width,
                              const int_fast8_t intra_mode,
                              const kvz_pixel *const in_ref_above,
                              const kvz_pixel *const in_ref_left,
                              kvz_pixel *const dst)
{
  const int_fast8_t width = 1 << log2_width;

  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 18;
  // Modes distance to horizontal or vertical mode.
  const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
  // Sample displacement per column in fractions of 32.
  const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

  // Temporary buffer for modes 11-25, indices from -width to width-1.
  kvz_pixel tmp_ref[2 * 32];
  const kvz_pixel *ref_main;

  if (sample_disp < 0) {
    // Negative sample_disp means, we need to use both references.
    const kvz_pixel *ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

    for (int_fast8_t x = -1; x < width; ++x) {
      tmp_ref[x + width] = ref_main[x];
    }
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
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
  }

  // Horizontal modes are computed in the vertical orientation and transposed.
  kvz_pixel tmp_out[32 * 32];
  kvz_pixel *buf = vertical_mode ? dst : tmp_out;

  if (sample_disp != 0) {
    int_fast16_t delta_pos = 0;
    for (int_fast8_t y = 0; y < width; ++y) {
      delta_pos += sample_disp;
      int_fast8_t delta_int = delta_pos >> 5;
      int_fast8_t delta_fract = delta_pos & (32 - 1);
      kvz_pixel *row = &buf[y * width];

      if (delta_fract) {
        angular_filter_row(ref_main + delta_int, width, delta_fract, row);
      }
      else {
        memcpy(row, ref_main + delta_int, (size_t)width);
      }
    }
  }
  else {
    for (int_fast8_t y = 0; y < width; ++y) {
      memcpy(&buf[y * width], ref_main, (size_t)width);
    }
  }

  if (!vertical_mode) {
    for (int_fast8_t y = 0; y < width; ++y) {
      for (int_fast8_t x = 0; x < width; ++x) {
        dst[y * width + x] = tmp_out[x * width + y];
      }
    }
  }
}

static void intra_pred_planar_neon(const int_fast8_t log2_width,
                                   const kvz_pixel *const ref_top,
                                   const kvz_pixel *const ref_left,
                                   kvz_pixel *const dst)
{
  const int width = 1 << log2_width;
  const int32_t top_right = ref_top[width + 1];
  const int32_t bottom_left = ref_left[width + 1];

  if (width < 8) {
    // Width 4 is too narrow for the 8-lane loop, so use the scalar form.
    for (int y = 0; y < width; ++y) {
      for (int x = 0; x < width; ++x) {
        int32_t hor = (width - 1 - x) * ref_left[y + 1] + (x + 1) * top_right;
        int32_t ver = (width - 1 - y) * ref_top[x + 1] + (y + 1) * bottom_left;
        dst[y * width + x] = (kvz_pixel)((ver + hor + width) >> (log2_width + 1));
      }
    }
    return;
  }

  // Per-column constants: t = top[x], a = width-1-x, b = x+1.
  int32_t tv[32], av[32], bv[32];
  for (int x = 0; x < width; ++x) {
    tv[x] = ref_top[x + 1];
    av[x] = width - 1 - x;
    bv[x] = x + 1;
  }

  const int32x4_t shift = vdupq_n_s32(-(log2_width + 1));
  const int32x4_t vtr = vdupq_n_s32(top_right);

  for (int y = 0; y < width; ++y) {
    const int32x4_t vl = vdupq_n_s32(ref_left[y + 1]);
    const int32x4_t vwy = vdupq_n_s32(width - 1 - y);
    const int32x4_t vadd = vdupq_n_s32((y + 1) * bottom_left + width);

    for (int x = 0; x < width; x += 8) {
      int32x4_t t0 = vld1q_s32(&tv[x]), t1 = vld1q_s32(&tv[x + 4]);
      int32x4_t a0 = vld1q_s32(&av[x]), a1 = vld1q_s32(&av[x + 4]);
      int32x4_t b0 = vld1q_s32(&bv[x]), b1 = vld1q_s32(&bv[x + 4]);

      int32x4_t acc0 = vmulq_s32(a0, vl);
      acc0 = vmlaq_s32(acc0, b0, vtr);
      acc0 = vmlaq_s32(acc0, t0, vwy);
      acc0 = vaddq_s32(acc0, vadd);

      int32x4_t acc1 = vmulq_s32(a1, vl);
      acc1 = vmlaq_s32(acc1, b1, vtr);
      acc1 = vmlaq_s32(acc1, t1, vwy);
      acc1 = vaddq_s32(acc1, vadd);

      acc0 = vshlq_s32(acc0, shift);
      acc1 = vshlq_s32(acc1, shift);

      vst1_u8(&dst[y * width + x],
              vqmovun_s16(vcombine_s16(vqmovn_s32(acc0), vqmovn_s32(acc1))));
    }
  }
}

static void intra_pred_filtered_dc_neon(const int_fast8_t log2_width,
                                        const kvz_pixel *const ref_top,
                                        const kvz_pixel *const ref_left,
                                        kvz_pixel *const out_block)
{
  const int width = 1 << log2_width;

  int sum = 0;
  for (int i = 0; i < width; ++i) {
    sum += ref_top[i + 1];
    sum += ref_left[i + 1];
  }
  const kvz_pixel dc_val = (kvz_pixel)((sum + width) >> (log2_width + 1));

  // Fill the whole block with the DC value, then patch the first row and column.
  for (int y = 0; y < width; ++y) {
    kvz_pixel *row = &out_block[y * width];
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      vst1_u8(row + x, vdup_n_u8(dc_val));
    }
    for (; x < width; ++x) {
      row[x] = dc_val;
    }
  }

  // Filter top-left with ([1 2 1] / 4)
  out_block[0] = (kvz_pixel)((ref_left[1] + 2 * dc_val + ref_top[1] + 2) / 4);

  // Filter rest of the top row with ([1 3] / 4)
  for (int x = 1; x < width; ++x) {
    out_block[x] = (kvz_pixel)((ref_top[x + 1] + 3 * dc_val + 2) / 4);
  }
  // Filter rest of the left column with ([1 3] / 4)
  for (int y = 1; y < width; ++y) {
    out_block[y * width] = (kvz_pixel)((ref_left[y + 1] + 3 * dc_val + 2) / 4);
  }
}

int kvz_strategy_register_intra_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "angular_pred", "neon", 40, &angular_pred_neon);
    success &= kvz_strategyselector_register(opaque, "intra_pred_planar", "neon", 40, &intra_pred_planar_neon);
    success &= kvz_strategyselector_register(opaque, "intra_pred_filtered_dc", "neon", 40, &intra_pred_filtered_dc_neon);
  }

  return success;
}

#else

int kvz_strategy_register_intra_neon(void* opaque, uint8_t bitdepth)
{
  (void)opaque;
  (void)bitdepth;
  return 1;
}

#endif
