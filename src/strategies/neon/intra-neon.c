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
  /* The weights fit in 8 bits, so widening u8 multiplies replace the separate
   * widening moves. */
  const uint8x8_t wa = vdup_n_u8((uint8_t)(32 - f));
  const uint8x8_t wb = vdup_n_u8((uint8_t)f);
  int x = 0;
  for (; x + 16 <= width; x += 16) {
    const uint8x16_t a = vld1q_u8(r + x);
    const uint8x16_t b = vld1q_u8(r + x + 1);
    const uint16x8_t lo = vmlal_u8(vmull_u8(vget_low_u8(a), wa), vget_low_u8(b), wb);
    const uint16x8_t hi = vmlal_high_u8(vmull_high_u8(a, vcombine_u8(wa, wa)), b, vcombine_u8(wb, wb));
    vst1q_u8(out + x, vrshrn_high_n_u16(vrshrn_n_u16(lo, 5), hi, 5));
  }
  for (; x + 8 <= width; x += 8) {
    const uint16x8_t sum = vmlal_u8(vmull_u8(vld1_u8(r + x), wa), vld1_u8(r + x + 1), wb);
    vst1_u8(out + x, vrshrn_n_u16(sum, 5));
  }
  for (; x < width; ++x) {
    out[x] = (kvz_pixel)(((32 - f) * r[x] + f * r[x + 1] + 16) >> 5);
  }
}

/* Two 4-wide angular rows at once: lanes 0-3 are row y (offset d0, fraction f0)
 * and lanes 4-7 are row y + 1. Only r[d .. d + 4] is read, like generic. A
 * fraction of 0 gives ((32 * r[x] + 16) >> 5) == r[x], the plain copy. */
static inline uint8x8_t angular_rows_w4(const kvz_pixel *r, int d0, int f0, int d1, int f1)
{
  uint32_t a0, a1, b0, b1;
  memcpy(&a0, r + d0, 4);
  memcpy(&b0, r + d0 + 1, 4);
  memcpy(&a1, r + d1, 4);
  memcpy(&b1, r + d1 + 1, 4);
  const uint8x8_t a = vreinterpret_u8_u64(vcreate_u64((uint64_t)a0 | ((uint64_t)a1 << 32)));
  const uint8x8_t b = vreinterpret_u8_u64(vcreate_u64((uint64_t)b0 | ((uint64_t)b1 << 32)));
  const uint16x4_t wa0 = vdup_n_u16((uint16_t)(32 - f0)), wa1 = vdup_n_u16((uint16_t)(32 - f1));
  const uint16x4_t wb0 = vdup_n_u16((uint16_t)f0), wb1 = vdup_n_u16((uint16_t)f1);
  uint16x8_t s = vmulq_u16(vmovl_u8(a), vcombine_u16(wa0, wa1));
  s = vmlaq_u16(s, vmovl_u8(b), vcombine_u16(wb0, wb1));
  return vrshrn_n_u16(s, 5);
}

/* In-register transpose of an 8x8 byte block. */
static inline void transpose8x8_u8(uint8x8_t r[8])
{
  const uint8x8x2_t t0 = vtrn_u8(r[0], r[1]);
  const uint8x8x2_t t1 = vtrn_u8(r[2], r[3]);
  const uint8x8x2_t t2 = vtrn_u8(r[4], r[5]);
  const uint8x8x2_t t3 = vtrn_u8(r[6], r[7]);
  const uint16x4x2_t u0 = vtrn_u16(vreinterpret_u16_u8(t0.val[0]), vreinterpret_u16_u8(t1.val[0]));
  const uint16x4x2_t u1 = vtrn_u16(vreinterpret_u16_u8(t0.val[1]), vreinterpret_u16_u8(t1.val[1]));
  const uint16x4x2_t u2 = vtrn_u16(vreinterpret_u16_u8(t2.val[0]), vreinterpret_u16_u8(t3.val[0]));
  const uint16x4x2_t u3 = vtrn_u16(vreinterpret_u16_u8(t2.val[1]), vreinterpret_u16_u8(t3.val[1]));
  const uint32x2x2_t v0 = vtrn_u32(vreinterpret_u32_u16(u0.val[0]), vreinterpret_u32_u16(u2.val[0]));
  const uint32x2x2_t v1 = vtrn_u32(vreinterpret_u32_u16(u1.val[0]), vreinterpret_u32_u16(u3.val[0]));
  const uint32x2x2_t v2 = vtrn_u32(vreinterpret_u32_u16(u0.val[1]), vreinterpret_u32_u16(u2.val[1]));
  const uint32x2x2_t v3 = vtrn_u32(vreinterpret_u32_u16(u1.val[1]), vreinterpret_u32_u16(u3.val[1]));
  r[0] = vreinterpret_u8_u32(v0.val[0]);
  r[4] = vreinterpret_u8_u32(v0.val[1]);
  r[1] = vreinterpret_u8_u32(v1.val[0]);
  r[5] = vreinterpret_u8_u32(v1.val[1]);
  r[2] = vreinterpret_u8_u32(v2.val[0]);
  r[6] = vreinterpret_u8_u32(v2.val[1]);
  r[3] = vreinterpret_u8_u32(v3.val[0]);
  r[7] = vreinterpret_u8_u32(v3.val[1]);
}

/* dst[y * width + x] = src[x * width + y] for a square block. */
static inline void transpose_block(const kvz_pixel *src, kvz_pixel *dst, int width)
{
  if (width == 4) {
    static const uint8_t idx[16] = { 0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15 };
    vst1q_u8(dst, vqtbl1q_u8(vld1q_u8(src), vld1q_u8(idx)));
    return;
  }
  for (int by = 0; by < width; by += 8) {
    for (int bx = 0; bx < width; bx += 8) {
      uint8x8_t r[8];
      for (int k = 0; k < 8; ++k) {
        r[k] = vld1_u8(src + (by + k) * width + bx);
      }
      transpose8x8_u8(r);
      for (int k = 0; k < 8; ++k) {
        vst1_u8(dst + (bx + k) * width + by, r[k]);
      }
    }
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

  if (sample_disp != 0 && width == 4) {
    int_fast16_t delta_pos = 0;
    for (int_fast8_t y = 0; y < 4; y += 2) {
      const int p0 = delta_pos + sample_disp;
      const int p1 = p0 + sample_disp;
      delta_pos = p1;
      vst1_u8(&buf[y * 4], angular_rows_w4(ref_main, p0 >> 5, p0 & 31, p1 >> 5, p1 & 31));
    }
  }
  else if (sample_disp != 0) {
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
    // Pure horizontal or vertical: every row is the reference row.
    if (width >= 16) {
      const uint8x16_t r0 = vld1q_u8(ref_main);
      const uint8x16_t r1 = width == 32 ? vld1q_u8(ref_main + 16) : r0;
      for (int_fast8_t y = 0; y < width; ++y) {
        vst1q_u8(&buf[y * width], r0);
        if (width == 32) vst1q_u8(&buf[y * width + 16], r1);
      }
    } else if (width == 8) {
      const uint8x8_t r0 = vld1_u8(ref_main);
      for (int_fast8_t y = 0; y < 8; ++y) vst1_u8(&buf[y * 8], r0);
    } else {
      uint32_t r0;
      memcpy(&r0, ref_main, 4);
      vst1q_u8(buf, vreinterpretq_u8_u32(vdupq_n_u32(r0)));
    }
  }

  if (!vertical_mode) {
    transpose_block(tmp_out, dst, width);
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
  const kvz_pixel *const top = ref_top + 1;
  const kvz_pixel *const left = ref_left + 1;

  // DC value: the rounded mean of the top and left references.
  uint32_t sum;
  if (width == 4) {
    uint32_t t, l;
    memcpy(&t, top, 4);
    memcpy(&l, left, 4);
    sum = vaddlv_u8(vreinterpret_u8_u64(vcreate_u64((uint64_t)t | ((uint64_t)l << 32))));
  } else {
    uint16x8_t acc = vdupq_n_u16(0);
    for (int i = 0; i < width; i += 8) {
      acc = vaddw_u8(acc, vld1_u8(top + i));
      acc = vaddw_u8(acc, vld1_u8(left + i));
    }
    sum = vaddlvq_u16(acc);
  }
  const kvz_pixel dc_val = (kvz_pixel)((sum + width) >> (log2_width + 1));

  // Fill the block with the DC value. width * width is a multiple of 16.
  const uint8x16_t fill = vdupq_n_u8(dc_val);
  for (int i = 0; i < width * width; i += 16) {
    vst1q_u8(out_block + i, fill);
  }

  // Top row: (top[x] + 3 * dc + 2) >> 2, in 16-bit lanes.
  const uint16x8_t bias = vdupq_n_u16((uint16_t)(3 * dc_val + 2));
  if (width == 4) {
    uint32_t t;
    memcpy(&t, top, 4);
    const uint8x8_t r = vshrn_n_u16(vaddw_u8(bias, vreinterpret_u8_u32(vdup_n_u32(t))), 2);
    vst1_lane_u32((uint32_t *)(void *)out_block, vreinterpret_u32_u8(r), 0);
  } else {
    for (int x = 0; x < width; x += 8) {
      vst1_u8(out_block + x, vshrn_n_u16(vaddw_u8(bias, vld1_u8(top + x)), 2));
    }
  }
  // Left column, same filter.
  for (int y = 1; y < width; ++y) {
    out_block[y * width] = (kvz_pixel)((left[y] + 3 * dc_val + 2) >> 2);
  }
  // Top-left corner with ([1 2 1] / 4).
  out_block[0] = (kvz_pixel)((left[0] + 2 * dc_val + top[0] + 2) >> 2);
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
