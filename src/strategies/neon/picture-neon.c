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

#include "strategies/neon/picture-neon.h"

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

#include <stdlib.h>
#include <string.h>

#include "kvazaar.h"
#include "strategies/neon/neon-common.h"
#include "strategies/optimized_sad_func_ptr_t.h"
#include "strategies/strategies-picture.h"
#include "strategyselector.h"

/*
 * 8-bit NEON versions of the functions of picture-generic.c. Every function
 * computes exactly the same result as the generic version, including the
 * block-splitting order of the any-size SATD functions.
 */

//////////////////////////////////////////////////////////////////////////
// SAD

/**
 * Sum of absolute differences of a width x height block. The absolute
 * differences are accumulated in 32-bit lanes, so there is no overflow.
 */
static unsigned reg_sad_neon(const kvz_pixel *const data1, const kvz_pixel *const data2,
                             const int width, const int height,
                             const unsigned stride1, const unsigned stride2)
{
  if (width < 8) {
    /* Too narrow for a vector to pay for its setup and reduction. */
    unsigned sum = 0;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        sum += abs((int)data1[y * stride1 + x] - (int)data2[y * stride2 + x]);
      }
    }
    return sum;
  }

  uint32x4_t acc = vdupq_n_u32(0);
  unsigned tail = 0;

  for (int y = 0; y < height; ++y) {
    const kvz_pixel *p1 = data1 + y * stride1;
    const kvz_pixel *p2 = data2 + y * stride2;
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      uint8x16_t d = vabdq_u8(vld1q_u8(p1 + x), vld1q_u8(p2 + x));
      acc = vpadalq_u16(acc, vpaddlq_u8(d));
    }
    if (x + 8 <= width) {
      uint8x8_t d = vabd_u8(vld1_u8(p1 + x), vld1_u8(p2 + x));
      acc = vaddw_u16(acc, vpaddl_u8(d));
      x += 8;
    }
    for (; x < width; ++x) {
      tail += abs((int)p1[x] - (int)p2[x]);
    }
  }
  return vaddvq_u32(acc) + tail;
}

/* Dedicated 4-wide and 8-wide SAD. The generic row loop above handles narrow
 * blocks with a scalar tail, which made 4x4 and 8x8 slower than generic. Here a
 * 4-pixel row pair fills one 64-bit vector and an 8-pixel row fills one vector. */
static inline uint8x8_t load_row_pair4(const kvz_pixel *p, unsigned stride)
{
  uint32_t r0, r1;
  memcpy(&r0, p, 4);
  memcpy(&r1, p + stride, 4);
  return vreinterpret_u8_u64(vcreate_u64((uint64_t)r0 | ((uint64_t)r1 << 32)));
}

static inline unsigned sad_w4_neon(const kvz_pixel *p1, unsigned s1,
                                   const kvz_pixel *p2, unsigned s2, int height)
{
  uint16x8_t acc = vdupq_n_u16(0);
  unsigned tail = 0;
  int y = 0;
  for (; y + 2 <= height; y += 2) {
    acc = vaddw_u8(acc, vabd_u8(load_row_pair4(p1 + y * s1, s1),
                                load_row_pair4(p2 + y * s2, s2)));
  }
  for (; y < height; ++y) {
    for (int x = 0; x < 4; ++x) {
      tail += abs((int)p1[y * s1 + x] - (int)p2[y * s2 + x]);
    }
  }
  return vaddlvq_u16(acc) + tail;
}

static inline unsigned sad_w8_neon(const kvz_pixel *p1, unsigned s1,
                                   const kvz_pixel *p2, unsigned s2, int height)
{
  uint16x8_t acc = vdupq_n_u16(0);
  for (int y = 0; y < height; ++y) {
    acc = vaddw_u8(acc, vabd_u8(vld1_u8(p1 + y * s1), vld1_u8(p2 + y * s2)));
  }
  return vaddlvq_u16(acc);
}

static unsigned sad_4x4_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return sad_w4_neon(block1, 4, block2, 4, 4);
}
static unsigned sad_8x8_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return sad_w8_neon(block1, 8, block2, 8, 8);
}
static unsigned sad_16x16_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return reg_sad_neon(block1, block2, 16, 16, 16, 16);
}
static unsigned sad_32x32_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return reg_sad_neon(block1, block2, 32, 32, 32, 32);
}
static unsigned sad_64x64_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return reg_sad_neon(block1, block2, 64, 64, 64, 64);
}

static void sad_4x4_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = sad_4x4_neon(preds[0], orig);
  costs_out[1] = sad_4x4_neon(preds[1], orig);
}
static void sad_8x8_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = sad_8x8_neon(preds[0], orig);
  costs_out[1] = sad_8x8_neon(preds[1], orig);
}
static void sad_16x16_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = reg_sad_neon(preds[0], orig, 16, 16, 16, 16);
  costs_out[1] = reg_sad_neon(preds[1], orig, 16, 16, 16, 16);
}
static void sad_32x32_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = reg_sad_neon(preds[0], orig, 32, 32, 32, 32);
  costs_out[1] = reg_sad_neon(preds[1], orig, 32, 32, 32, 32);
}
static void sad_64x64_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = reg_sad_neon(preds[0], orig, 64, 64, 64, 64);
  costs_out[1] = reg_sad_neon(preds[1], orig, 64, 64, 64, 64);
}

/* Fixed-width SAD entry points for get_optimized_sad. The width is a
 * constant, so the compiler can unroll the row loop. */
#define REG_SAD_FIXED_W(w) \
static uint32_t reg_sad_w ## w ## _neon(const kvz_pixel * const data1, const kvz_pixel * const data2, \
                                        const int32_t height, const uint32_t stride1, const uint32_t stride2) \
{ \
  return reg_sad_neon(data1, data2, (w), height, stride1, stride2); \
}
static uint32_t reg_sad_w4_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                                const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return sad_w4_neon(data1, stride1, data2, stride2, height);
}
static uint32_t reg_sad_w8_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                                const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return sad_w8_neon(data1, stride1, data2, stride2, height);
}
REG_SAD_FIXED_W(12)
REG_SAD_FIXED_W(16)
REG_SAD_FIXED_W(24)
REG_SAD_FIXED_W(32)
REG_SAD_FIXED_W(64)

static uint32_t reg_sad_w0_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                                const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return 0;
}

static optimized_sad_func_ptr_t get_optimized_sad_neon(int32_t width)
{
  switch (width) {
    case 0:  return reg_sad_w0_neon;
    case 4:  return reg_sad_w4_neon;
    case 8:  return reg_sad_w8_neon;
    case 12: return reg_sad_w12_neon;
    case 16: return reg_sad_w16_neon;
    case 24: return reg_sad_w24_neon;
    case 32: return reg_sad_w32_neon;
    case 64: return reg_sad_w64_neon;
    default: return NULL;
  }
}

/* Vertically interpolated SAD: the reference row is shared by every row of
 * the picture block, so the second stride is zero. */
static uint32_t ver_sad_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                             int block_width, int block_height, unsigned pic_stride)
{
  return reg_sad_neon(pic_data, ref_data, block_width, block_height, pic_stride, 0);
}

/* SAD of a block against one reference sample per row, broadcast across the
 * row. */
static unsigned hor_sad_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                             int block_width, int block_height, unsigned pic_stride, unsigned ref_stride)
{
  uint32x4_t acc = vdupq_n_u32(0);
  unsigned tail = 0;

  for (int y = 0; y < block_height; ++y) {
    const kvz_pixel *p = pic_data + y * pic_stride;
    const kvz_pixel ref = ref_data[y * ref_stride];
    const uint8x16_t vref16 = vdupq_n_u8(ref);
    const uint8x8_t vref8 = vdup_n_u8(ref);
    int x = 0;
    for (; x + 16 <= block_width; x += 16) {
      uint8x16_t d = vabdq_u8(vld1q_u8(p + x), vref16);
      acc = vpadalq_u16(acc, vpaddlq_u8(d));
    }
    if (x + 8 <= block_width) {
      uint8x8_t d = vabd_u8(vld1_u8(p + x), vref8);
      acc = vaddw_u16(acc, vpaddl_u8(d));
      x += 8;
    }
    for (; x < block_width; ++x) {
      tail += abs((int)p[x] - (int)ref);
    }
  }
  return vaddvq_u32(acc) + tail;
}

static uint32_t hor_sad_neon_wrapper(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                     int32_t width, int32_t height, uint32_t pic_stride,
                                     uint32_t ref_stride, uint32_t left, uint32_t right)
{
  uint32_t result = 0;
  if (left) {
    result += hor_sad_neon(pic_data, ref_data + left, left,
                           height, pic_stride, ref_stride);

    result += reg_sad_neon(pic_data + left, ref_data + left, width - left,
                           height, pic_stride, ref_stride);
  } else if (right) {
    result += reg_sad_neon(pic_data, ref_data, width - right,
                           height, pic_stride, ref_stride);

    result += hor_sad_neon(pic_data + width - right,
                           ref_data + width - right - 1,
                           right, height, pic_stride, ref_stride);
  } else {
    result += reg_sad_neon(pic_data, ref_data, width,
                           height, pic_stride, ref_stride);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////
// SATD

/* Loads 8 pixels of a and b and returns a - b as signed 16-bit lanes. */
static INLINE int16x8_t load_diff8(const kvz_pixel *a, const kvz_pixel *b)
{
  return vreinterpretq_s16_u16(vsubl_u8(vld1_u8(a), vld1_u8(b)));
}

/* Loads 4 pixels of a and b and returns a - b in the low four 16-bit lanes. */
static INLINE int16x4_t load_diff4(const kvz_pixel *a, const kvz_pixel *b)
{
  uint32_t va, vb;
  memcpy(&va, a, 4);
  memcpy(&vb, b, 4);
  uint8x8_t ua = vcreate_u8((uint64_t)va);
  uint8x8_t ub = vcreate_u8((uint64_t)vb);
  return vget_low_s16(vreinterpretq_s16_u16(vsubl_u8(ua, ub)));
}

/*
 * One pass of the 4-point Hadamard transform over four rows. Each vector holds
 * one row, so the transform runs down the columns. The result matches the
 * first two stages of hadamard_4x4_generic.
 */
static INLINE void hadamard4_pass(int16x4_t r[4])
{
  int16x4_t s03p = vadd_s16(r[0], r[3]);
  int16x4_t s03m = vsub_s16(r[0], r[3]);
  int16x4_t s12p = vadd_s16(r[1], r[2]);
  int16x4_t s12m = vsub_s16(r[1], r[2]);
  r[0] = vadd_s16(s03p, s12p);
  r[1] = vadd_s16(s03m, s12m);
  r[2] = vsub_s16(s03p, s12p);
  r[3] = vsub_s16(s03m, s12m);
}

static INLINE void transpose4x4_s16(int16x4_t r[4])
{
  int16x4x2_t t01 = vtrn_s16(r[0], r[1]);
  int16x4x2_t t23 = vtrn_s16(r[2], r[3]);
  int32x2x2_t c02 = vtrn_s32(vreinterpret_s32_s16(t01.val[0]), vreinterpret_s32_s16(t23.val[0]));
  int32x2x2_t c13 = vtrn_s32(vreinterpret_s32_s16(t01.val[1]), vreinterpret_s32_s16(t23.val[1]));
  r[0] = vreinterpret_s16_s32(c02.val[0]);
  r[1] = vreinterpret_s16_s32(c13.val[0]);
  r[2] = vreinterpret_s16_s32(c02.val[1]);
  r[3] = vreinterpret_s16_s32(c13.val[1]);
}

/*
 * One pass of the 8-point Hadamard transform over eight rows, the same
 * butterfly stages as the vertical pass of satd_8x8_subblock_generic.
 */
static INLINE void hadamard8_pass(int16x8_t m[8])
{
  int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  m[0] = vaddq_s16(b0, b1); m[1] = vsubq_s16(b0, b1);
  m[2] = vaddq_s16(b2, b3); m[3] = vsubq_s16(b2, b3);
  m[4] = vaddq_s16(b4, b5); m[5] = vsubq_s16(b4, b5);
  m[6] = vaddq_s16(b6, b7); m[7] = vsubq_s16(b6, b7);
}

/* 8x8 transpose of 16-bit lanes, built from 16-, 32- and 64-bit lane swaps. */
static INLINE void transpose8x8_s16(int16x8_t r[8])
{
  int16x8x2_t t01 = vtrnq_s16(r[0], r[1]);
  int16x8x2_t t23 = vtrnq_s16(r[2], r[3]);
  int16x8x2_t t45 = vtrnq_s16(r[4], r[5]);
  int16x8x2_t t67 = vtrnq_s16(r[6], r[7]);

  int32x4x2_t u02 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[0]), vreinterpretq_s32_s16(t23.val[0]));
  int32x4x2_t u13 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[1]), vreinterpretq_s32_s16(t23.val[1]));
  int32x4x2_t u46 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[0]), vreinterpretq_s32_s16(t67.val[0]));
  int32x4x2_t u57 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[1]), vreinterpretq_s32_s16(t67.val[1]));

  int16x8_t v00 = vreinterpretq_s16_s32(u02.val[0]);
  int16x8_t v01 = vreinterpretq_s16_s32(u02.val[1]);
  int16x8_t v10 = vreinterpretq_s16_s32(u13.val[0]);
  int16x8_t v11 = vreinterpretq_s16_s32(u13.val[1]);
  int16x8_t v20 = vreinterpretq_s16_s32(u46.val[0]);
  int16x8_t v21 = vreinterpretq_s16_s32(u46.val[1]);
  int16x8_t v30 = vreinterpretq_s16_s32(u57.val[0]);
  int16x8_t v31 = vreinterpretq_s16_s32(u57.val[1]);

  r[0] = vcombine_s16(vget_low_s16(v00), vget_low_s16(v20));
  r[4] = vcombine_s16(vget_high_s16(v00), vget_high_s16(v20));
  r[2] = vcombine_s16(vget_low_s16(v01), vget_low_s16(v21));
  r[6] = vcombine_s16(vget_high_s16(v01), vget_high_s16(v21));
  r[1] = vcombine_s16(vget_low_s16(v10), vget_low_s16(v30));
  r[5] = vcombine_s16(vget_high_s16(v10), vget_high_s16(v30));
  r[3] = vcombine_s16(vget_low_s16(v11), vget_low_s16(v31));
  r[7] = vcombine_s16(vget_high_s16(v11), vget_high_s16(v31));
}

/* SATD of one 4x4 block. Matches hadamard_4x4_generic, including the final
 * halving with rounding. */
static unsigned satd_4x4_subblock_neon(const kvz_pixel *buf1, int32_t stride1,
                                       const kvz_pixel *buf2, int32_t stride2)
{
  int16x4_t r[4];
  for (int i = 0; i < 4; ++i) {
    r[i] = load_diff4(buf1 + i * stride1, buf2 + i * stride2);
  }
  hadamard4_pass(r);
  transpose4x4_s16(r);
  hadamard4_pass(r);

  int32_t sum = 0;
  for (int i = 0; i < 4; ++i) {
    sum += vaddlv_s16(vabs_s16(r[i]));
  }
  return (sum + 1) >> 1;
}

/* SATD of one 8x8 block, matching satd_8x8_subblock_generic. */
static unsigned satd_8x8_subblock_neon(const kvz_pixel *org, int32_t stride_org,
                                       const kvz_pixel *cur, int32_t stride_cur)
{
  int16x8_t m[8];
  for (int i = 0; i < 8; ++i) {
    m[i] = load_diff8(org + i * stride_org, cur + i * stride_cur);
  }
  hadamard8_pass(m);
  transpose8x8_s16(m);
  hadamard8_pass(m);

  int32_t sum = 0;
  for (int i = 0; i < 8; ++i) {
    sum += vaddlvq_s16(vabsq_s16(m[i]));
  }
  return (sum + 2) >> 2;
}

static unsigned satd_4x4_neon(const kvz_pixel *piOrg, const kvz_pixel *piCur)
{
  return satd_4x4_subblock_neon(piOrg, 4, piCur, 4);
}

static void satd_4x4_subblock_quad_neon(const kvz_pixel *preds[4], const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned costs[4])
{
  for (int k = 0; k < 4; ++k) {
    costs[k] = satd_4x4_subblock_neon(orig, orig_stride, preds[k], stride);
  }
}

static void satd_8x8_subblock_quad_neon(const kvz_pixel **preds, const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned *costs)
{
  for (int k = 0; k < 4; ++k) {
    costs[k] = satd_8x8_subblock_neon(orig, orig_stride, preds[k], stride);
  }
}

#define SATD_NXN_NEON(n) \
static unsigned satd_ ## n ## x ## n ## _neon(const kvz_pixel * const block1, \
                                              const kvz_pixel * const block2) \
{ \
  unsigned sum = 0; \
  for (unsigned y = 0; y < (n); y += 8) { \
    unsigned row = y * (n); \
    for (unsigned x = 0; x < (n); x += 8) { \
      sum += satd_8x8_subblock_neon(&block1[row + x], (n), &block2[row + x], (n)); \
    } \
  } \
  return sum; \
} \
\
static void satd_ ## n ## x ## n ## _dual_neon(const pred_buffer preds, const kvz_pixel * const orig, \
                                               unsigned num_modes, unsigned *costs_out) \
{ \
  for (int k = 0; k < 2; ++k) { \
    unsigned sum = 0; \
    for (unsigned y = 0; y < (n); y += 8) { \
      unsigned row = y * (n); \
      for (unsigned x = 0; x < (n); x += 8) { \
        sum += satd_8x8_subblock_neon(&preds[k][row + x], (n), &orig[row + x], (n)); \
      } \
    } \
    costs_out[k] = sum; \
  } \
}

SATD_NXN_NEON(8)
SATD_NXN_NEON(16)
SATD_NXN_NEON(32)
SATD_NXN_NEON(64)

static void satd_4x4_dual_neon(const pred_buffer preds, const kvz_pixel * const orig,
                               unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = satd_4x4_neon(orig, preds[0]);
  costs_out[1] = satd_4x4_neon(orig, preds[1]);
}

static unsigned satd_any_size_neon(int width, int height,
                                   const kvz_pixel *block1, int stride1,
                                   const kvz_pixel *block2, int stride2)
{
  unsigned sum = 0;
  if (width % 8 != 0) {
    /* Process the first column using 4x4 blocks. */
    for (int y = 0; y < height; y += 4) {
      sum += satd_4x4_subblock_neon(&block1[y * stride1], stride1,
                                    &block2[y * stride2], stride2);
    }
    block1 += 4;
    block2 += 4;
    width -= 4;
  }
  if (height % 8 != 0) {
    /* Process the first row using 4x4 blocks. */
    for (int x = 0; x < width; x += 4) {
      sum += satd_4x4_subblock_neon(&block1[x], stride1,
                                    &block2[x], stride2);
    }
    block1 += 4 * stride1;
    block2 += 4 * stride2;
    height -= 4;
  }
  /* The rest can now be processed with 8x8 blocks. */
  for (int y = 0; y < height; y += 8) {
    const kvz_pixel *row1 = &block1[y * stride1];
    const kvz_pixel *row2 = &block2[y * stride2];
    for (int x = 0; x < width; x += 8) {
      sum += satd_8x8_subblock_neon(&row1[x], stride1, &row2[x], stride2);
    }
  }
  return sum;
}

/* Quad version of satd_any_size. It follows satd_any_size_quad_generic step by
 * step, so the results are identical for every input the generic accepts. */
static void satd_any_size_quad_neon(int width, int height,
                                    const kvz_pixel **preds, const int stride,
                                    const kvz_pixel *orig, const int orig_stride,
                                    unsigned num_modes, unsigned *costs_out, int8_t *valid)
{
  unsigned sums[4] = { 0 };
  const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };
  const kvz_pixel *orig_ptr = orig;
  costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0;

  const int width_mod_8 = width % 8;
  if (width_mod_8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    /* Process the first column using 4x4 blocks. */
    for (int y = 0; y < height; y += 4) {
      satd_4x4_subblock_quad_neon(pred_ptrs_tmp, stride, &orig[y * orig_stride], orig_stride, sums);
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs_tmp[blk] += 4 * stride;
        costs_out[blk] += sums[blk];
      }
    }
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] += 4;
    }
    width -= 4;
  }
  if (height % 8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    /* Process the first row using 4x4 blocks. */
    for (int x = 0; x < width; x += 4) {
      satd_4x4_subblock_quad_neon(pred_ptrs_tmp, stride, &orig_ptr[x], orig_stride, sums);
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs_tmp[blk] += 4;
        costs_out[blk] += sums[blk];
      }
    }
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] += 4 * stride;
    }
    height -= 4;
  }
  /* The rest can now be processed with 8x8 blocks. */
  for (int y = height % 8; y < height; y += 8) {
    orig_ptr = &orig[y * orig_stride + width_mod_8];
    pred_ptrs[0] = &preds[0][y * stride + width_mod_8];
    pred_ptrs[1] = &preds[1][y * stride + width_mod_8];
    pred_ptrs[2] = &preds[2][y * stride + width_mod_8];
    pred_ptrs[3] = &preds[3][y * stride + width_mod_8];
    for (int x = width_mod_8; x < width; x += 8) {
      satd_8x8_subblock_quad_neon(pred_ptrs, stride, orig_ptr, orig_stride, sums);
      orig_ptr += 8;
      pred_ptrs[0] += 8;
      pred_ptrs[1] += 8;
      pred_ptrs[2] += 8;
      pred_ptrs[3] += 8;
      costs_out[0] += sums[0];
      costs_out[1] += sums[1];
      costs_out[2] += sums[2];
      costs_out[3] += sums[3];
    }
  }
}

//////////////////////////////////////////////////////////////////////////
// SSD

static unsigned pixels_calc_ssd_neon(const kvz_pixel *const ref, const kvz_pixel *const rec,
                                     const int ref_stride, const int rec_stride,
                                     const int width)
{
  int32x4_t acc = vdupq_n_s32(0);
  int tail = 0;

  for (int y = 0; y < width; ++y) {
    const kvz_pixel *r = ref + y * ref_stride;
    const kvz_pixel *c = rec + y * rec_stride;
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int16x8_t d = load_diff8(r + x, c + x);
      acc = vmlal_s16(acc, vget_low_s16(d), vget_low_s16(d));
      acc = vmlal_high_s16(acc, d, d);
    }
    for (; x < width; ++x) {
      int diff = (int)r[x] - (int)c[x];
      tail += diff * diff;
    }
  }
  return (unsigned)(vaddvq_s32(acc) + tail);
}

//////////////////////////////////////////////////////////////////////////
// Bi-prediction

/* Rounds two 16-bit inputs with the bi-prediction offset and shift, and
 * saturates the result to 8 bits. Each input is widened to 32 bits first, so
 * the sum cannot overflow. */
static INLINE uint8x8_t bipred_round8(int16x8_t a, int16x8_t b)
{
  const int32x4_t off = vdupq_n_s32(1 << 6);
  int32x4_t lo = vaddq_s32(vaddl_s16(vget_low_s16(a), vget_low_s16(b)), off);
  int32x4_t hi = vaddq_s32(vaddl_high_s16(a, b), off);
  lo = vshrq_n_s32(lo, 7);
  hi = vshrq_n_s32(hi, 7);
  return vqmovun_s16(vcombine_s16(vqmovn_s32(lo), vqmovn_s32(hi)));
}

/* Converts 8-bit samples to the 14-bit intermediate precision. */
static INLINE int16x8_t px_to_im8(const kvz_pixel *p)
{
  return vreinterpretq_s16_u16(vshlq_n_u16(vmovl_u8(vld1_u8(p)), 6));
}

static void bipred_average_px_px(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = px_to_im8(px_L0 + y * pu_w + x);
      int16x8_t b = px_to_im8(px_L1 + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_L0 = px_L0[y * pu_w + x] << 6;
      int16_t sample_L1 = px_L1[y * pu_w + x] << 6;
      int32_t rounded = (sample_L0 + sample_L1 + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_im_im(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = vld1q_s16(im_L0 + y * pu_w + x);
      int16x8_t b = vld1q_s16(im_L1 + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_L0 = im_L0[y * pu_w + x];
      int16_t sample_L1 = im_L1[y * pu_w + x];
      int32_t rounded = (sample_L0 + sample_L1 + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_px_im(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = px_to_im8(px + y * pu_w + x);
      int16x8_t b = vld1q_s16(im + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_px = px[y * pu_w + x] << 6;
      int16_t sample_im = im[y * pu_w + x];
      int32_t rounded = (sample_px + sample_im + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_neon(lcu_t *const lcu,
  const yuv_t *const px_L0,
  const yuv_t *const px_L1,
  const yuv_im_t *const im_L0,
  const yuv_im_t *const im_L1,
  const unsigned pu_x,
  const unsigned pu_y,
  const unsigned pu_w,
  const unsigned pu_h,
  const unsigned im_flags_L0,
  const unsigned im_flags_L1,
  const bool predict_luma,
  const bool predict_chroma,
  const uint8_t chroma_shift_w,
  const uint8_t chroma_shift_h)
{
  // After reconstruction, merge the predictors by taking an average of each pixel.
  if (predict_luma) {
    unsigned pb_offset = SUB_SCU(pu_y) * LCU_WIDTH + SUB_SCU(pu_x);

    if (!(im_flags_L0 & 1) && !(im_flags_L1 & 1)) {
      bipred_average_px_px(lcu->rec.y + pb_offset, px_L0->y, px_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else if ((im_flags_L0 & 1) && (im_flags_L1 & 1)) {
      bipred_average_im_im(lcu->rec.y + pb_offset, im_L0->y, im_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else {
      kvz_pixel    *src_px = (im_flags_L0 & 1) ? px_L1->y : px_L0->y;
      kvz_pixel_im *src_im = (im_flags_L0 & 1) ? im_L0->y : im_L1->y;
      bipred_average_px_im(lcu->rec.y + pb_offset, src_px, src_im, pu_w, pu_h, LCU_WIDTH);
    }
  }
  if (predict_chroma) {
    unsigned pb_offset = ((SUB_SCU(pu_y) >> chroma_shift_h) * (LCU_WIDTH >> chroma_shift_w)) + (SUB_SCU(pu_x) >> chroma_shift_w);
    unsigned pb_w = pu_w >> chroma_shift_w;
    unsigned pb_h = pu_h >> chroma_shift_h;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_px_px(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_px(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_im_im(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_im_im(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    } else {
      kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_px_im(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_im(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    }
  }
}

//////////////////////////////////////////////////////////////////////////
// Variance

/*
 * Variance of an array of 8-bit samples. The sums are exact integers, and the
 * variance is n * sum(x^2) - (sum x)^2 over n^2, which equals the two-pass
 * definition of pixel_var_generic.
 */
static double pixel_var_neon(const kvz_pixel *arr, const uint32_t len)
{
  uint32x4_t sum_acc = vdupq_n_u32(0);
  uint32x4_t sq_acc = vdupq_n_u32(0);
  uint32_t sum = 0;
  uint32_t sumsq = 0;
  uint32_t i = 0;

  for (; i + 16 <= len; i += 16) {
    uint8x16_t v = vld1q_u8(arr + i);
    sum_acc = vpadalq_u16(sum_acc, vpaddlq_u8(v));
    uint16x8_t sq_lo = vmull_u8(vget_low_u8(v), vget_low_u8(v));
    uint16x8_t sq_hi = vmull_u8(vget_high_u8(v), vget_high_u8(v));
    sq_acc = vpadalq_u16(sq_acc, sq_lo);
    sq_acc = vpadalq_u16(sq_acc, sq_hi);
  }
  sum = vaddvq_u32(sum_acc);
  sumsq = vaddvq_u32(sq_acc);
  for (; i < len; ++i) {
    sum += arr[i];
    sumsq += (uint32_t)arr[i] * arr[i];
  }

  const double n = (double)len;
  return (n * (double)sumsq - (double)sum * (double)sum) / (n * n);
}

#endif // COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

int kvz_strategy_register_picture_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "reg_sad", "neon", 40, &reg_sad_neon);

    success &= kvz_strategyselector_register(opaque, "sad_4x4", "neon", 40, &sad_4x4_neon);
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "neon", 40, &sad_8x8_neon);
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "neon", 40, &sad_16x16_neon);
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "neon", 40, &sad_32x32_neon);
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "neon", 40, &sad_64x64_neon);

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "neon", 40, &satd_4x4_neon);
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "neon", 40, &satd_8x8_neon);
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "neon", 40, &satd_16x16_neon);
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "neon", 40, &satd_32x32_neon);
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "neon", 40, &satd_64x64_neon);

    success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "neon", 40, &sad_4x4_dual_neon);
    success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "neon", 40, &sad_8x8_dual_neon);
    success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "neon", 40, &sad_16x16_dual_neon);
    success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "neon", 40, &sad_32x32_dual_neon);
    success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "neon", 40, &sad_64x64_dual_neon);

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "neon", 40, &satd_4x4_dual_neon);
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "neon", 40, &satd_8x8_dual_neon);
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "neon", 40, &satd_16x16_dual_neon);
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "neon", 40, &satd_32x32_dual_neon);
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "neon", 40, &satd_64x64_dual_neon);
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "neon", 40, &satd_any_size_neon);
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "neon", 40, &satd_any_size_quad_neon);

    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "neon", 40, &pixels_calc_ssd_neon);
    success &= kvz_strategyselector_register(opaque, "bipred_average", "neon", 40, &bipred_average_neon);

    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "neon", 40, &get_optimized_sad_neon);
    success &= kvz_strategyselector_register(opaque, "ver_sad", "neon", 40, &ver_sad_neon);
    success &= kvz_strategyselector_register(opaque, "hor_sad", "neon", 40, &hor_sad_neon_wrapper);

    success &= kvz_strategyselector_register(opaque, "pixel_var", "neon", 40, &pixel_var_neon);
  }
#endif
  return success;
}
