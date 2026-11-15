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

#include "strategies/strategies-cpp.h"

extern "C" {
#include "strategies/neon/picture-neon.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "strategies/neon/neon-common.h"
#include "strategies/optimized_sad_func_ptr_t.h"
#include "strategies/strategies-picture.h"
#include "strategyselector.h"
}

#if COMPILE_ARM_NEON
extern "C" {
#include "strategies/generic/picture-generic.h"
}

/*
 * SATD strategies for 10 and 12-bit pixels and temporal filter strategies for
 * all bit depths. The results are identical to the generic functions.
 */
namespace {

#if KVZ_BIT_DEPTH > 8
//////////////////////////////////////////////////////////////////////////
// SATD, 10 and 12-bit
//
// The 2D Hadamard transform of an 8x8 block multiplies the magnitude by 8 in
// each direction. The differences have at most KVZ_BIT_DEPTH + 1 bits, so the
// first pass fits 16-bit lanes at 10 and 12 bits (at most 32760). The last
// butterfly stage feeds straight into the absolute sum and
// |a + b| + |a - b| == 2 * max(|a|, |b|) for integers, so it is replaced by a
// maximum. At 10 bits the two remaining stages of the second pass also fit 16
// bits (at most 32736). At 12 bits they are done in 32-bit lanes. In a 4x4
// block every value before the last stage fits 16 bits (at most 32760).

static INLINE int16x8_t hbd_load_diff8(const kvz_pixel *a, const kvz_pixel *b)
{
  return vreinterpretq_s16_u16(vsubq_u16(vld1q_u16(a), vld1q_u16(b)));
}

static INLINE int16x4_t hbd_load_diff4(const kvz_pixel *a, const kvz_pixel *b)
{
  return vreinterpret_s16_u16(vsub_u16(vld1_u16(a), vld1_u16(b)));
}

/* One pass of the 4-point transform over the four rows r, for two 4x4 blocks
 * at a time (one in each 64-bit half). */
static INLINE void hbd_hadamard4_pass_x2(int16x8_t r[4])
{
  const int16x8_t s03p = vaddq_s16(r[0], r[3]);
  const int16x8_t s03m = vsubq_s16(r[0], r[3]);
  const int16x8_t s12p = vaddq_s16(r[1], r[2]);
  const int16x8_t s12m = vsubq_s16(r[1], r[2]);
  r[0] = vaddq_s16(s03p, s12p);
  r[1] = vaddq_s16(s03m, s12m);
  r[2] = vsubq_s16(s03p, s12p);
  r[3] = vsubq_s16(s03m, s12m);
}

/* Transposes the 4x4 block in each 64-bit half. */
static INLINE void hbd_transpose4x4_x2(int16x8_t r[4])
{
  const int16x8x2_t t01 = vtrnq_s16(r[0], r[1]);
  const int16x8x2_t t23 = vtrnq_s16(r[2], r[3]);
  const int32x4x2_t c02 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[0]), vreinterpretq_s32_s16(t23.val[0]));
  const int32x4x2_t c13 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[1]), vreinterpretq_s32_s16(t23.val[1]));
  r[0] = vreinterpretq_s16_s32(c02.val[0]);
  r[1] = vreinterpretq_s16_s32(c13.val[0]);
  r[2] = vreinterpretq_s16_s32(c02.val[1]);
  r[3] = vreinterpretq_s16_s32(c13.val[1]);
}

/*
 * SATD of two 4x4 blocks: r[i] holds row i of the differences of the first
 * block in the low half and of the second block in the high half. Returns
 * half of the absolute sum of each block, which is what hadamard_4x4_generic
 * returns.
 */
static INLINE uint64x2_t hbd_satd_4x4_x2(int16x8_t r[4])
{
  hbd_hadamard4_pass_x2(r);
  hbd_transpose4x4_x2(r);
  const int16x8_t s03p = vaddq_s16(r[0], r[3]);
  const int16x8_t s03m = vsubq_s16(r[0], r[3]);
  const int16x8_t s12p = vaddq_s16(r[1], r[2]);
  const int16x8_t s12m = vsubq_s16(r[1], r[2]);
  const uint16x8_t mp = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(s03p), vabsq_s16(s12p)));
  const uint16x8_t mm = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(s03m), vabsq_s16(s12m)));
  // Each maximum is at most 32760, so their sum fits 16 unsigned bits.
  return vpaddlq_u32(vpaddlq_u16(vaddq_u16(mp, mm)));
}

static INLINE unsigned hbd_satd_4x4(int16x4_t d0, int16x4_t d1, int16x4_t d2, int16x4_t d3)
{
  // Two copies of the block, as the transpose works on 64-bit halves.
  int16x8_t r[4] = { vcombine_s16(d0, d0), vcombine_s16(d1, d1), vcombine_s16(d2, d2), vcombine_s16(d3, d3) };
  return (unsigned)vgetq_lane_u64(hbd_satd_4x4_x2(r), 0);
}

static unsigned satd_4x4_subblock_hbd_neon(const kvz_pixel *buf1, int32_t stride1,
                                           const kvz_pixel *buf2, int32_t stride2)
{
  return hbd_satd_4x4(hbd_load_diff4(buf1, buf2),
                      hbd_load_diff4(buf1 + stride1, buf2 + stride2),
                      hbd_load_diff4(buf1 + 2 * stride1, buf2 + 2 * stride2),
                      hbd_load_diff4(buf1 + 3 * stride1, buf2 + 3 * stride2));
}

/* Sum of the SATDs of two 4x4 blocks, the second at (dx, dy) from the first. */
static unsigned satd_4x4_subblock_pair_hbd_neon(const kvz_pixel *buf1, int32_t stride1,
                                                const kvz_pixel *buf2, int32_t stride2,
                                                int dx, int dy)
{
  const int o1 = dy * stride1 + dx, o2 = dy * stride2 + dx;
  int16x8_t r[4];
  for (int i = 0; i < 4; ++i) {
    r[i] = vcombine_s16(hbd_load_diff4(buf1 + i * stride1, buf2 + i * stride2),
                        hbd_load_diff4(buf1 + o1 + i * stride1, buf2 + o2 + i * stride2));
  }
  const uint64x2_t s = hbd_satd_4x4_x2(r);
  return (unsigned)(vgetq_lane_u64(s, 0) + vgetq_lane_u64(s, 1));
}

static void satd_4x4_subblock_quad_hbd_neon(const kvz_pixel *preds[4], const int stride,
                                            const kvz_pixel *orig, const int orig_stride,
                                            unsigned costs[4])
{
  uint16x4_t o[4];
  for (int i = 0; i < 4; ++i) o[i] = vld1_u16(orig + i * orig_stride);
  for (int k = 0; k < 4; k += 2) {
    int16x8_t r[4];
    for (int i = 0; i < 4; ++i) {
      const uint16x8_t oo = vcombine_u16(o[i], o[i]);
      const uint16x8_t p = vcombine_u16(vld1_u16(preds[k] + i * stride), vld1_u16(preds[k + 1] + i * stride));
      r[i] = vreinterpretq_s16_u16(vsubq_u16(oo, p));
    }
    const uint64x2_t s = hbd_satd_4x4_x2(r);
    costs[k] = (unsigned)vgetq_lane_u64(s, 0);
    costs[k + 1] = (unsigned)vgetq_lane_u64(s, 1);
  }
}

/* One pass of the 8-point transform over the eight rows m. */
static INLINE void hbd_hadamard8_pass(int16x8_t m[8])
{
  const int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  const int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  const int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  const int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  const int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  const int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  const int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  const int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  m[0] = vaddq_s16(b0, b1); m[1] = vsubq_s16(b0, b1);
  m[2] = vaddq_s16(b2, b3); m[3] = vsubq_s16(b2, b3);
  m[4] = vaddq_s16(b4, b5); m[5] = vsubq_s16(b4, b5);
  m[6] = vaddq_s16(b6, b7); m[7] = vsubq_s16(b6, b7);
}

static INLINE void hbd_transpose8x8(int16x8_t r[8])
{
  const int16x8x2_t t01 = vtrnq_s16(r[0], r[1]);
  const int16x8x2_t t23 = vtrnq_s16(r[2], r[3]);
  const int16x8x2_t t45 = vtrnq_s16(r[4], r[5]);
  const int16x8x2_t t67 = vtrnq_s16(r[6], r[7]);

  const int32x4x2_t u02 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[0]), vreinterpretq_s32_s16(t23.val[0]));
  const int32x4x2_t u13 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[1]), vreinterpretq_s32_s16(t23.val[1]));
  const int32x4x2_t u46 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[0]), vreinterpretq_s32_s16(t67.val[0]));
  const int32x4x2_t u57 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[1]), vreinterpretq_s32_s16(t67.val[1]));

  const int16x8_t v00 = vreinterpretq_s16_s32(u02.val[0]);
  const int16x8_t v01 = vreinterpretq_s16_s32(u02.val[1]);
  const int16x8_t v10 = vreinterpretq_s16_s32(u13.val[0]);
  const int16x8_t v11 = vreinterpretq_s16_s32(u13.val[1]);
  const int16x8_t v20 = vreinterpretq_s16_s32(u46.val[0]);
  const int16x8_t v21 = vreinterpretq_s16_s32(u46.val[1]);
  const int16x8_t v30 = vreinterpretq_s16_s32(u57.val[0]);
  const int16x8_t v31 = vreinterpretq_s16_s32(u57.val[1]);

  r[0] = vcombine_s16(vget_low_s16(v00), vget_low_s16(v20));
  r[4] = vcombine_s16(vget_high_s16(v00), vget_high_s16(v20));
  r[2] = vcombine_s16(vget_low_s16(v01), vget_low_s16(v21));
  r[6] = vcombine_s16(vget_high_s16(v01), vget_high_s16(v21));
  r[1] = vcombine_s16(vget_low_s16(v10), vget_low_s16(v30));
  r[5] = vcombine_s16(vget_high_s16(v10), vget_high_s16(v30));
  r[3] = vcombine_s16(vget_low_s16(v11), vget_low_s16(v31));
  r[7] = vcombine_s16(vget_high_s16(v11), vget_high_s16(v31));
}

#if KVZ_BIT_DEPTH <= 10
/* The second pass fused with the absolute sum, in 16-bit lanes. Returns half
 * of the absolute sum. */
static INLINE uint32_t hbd_hadamard8_pass_absmax_sum(const int16x8_t m[8])
{
  const int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  const int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  const int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  const int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  const int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  const int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  const int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  const int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  const uint16x8_t m01 = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(b0), vabsq_s16(b1)));
  const uint16x8_t m23 = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(b2), vabsq_s16(b3)));
  const uint16x8_t m45 = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(b4), vabsq_s16(b5)));
  const uint16x8_t m67 = vreinterpretq_u16_s16(vmaxq_s16(vabsq_s16(b6), vabsq_s16(b7)));
  // Each maximum is at most 32736, so the sum of two fits 16 unsigned bits.
  return vaddvq_u32(vpadalq_u16(vpaddlq_u16(vaddq_u16(m01, m23)), vaddq_u16(m45, m67)));
}
#else
/* The last two stages of the second pass fused with the absolute sum, in
 * 32-bit lanes, for four columns. a holds the outputs of the first stage. */
static INLINE uint32x4_t hbd_ver8_absmax_32(const int32x4_t a[8])
{
  const int32x4_t b0 = vaddq_s32(a[0], a[2]), b2 = vsubq_s32(a[0], a[2]);
  const int32x4_t b1 = vaddq_s32(a[1], a[3]), b3 = vsubq_s32(a[1], a[3]);
  const int32x4_t b4 = vaddq_s32(a[4], a[6]), b6 = vsubq_s32(a[4], a[6]);
  const int32x4_t b5 = vaddq_s32(a[5], a[7]), b7 = vsubq_s32(a[5], a[7]);
  const int32x4_t m01 = vmaxq_s32(vabsq_s32(b0), vabsq_s32(b1));
  const int32x4_t m23 = vmaxq_s32(vabsq_s32(b2), vabsq_s32(b3));
  const int32x4_t m45 = vmaxq_s32(vabsq_s32(b4), vabsq_s32(b5));
  const int32x4_t m67 = vmaxq_s32(vabsq_s32(b6), vabsq_s32(b7));
  return vreinterpretq_u32_s32(vaddq_s32(vaddq_s32(m01, m23), vaddq_s32(m45, m67)));
}

/* The second pass fused with the absolute sum. The first stage widens to 32
 * bits. Returns half of the absolute sum. */
static INLINE uint32_t hbd_hadamard8_pass_absmax_sum(const int16x8_t m[8])
{
  int32x4_t lo[8], hi[8];
  for (int i = 0; i < 4; ++i) {
    lo[i] = vaddl_s16(vget_low_s16(m[i]), vget_low_s16(m[i + 4]));
    lo[i + 4] = vsubl_s16(vget_low_s16(m[i]), vget_low_s16(m[i + 4]));
    hi[i] = vaddl_high_s16(m[i], m[i + 4]);
    hi[i + 4] = vsubl_high_s16(m[i], m[i + 4]);
  }
  return vaddvq_u32(vaddq_u32(hbd_ver8_absmax_32(lo), hbd_ver8_absmax_32(hi)));
}
#endif

/* SATD of one 8x8 block given its differences, as satd_8x8_subblock_generic. */
static INLINE unsigned hbd_satd_8x8(int16x8_t m[8])
{
  hbd_hadamard8_pass(m);
  hbd_transpose8x8(m);
  // (sum + 2) >> 2 with sum = 2 * half.
  return (hbd_hadamard8_pass_absmax_sum(m) + 1) >> 1;
}

static unsigned satd_8x8_subblock_hbd_neon(const kvz_pixel *buf1, int32_t stride1,
                                           const kvz_pixel *buf2, int32_t stride2)
{
  int16x8_t m[8];
  for (int i = 0; i < 8; ++i) {
    m[i] = hbd_load_diff8(buf1 + i * stride1, buf2 + i * stride2);
  }
  return hbd_satd_8x8(m);
}

static void satd_8x8_subblock_quad_hbd_neon(const kvz_pixel **preds, const int stride,
                                            const kvz_pixel *orig, const int orig_stride,
                                            unsigned *costs)
{
  uint16x8_t o[8];
  for (int i = 0; i < 8; ++i) o[i] = vld1q_u16(orig + i * orig_stride);
  for (int k = 0; k < 4; ++k) {
    int16x8_t m[8];
    for (int i = 0; i < 8; ++i) {
      m[i] = vreinterpretq_s16_u16(vsubq_u16(o[i], vld1q_u16(preds[k] + i * stride)));
    }
    costs[k] = hbd_satd_8x8(m);
  }
}

static unsigned satd_4x4_hbd_neon(const kvz_pixel *org, const kvz_pixel *cur)
{
  return satd_4x4_subblock_hbd_neon(org, 4, cur, 4);
}

static void satd_4x4_dual_hbd_neon(const pred_buffer preds, const kvz_pixel *const orig,
                                   unsigned /* num_modes */, unsigned *costs_out)
{
  int16x8_t r[4];
  for (int i = 0; i < 4; ++i) {
    const uint16x4_t o = vld1_u16(orig + 4 * i);
    r[i] = vreinterpretq_s16_u16(vsubq_u16(vcombine_u16(o, o),
                                           vcombine_u16(vld1_u16(preds[0] + 4 * i), vld1_u16(preds[1] + 4 * i))));
  }
  const uint64x2_t s = hbd_satd_4x4_x2(r);
  costs_out[0] = (unsigned)vgetq_lane_u64(s, 0);
  costs_out[1] = (unsigned)vgetq_lane_u64(s, 1);
}

/* Sum of the 8x8 SATDs of an area of width x height pixels, both multiples of 8. */
static INLINE unsigned satd_8x8_blocks_hbd_neon(int width, int height,
                                                const kvz_pixel *buf1, int stride1,
                                                const kvz_pixel *buf2, int stride2)
{
  unsigned sum = 0;
  for (int y = 0; y < height; y += 8) {
    for (int x = 0; x < width; x += 8) {
      sum += satd_8x8_subblock_hbd_neon(buf1 + y * stride1 + x, stride1, buf2 + y * stride2 + x, stride2);
    }
  }
  return sum;
}

template <int N>
unsigned satd_nxn_hbd_neon(const kvz_pixel *const block1, const kvz_pixel *const block2)
{
  return satd_8x8_blocks_hbd_neon(N, N, block1, N, block2, N) >> (KVZ_BIT_DEPTH - 8);
}

template <int N>
void satd_nxn_dual_hbd_neon(const pred_buffer preds, const kvz_pixel *const orig,
                            unsigned /* num_modes */, unsigned *costs_out)
{
  costs_out[0] = satd_8x8_blocks_hbd_neon(N, N, preds[0], N, orig, N) >> (KVZ_BIT_DEPTH - 8);
  costs_out[1] = satd_8x8_blocks_hbd_neon(N, N, preds[1], N, orig, N) >> (KVZ_BIT_DEPTH - 8);
}

/* Follows the SATD_ANY_SIZE macro of strategies-picture.h. */
static unsigned satd_any_size_hbd_neon(int width, int height,
                                       const kvz_pixel *block1, int stride1,
                                       const kvz_pixel *block2, int stride2)
{
  unsigned sum = 0;
  if (width % 8 != 0) {
    // Process the first column using 4x4 blocks, two at a time.
    int y = 0;
    for (; y + 8 <= height; y += 8) {
      sum += satd_4x4_subblock_pair_hbd_neon(&block1[y * stride1], stride1,
                                             &block2[y * stride2], stride2, 0, 4);
    }
    for (; y < height; y += 4) {
      sum += satd_4x4_subblock_hbd_neon(&block1[y * stride1], stride1,
                                        &block2[y * stride2], stride2);
    }
    block1 += 4;
    block2 += 4;
    width -= 4;
  }
  if (height % 8 != 0) {
    // Process the first row using 4x4 blocks, two at a time.
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      sum += satd_4x4_subblock_pair_hbd_neon(&block1[x], stride1, &block2[x], stride2, 4, 0);
    }
    for (; x < width; x += 4) {
      sum += satd_4x4_subblock_hbd_neon(&block1[x], stride1, &block2[x], stride2);
    }
    block1 += 4 * stride1;
    block2 += 4 * stride2;
    height -= 4;
  }
  // The rest can now be processed with 8x8 blocks.
  sum += satd_8x8_blocks_hbd_neon(width, height, block1, stride1, block2, stride2);
  return sum >> (KVZ_BIT_DEPTH - 8);
}

/* Follows satd_any_size_quad_generic step by step, so the results are
 * identical for every input the generic accepts. */
static void satd_any_size_quad_hbd_neon(int width, int height,
                                        const kvz_pixel **preds, const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned /* num_modes */, unsigned *costs_out, int8_t * /* valid */)
{
  unsigned sums[4] = { 0 };
  const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };
  const kvz_pixel *orig_ptr = orig;
  costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0;
  const int width_mod_8 = width % 8;
  if (width_mod_8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    // Process the first column using 4x4 blocks.
    for (int y = 0; y < height; y += 4) {
      satd_4x4_subblock_quad_hbd_neon(pred_ptrs_tmp, stride, &orig[y * orig_stride], orig_stride, sums);
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
    // Process the first row using 4x4 blocks.
    for (int x = 0; x < width; x += 4) {
      satd_4x4_subblock_quad_hbd_neon(pred_ptrs_tmp, stride, &orig_ptr[x], orig_stride, sums);
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
  // The rest can now be processed with 8x8 blocks.
  for (int y = height % 8; y < height; y += 8) {
    orig_ptr = &orig[y * orig_stride + width_mod_8];
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] = &preds[blk][y * stride + width_mod_8];
    }
    for (int x = width_mod_8; x < width; x += 8) {
      satd_8x8_subblock_quad_hbd_neon(pred_ptrs, stride, orig_ptr, orig_stride, sums);
      orig_ptr += 8;
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs[blk] += 8;
        costs_out[blk] += sums[blk];
      }
    }
  }
  for (int blk = 0; blk < 4; ++blk) {
    costs_out[blk] >>= (KVZ_BIT_DEPTH - 8);
  }
}
#endif // KVZ_BIT_DEPTH > 8

//////////////////////////////////////////////////////////////////////////
// Temporal filter, all bit depths

/* Loads 8 or 4 pixels as unsigned 16-bit lanes. */
static INLINE uint16x8_t tf_load8_neon(const kvz_pixel *p)
{
#if KVZ_BIT_DEPTH == 8
  return vmovl_u8(vld1_u8(p));
#else
  return vld1q_u16(p);
#endif
}

static INLINE uint16x4_t tf_load4_neon(const kvz_pixel *p)
{
#if KVZ_BIT_DEPTH == 8
  uint32_t v;
  memcpy(&v, p, 4);
  return vget_low_u16(vmovl_u8(vcreate_u8((uint64_t)v)));
#else
  return vld1_u16(p);
#endif
}

/*
 * Quarter-pel bilinear interpolation of 8 pixels. The weights sum to 16, so
 * the result before the shift is at most 16 * 4095 + 8: the arithmetic is
 * exact in unsigned 16-bit lanes.
 */
static INLINE uint16x8_t tf_bilinear8_neon(uint16x8_t a, uint16x8_t b, uint16x8_t c, uint16x8_t d,
                                           uint16_t w00, uint16_t w10, uint16_t w01, uint16_t w11)
{
  uint16x8_t p = vmulq_n_u16(a, w00);
  p = vmlaq_n_u16(p, b, w10);
  p = vmlaq_n_u16(p, c, w01);
  p = vmlaq_n_u16(p, d, w11);
  return vrshrq_n_u16(p, 4);  // (p + 8) >> 4
}

/*
 * Sum of squared differences of a size x size block (a multiple of 8) and the
 * bilinear interpolation of ref with the fraction (fx, fy), or ref itself
 * when both are 0. The sum wraps like the unsigned sum of the generic code.
 */
static uint32_t tf_bilinear_ssd_block_neon(const kvz_pixel *org, int org_stride,
                                           const kvz_pixel *ref, int ref_stride,
                                           int size, int fx, int fy)
{
  const uint16_t w00 = (uint16_t)((4 - fx) * (4 - fy)), w10 = (uint16_t)(fx * (4 - fy));
  const uint16_t w01 = (uint16_t)((4 - fx) * fy), w11 = (uint16_t)(fx * fy);
  const bool interp = fx || fy;
  uint32x4_t acc = vdupq_n_u32(0);
  for (int y = 0; y < size; y++) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    const kvz_pixel *r1 = r0 + ref_stride;
    const kvz_pixel *o = &org[y * org_stride];
    for (int x = 0; x < size; x += 8) {
      uint16x8_t p = tf_load8_neon(r0 + x);
      if (interp) {
        p = tf_bilinear8_neon(p, tf_load8_neon(r0 + x + 1), tf_load8_neon(r1 + x), tf_load8_neon(r1 + x + 1),
                              w00, w10, w01, w11);
      }
      const uint16x8_t ad = vabdq_u16(tf_load8_neon(o + x), p);
      acc = vmlal_u16(acc, vget_low_u16(ad), vget_low_u16(ad));
      acc = vmlal_high_u16(acc, ad, ad);
    }
  }
  return vaddvq_u32(acc);
}

static unsigned bilinear_ssd_neon(const kvz_pixel *org, int org_stride,
                                  const kvz_pixel *ref, int ref_stride,
                                  int size, int fx, int fy)
{
  if (size % 8 == 0) {
    return tf_bilinear_ssd_block_neon(org, org_stride, ref, ref_stride, size, fx, fy);
  }
  // The generic code, for other block sizes.
  const int w00 = (4 - fx) * (4 - fy), w10 = fx * (4 - fy), w01 = (4 - fx) * fy, w11 = fx * fy;
  unsigned ssd = 0;
  for (int y = 0; y < size; y++) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    const kvz_pixel *r1 = r0 + ref_stride;
    const kvz_pixel *o = &org[y * org_stride];
    for (int x = 0; x < size; x++) {
      const int pred = (w00 * r0[x] + w10 * r0[x + 1] + w01 * r1[x] + w11 * r1[x + 1] + 8) >> 4;
      const int diff = o[x] - pred;
      ssd += diff * diff;
    }
  }
  return ssd;
}

static void bilinear_ssd_multi_neon(const kvz_pixel *org, int org_stride,
                                    const kvz_pixel *ref, int ref_stride,
                                    int size, const int *mv_x, const int *mv_y,
                                    int num, unsigned *out)
{
  if (size % 8 != 0) {
    kvz_bilinear_ssd_multi_generic(org, org_stride, ref, ref_stride, size, mv_x, mv_y, num, out);
    return;
  }
  for (int i = 0; i < num; i++) {
    const int fx = mv_x[i] & 3, fy = mv_y[i] & 3;
    const kvz_pixel *r = ref + (mv_y[i] >> 2) * ref_stride + (mv_x[i] >> 2);
    const uint32_t ssd = tf_bilinear_ssd_block_neon(org, org_stride, r, ref_stride, size, fx, fy);
    // Integer positions are scaled to 8 bits like kvz_pixels_calc_ssd.
    out[i] = (fx || fy) ? ssd : ssd >> (2 * (KVZ_BIT_DEPTH - 8));
  }
}

/*
 * Bilinear interpolation of 8 (or the low 4) pixels with weights that sum to
 * 1 << shift (at most 64). At 8 bits the sums fit unsigned 16-bit lanes, at
 * 10 and 12 bits they are computed in 32-bit lanes. The rounding shift is
 * (v + (1 << (shift - 1))) >> shift, and the narrowing truncates like the
 * cast of the generic code.
 */
#if KVZ_BIT_DEPTH == 8
typedef uint8x8_t tf_px8_t;
static INLINE tf_px8_t tf_bilinear_pred8_neon(uint16x8_t a, uint16x8_t b, uint16x8_t c, uint16x8_t d,
                                              const uint16_t w[4], int16x8_t neg_shift)
{
  uint16x8_t p = vmulq_n_u16(a, w[0]);
  p = vmlaq_n_u16(p, b, w[1]);
  p = vmlaq_n_u16(p, c, w[2]);
  p = vmlaq_n_u16(p, d, w[3]);
  return vmovn_u16(vrshlq_u16(p, neg_shift));
}
#else
typedef uint16x8_t tf_px8_t;
static INLINE uint16x4_t tf_bilinear_pred4_neon(uint16x4_t a, uint16x4_t b, uint16x4_t c, uint16x4_t d,
                                                const uint16_t w[4], int32x4_t neg_shift)
{
  uint32x4_t p = vmull_n_u16(a, w[0]);
  p = vmlal_n_u16(p, b, w[1]);
  p = vmlal_n_u16(p, c, w[2]);
  p = vmlal_n_u16(p, d, w[3]);
  return vmovn_u32(vrshlq_u32(p, neg_shift));
}

static INLINE tf_px8_t tf_bilinear_pred8_neon(uint16x8_t a, uint16x8_t b, uint16x8_t c, uint16x8_t d,
                                              const uint16_t w[4], int32x4_t neg_shift)
{
  return vcombine_u16(tf_bilinear_pred4_neon(vget_low_u16(a), vget_low_u16(b), vget_low_u16(c), vget_low_u16(d),
                                             w, neg_shift),
                      tf_bilinear_pred4_neon(vget_high_u16(a), vget_high_u16(b), vget_high_u16(c), vget_high_u16(d),
                                             w, neg_shift));
}
#endif

static void tf_bilinear_pred_neon(kvz_pixel *dst, int dst_stride,
                                  const kvz_pixel *ref, int ref_stride,
                                  int width, int height,
                                  int w00, int w10, int w01, int w11, int shift)
{
  const uint16_t w[4] = { (uint16_t)w00, (uint16_t)w10, (uint16_t)w01, (uint16_t)w11 };
#if KVZ_BIT_DEPTH == 8
  const int16x8_t neg_shift = vdupq_n_s16((int16_t)-shift);
#else
  const int32x4_t neg_shift = vdupq_n_s32(-shift);
#endif
  const int round = 1 << (shift - 1);
  for (int y = 0; y < height; y++) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    const kvz_pixel *r1 = r0 + ref_stride;
    kvz_pixel *d = &dst[y * dst_stride];
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      const tf_px8_t p = tf_bilinear_pred8_neon(tf_load8_neon(r0 + x), tf_load8_neon(r0 + x + 1),
                                                tf_load8_neon(r1 + x), tf_load8_neon(r1 + x + 1), w, neg_shift);
#if KVZ_BIT_DEPTH == 8
      vst1_u8(d + x, p);
#else
      vst1q_u16(d + x, p);
#endif
    }
    if (x + 4 <= width) {
#if KVZ_BIT_DEPTH == 8
      const uint16x4_t z = vdup_n_u16(0);
      const tf_px8_t p = tf_bilinear_pred8_neon(vcombine_u16(tf_load4_neon(r0 + x), z),
                                                vcombine_u16(tf_load4_neon(r0 + x + 1), z),
                                                vcombine_u16(tf_load4_neon(r1 + x), z),
                                                vcombine_u16(tf_load4_neon(r1 + x + 1), z), w, neg_shift);
      vst1_lane_u32((uint32_t *)(void *)(d + x), vreinterpret_u32_u8(p), 0);
#else
      vst1_u16(d + x, tf_bilinear_pred4_neon(tf_load4_neon(r0 + x), tf_load4_neon(r0 + x + 1),
                                             tf_load4_neon(r1 + x), tf_load4_neon(r1 + x + 1), w, neg_shift));
#endif
      x += 4;
    }
    for (; x < width; x++) {
      d[x] = (kvz_pixel)((w00 * r0[x] + w10 * r0[x + 1] + w01 * r1[x] + w11 * r1[x + 1] + round) >> shift);
    }
  }
}

/*
 * 8x8 blocks only. Per 32-bit lane the largest sum, that of the squared
 * gradients, is at most 7 rows * 4 * 8190^2 < 2^32 at 12 bits.
 */
static void tf_block_stats_neon(const kvz_pixel *org, int org_stride,
                                const kvz_pixel *pred, int pred_stride,
                                int width, int height, int64_t out[5])
{
  if (width != 8 || height != 8) {
    kvz_tf_block_stats_generic(org, org_stride, pred, pred_stride, width, height, out);
    return;
  }
  // Lanes 0..6: the pixels that have a right neighbour.
  static const uint16_t mask7_arr[8] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0 };
  const uint16x8_t mask7 = vld1q_u16(mask7_arr);
  uint16x8_t o[8];
  int16x8_t d[8];
  for (int y = 0; y < 8; y++) {
    o[y] = tf_load8_neon(&org[y * org_stride]);
    d[y] = vsubq_s16(vreinterpretq_s16_u16(o[y]), vreinterpretq_s16_u16(tf_load8_neon(&pred[y * pred_stride])));
  }
  uint32x4_t sum = vdupq_n_u32(0), sum2 = vdupq_n_u32(0), sse = vdupq_n_u32(0);
  uint32x4_t var_d = vdupq_n_u32(0), diff = vdupq_n_u32(0);
  for (int y = 0; y < 8; y++) {
    sum = vpadalq_u16(sum, o[y]);
    sum2 = vmlal_u16(sum2, vget_low_u16(o[y]), vget_low_u16(o[y]));
    sum2 = vmlal_high_u16(sum2, o[y], o[y]);
    const uint16x8_t ad = vreinterpretq_u16_s16(vabsq_s16(d[y]));
    sse = vmlal_u16(sse, vget_low_u16(ad), vget_low_u16(ad));
    sse = vmlal_high_u16(sse, ad, ad);
    if (y < 7) {
      const uint16x8_t adm = vandq_u16(ad, mask7);
      var_d = vmlal_u16(var_d, vget_low_u16(adm), vget_low_u16(adm));
      var_d = vmlal_high_u16(var_d, adm, adm);
      const uint16x8_t dh = vandq_u16(vreinterpretq_u16_s16(vabdq_s16(vextq_s16(d[y], d[y], 1), d[y])), mask7);
      const uint16x8_t dv = vandq_u16(vreinterpretq_u16_s16(vabdq_s16(d[y + 1], d[y])), mask7);
      diff = vmlal_u16(diff, vget_low_u16(dh), vget_low_u16(dh));
      diff = vmlal_high_u16(diff, dh, dh);
      diff = vmlal_u16(diff, vget_low_u16(dv), vget_low_u16(dv));
      diff = vmlal_high_u16(diff, dv, dv);
    }
  }
  out[0] = (int64_t)vaddlvq_u32(sum);
  out[1] = (int64_t)vaddlvq_u32(sum2);
  out[2] = (int64_t)vaddlvq_u32(sse);
  out[3] = (int64_t)vaddlvq_u32(var_d);
  out[4] = (int64_t)vaddlvq_u32(diff);
}

/*
 * The operations and their order are those of kvz_tf_filter_block_generic,
 * in single precision without fused multiply-adds, four pixels at a time.
 */
static void tf_filter_block_neon(kvz_pixel *dst, int dst_stride,
                                 const kvz_pixel *src, int src_stride,
                                 const kvz_pixel *const *preds, int pred_stride,
                                 int num_refs, const float *weights,
                                 const float *const *luts,
                                 int width, int height, int max_val)
{
  const int32x4_t maxv = vdupq_n_s32(max_val);
  const int32x4_t zero = vdupq_n_s32(0);
  for (int y = 0; y < height; y++) {
    const kvz_pixel *s = &src[y * src_stride];
    kvz_pixel *d = &dst[y * dst_stride];
    int x = 0;
    for (; x + 4 <= width; x += 4) {
      const int32x4_t o = vreinterpretq_s32_u32(vmovl_u16(tf_load4_neon(s + x)));
      float32x4_t sum = vcvtq_f32_s32(o);
      float32x4_t wsum = vdupq_n_f32(1.0f);
      for (int r = 0; r < num_refs; r++) {
        const int32x4_t p = vreinterpretq_s32_u32(vmovl_u16(tf_load4_neon(&preds[r][y * pred_stride + x])));
        int32_t idx[4];
        vst1q_s32(idx, vsubq_s32(p, o));
        const float *lut = luts[r];
        const float lv[4] = { lut[idx[0]], lut[idx[1]], lut[idx[2]], lut[idx[3]] };
        const float32x4_t wgt = vmulq_f32(vdupq_n_f32(weights[r]), vld1q_f32(lv));
        sum = vaddq_f32(sum, vmulq_f32(wgt, vcvtq_f32_s32(p)));
        wsum = vaddq_f32(wsum, wgt);
      }
      int32x4_t res = vcvtq_s32_f32(vaddq_f32(vdivq_f32(sum, wsum), vdupq_n_f32(0.5f)));
      res = vminq_s32(vmaxq_s32(res, zero), maxv);
      const uint16x4_t res16 = vmovn_u32(vreinterpretq_u32_s32(res));
#if KVZ_BIT_DEPTH == 8
      vst1_lane_u32((uint32_t *)(void *)(d + x), vreinterpret_u32_u8(vmovn_u16(vcombine_u16(res16, res16))), 0);
#else
      vst1_u16(d + x, res16);
#endif
    }
    for (; x < width; x++) {
      const int org_val = s[x];
      float sum = (float)org_val;
      float wsum = 1.0f;
      for (int r = 0; r < num_refs; r++) {
        const int ref_val = preds[r][y * pred_stride + x];
        const float wgt = weights[r] * luts[r][ref_val - org_val];
        sum += wgt * ref_val;
        wsum += wgt;
      }
      d[x] = (kvz_pixel)CLIP(0, max_val, (int)(sum / wsum + 0.5f));
    }
  }
}

} // namespace
#endif // COMPILE_ARM_NEON

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

namespace {

/*
 * 8-bit NEON versions of the functions of picture-generic.c. Every function
 * computes exactly the same result as the generic version, including the
 * block-splitting order of the any-size SATD functions.
 */

//////////////////////////////////////////////////////////////////////////
// SAD

/* Loads a 4-pixel row pair (rows y and y + 1) into one 64-bit vector. */
inline uint8x8_t load_row_pair4(const kvz_pixel *p, unsigned stride)
{
  uint32_t r0, r1;
  memcpy(&r0, p, 4);
  memcpy(&r1, p + stride, 4);
  return vreinterpret_u8_u64(vcreate_u64((uint64_t)r0 | ((uint64_t)r1 << 32)));
}

/*
 * SAD of a W-wide block of `height` rows. |a - b| is accumulated with vabal
 * into 16-bit lanes. A row adds at most ceil(W / 8) * 255 to a lane, so the
 * lanes are flushed into 32-bit lanes every kRows rows, which keeps the sums
 * exact for any height. Only the bytes of the block are read: a 4-wide strip
 * (W = 4 or 12) is loaded as row pairs.
 */
template <int W>
unsigned sad_rows(const kvz_pixel *p1, unsigned s1, const kvz_pixel *p2, unsigned s2, int height)
{
  static_assert(W == 4 || W == 8 || W == 12 || W == 24 || (W % 16 == 0 && W <= 64), "unsupported width");
  constexpr int kPerLane = (W + 7) / 8;
  constexpr int kRows = (65535 / (255 * kPerLane)) & ~1;
  constexpr int kStrip = W % 8 ? W - 4 : -1;          // x of the 4-wide strip
  constexpr int kEight = W % 16 >= 8 ? W - 8 - (kStrip >= 0 ? 4 : 0) : -1; // x of a lone 8-wide part
  constexpr int kSixteen = W - (W % 16);              // columns covered by 16-wide parts

  uint32x4_t total = vdupq_n_u32(0);
  unsigned tail = 0;
  for (int y0 = 0; y0 < height; y0 += kRows) {
    const int y1 = height - y0 < kRows ? height : y0 + kRows;
    uint16x8_t acc = vdupq_n_u16(0);
    if constexpr (kStrip >= 0) {
      int y = y0;
      for (; y + 2 <= y1; y += 2) {
        acc = vabal_u8(acc, load_row_pair4(p1 + y * s1 + kStrip, s1), load_row_pair4(p2 + y * s2 + kStrip, s2));
      }
      for (; y < y1; ++y) {
        for (int x = kStrip; x < W; ++x) {
          tail += abs((int)p1[y * s1 + x] - (int)p2[y * s2 + x]);
        }
      }
    }
    for (int y = y0; y < y1; ++y) {
      const kvz_pixel *r1 = p1 + y * s1;
      const kvz_pixel *r2 = p2 + y * s2;
      if constexpr (kEight >= 0) {
        acc = vabal_u8(acc, vld1_u8(r1 + kEight), vld1_u8(r2 + kEight));
      }
      for (int x = 0; x < kSixteen; x += 16) {
        const uint8x16_t a = vld1q_u8(r1 + x);
        const uint8x16_t b = vld1q_u8(r2 + x);
        acc = vabal_u8(acc, vget_low_u8(a), vget_low_u8(b));
        acc = vabal_high_u8(acc, a, b);
      }
    }
    total = vpadalq_u16(total, acc);
  }
  return vaddvq_u32(total) + tail;
}

/**
 * Sum of absolute differences of a width x height block. Common widths go
 * to the width-specialized kernel; others use the generic row loop, which
 * accumulates in 32-bit lanes.
 */
unsigned reg_sad_neon(const kvz_pixel *const data1, const kvz_pixel *const data2,
                      const int width, const int height,
                      const unsigned stride1, const unsigned stride2)
{
  switch (width) {
    case 4:  return sad_rows<4>(data1, stride1, data2, stride2, height);
    case 8:  return sad_rows<8>(data1, stride1, data2, stride2, height);
    case 12: return sad_rows<12>(data1, stride1, data2, stride2, height);
    case 16: return sad_rows<16>(data1, stride1, data2, stride2, height);
    case 24: return sad_rows<24>(data1, stride1, data2, stride2, height);
    case 32: return sad_rows<32>(data1, stride1, data2, stride2, height);
    case 48: return sad_rows<48>(data1, stride1, data2, stride2, height);
    case 64: return sad_rows<64>(data1, stride1, data2, stride2, height);
    default: break;
  }

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

template <int N>
unsigned sad_nxn_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return sad_rows<N>(block1, N, block2, N, N);
}

template <int N>
void sad_nxn_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = sad_rows<N>(preds[0], N, orig, N, N);
  costs_out[1] = sad_rows<N>(preds[1], N, orig, N, N);
}

/* Fixed-width SAD entry points for get_optimized_sad. */
template <int W>
uint32_t reg_sad_w_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                        const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return sad_rows<W>(data1, stride1, data2, stride2, height);
}

uint32_t reg_sad_w0_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                         const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return 0;
}

optimized_sad_func_ptr_t get_optimized_sad_neon(int32_t width)
{
  switch (width) {
    case 0:  return reg_sad_w0_neon;
    case 4:  return reg_sad_w_neon<4>;
    case 8:  return reg_sad_w_neon<8>;
    case 12: return reg_sad_w_neon<12>;
    case 16: return reg_sad_w_neon<16>;
    case 24: return reg_sad_w_neon<24>;
    case 32: return reg_sad_w_neon<32>;
    case 64: return reg_sad_w_neon<64>;
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
  if (block_width < 8) {
    /* Narrow frame-border strips: a plain loop beats per-row vector setup. */
    unsigned sum = 0;
    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *p = pic_data + y * pic_stride;
      const int ref = ref_data[y * ref_stride];
      for (int x = 0; x < block_width; ++x) {
        sum += abs((int)p[x] - ref);
      }
    }
    return sum;
  }

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

/*
 * The last butterfly stage of the second Hadamard pass feeds straight into the
 * absolute sum, and |a + b| + |a - b| == 2 * max(|a|, |b|) for integers. So the
 * last stage and half of the absolute values are replaced by maxima. The sums
 * are exact: the 4x4 second-pass terms are at most 2040 and the 8x8 terms at
 * most 8160 (four of them, 32640, still fit a 16-bit lane).
 */
static INLINE uint32_t hadamard4_pass_absmax_sum(const int16x4_t r[4])
{
  const int16x4_t s03p = vadd_s16(r[0], r[3]);
  const int16x4_t s03m = vsub_s16(r[0], r[3]);
  const int16x4_t s12p = vadd_s16(r[1], r[2]);
  const int16x4_t s12m = vsub_s16(r[1], r[2]);
  const int16x4_t mp = vmax_s16(vabs_s16(s03p), vabs_s16(s12p));
  const int16x4_t mm = vmax_s16(vabs_s16(s03m), vabs_s16(s12m));
  return (uint32_t)vaddlv_s16(vadd_s16(mp, mm)); // half of the 4x4 absolute sum
}

static INLINE uint32_t hadamard8_pass_absmax_sum(const int16x8_t m[8])
{
  const int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  const int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  const int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  const int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  const int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  const int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  const int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  const int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  const int16x8_t m01 = vmaxq_s16(vabsq_s16(b0), vabsq_s16(b1));
  const int16x8_t m23 = vmaxq_s16(vabsq_s16(b2), vabsq_s16(b3));
  const int16x8_t m45 = vmaxq_s16(vabsq_s16(b4), vabsq_s16(b5));
  const int16x8_t m67 = vmaxq_s16(vabsq_s16(b6), vabsq_s16(b7));
  const uint16x8_t sum = vreinterpretq_u16_s16(vaddq_s16(vaddq_s16(m01, m23), vaddq_s16(m45, m67)));
  return vaddlvq_u16(sum); // half of the 8x8 absolute sum
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
  // Generic returns (sum + 1) >> 1 of the absolute sum, which is 2 * half.
  return hadamard4_pass_absmax_sum(r);
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
  // Generic returns (sum + 2) >> 2 of the absolute sum, with sum = 2 * half.
  return (2 * hadamard8_pass_absmax_sum(m) + 2) >> 2;
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

template <int N>
unsigned satd_nxn_neon(const kvz_pixel * const block1, const kvz_pixel * const block2)
{
  unsigned sum = 0;
  for (int y = 0; y < N; y += 8) {
    for (int x = 0; x < N; x += 8) {
      sum += satd_8x8_subblock_neon(&block1[y * N + x], N, &block2[y * N + x], N);
    }
  }
  return sum;
}

template <int N>
void satd_nxn_dual_neon(const pred_buffer preds, const kvz_pixel * const orig,
                        unsigned num_modes, unsigned *costs_out)
{
  for (int k = 0; k < 2; ++k) {
    costs_out[k] = satd_nxn_neon<N>(preds[k], orig);
  }
}

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
  /* |a - b| fits in 8 bits and its square in 16, so square the absolute
   * difference with a widening multiply and accumulate pairwise into 32 bits. */
  uint32x4_t acc = vdupq_n_u32(0);
  if (width == 4) {
    for (int y = 0; y < 4; y += 2) {
      const uint8x8_t d = vabd_u8(load_row_pair4(ref + y * ref_stride, ref_stride),
                                  load_row_pair4(rec + y * rec_stride, rec_stride));
      acc = vpadalq_u16(acc, vmull_u8(d, d));
    }
    return vaddvq_u32(acc);
  }
  unsigned tail = 0;
  for (int y = 0; y < width; ++y) {
    const kvz_pixel *r = ref + y * ref_stride;
    const kvz_pixel *c = rec + y * rec_stride;
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      const uint8x16_t d = vabdq_u8(vld1q_u8(r + x), vld1q_u8(c + x));
      acc = vpadalq_u16(acc, vmull_u8(vget_low_u8(d), vget_low_u8(d)));
      acc = vpadalq_u16(acc, vmull_high_u8(d, d));
    }
    if (x + 8 <= width) {
      const uint8x8_t d = vabd_u8(vld1_u8(r + x), vld1_u8(c + x));
      acc = vpadalq_u16(acc, vmull_u8(d, d));
      x += 8;
    }
    for (; x < width; ++x) {
      const int diff = (int)r[x] - (int)c[x];
      tail += (unsigned)(diff * diff);
    }
  }
  return vaddvq_u32(acc) + tail;
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
  /* The encoder calls this on whole frames (VAQ), so the sum of squares needs
   * 64 bits. Each 16-pixel step adds at most 4 * 255^2 to a 32-bit lane, so
   * the lanes are flushed into 64-bit accumulators every 4096 steps. */
  uint64x2_t sum64 = vdupq_n_u64(0);
  uint64x2_t sq64 = vdupq_n_u64(0);
  uint32_t i = 0;

  while (i + 16 <= len) {
    uint32x4_t sum_acc = vdupq_n_u32(0);
    uint32x4_t sq_acc = vdupq_n_u32(0);
    for (int step = 0; step < 4096 && i + 16 <= len; ++step, i += 16) {
      uint8x16_t v = vld1q_u8(arr + i);
      sum_acc = vpadalq_u16(sum_acc, vpaddlq_u8(v));
      sq_acc = vpadalq_u16(sq_acc, vmull_u8(vget_low_u8(v), vget_low_u8(v)));
      sq_acc = vpadalq_u16(sq_acc, vmull_high_u8(v, v));
    }
    sum64 = vpadalq_u32(sum64, sum_acc);
    sq64 = vpadalq_u32(sq64, sq_acc);
  }

  uint64_t sum = vaddvq_u64(sum64);
  uint64_t sumsq = vaddvq_u64(sq64);
  for (; i < len; ++i) {
    sum += arr[i];
    sumsq += (uint32_t)arr[i] * arr[i];
  }

  /* Population variance (n * sum(x^2) - sum(x)^2) / n^2, with the numerator
   * computed exactly. Like the AVX2 version this is not bit-identical to the
   * generic two-pass double loop, only more accurate. */
  const unsigned __int128 num = (unsigned __int128)len * sumsq - (unsigned __int128)sum * sum;
  const double n = (double)len;
  return (double)num / (n * n);
}

} // namespace

#endif // COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

int kvz_strategy_register_picture_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_ARM_NEON
  if (bitdepth == KVZ_BIT_DEPTH) {
#if KVZ_BIT_DEPTH > 8
    success &= kvz_strategyselector_register(opaque, "satd_4x4", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_neon<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_dual_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_neon<64>));
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_quad_hbd_neon));
#endif // KVZ_BIT_DEPTH > 8

    success &= kvz_strategyselector_register(opaque, "bilinear_ssd", "neon", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_neon));
    success &= kvz_strategyselector_register(opaque, "bilinear_ssd_multi", "neon", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_multi_neon));
    success &= kvz_strategyselector_register(opaque, "tf_bilinear_pred", "neon", 40, KVZ_STRATEGY_FPTR(tf_bilinear_pred_neon));
    success &= kvz_strategyselector_register(opaque, "tf_block_stats", "neon", 40, KVZ_STRATEGY_FPTR(tf_block_stats_neon));
    success &= kvz_strategyselector_register(opaque, "tf_filter_block", "neon", 40, KVZ_STRATEGY_FPTR(tf_filter_block_neon));
  }
#endif // COMPILE_ARM_NEON
#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "reg_sad", "neon", 40, KVZ_STRATEGY_FPTR(reg_sad_neon));

    success &= kvz_strategyselector_register(opaque, "sad_4x4", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<64>));

    success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_dual_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<64>));
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_neon));
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_quad_neon));

    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "neon", 40, KVZ_STRATEGY_FPTR(pixels_calc_ssd_neon));
    success &= kvz_strategyselector_register(opaque, "bipred_average", "neon", 40, KVZ_STRATEGY_FPTR(bipred_average_neon));

    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "neon", 40, KVZ_STRATEGY_FPTR(get_optimized_sad_neon));
    success &= kvz_strategyselector_register(opaque, "ver_sad", "neon", 40, KVZ_STRATEGY_FPTR(ver_sad_neon));
    success &= kvz_strategyselector_register(opaque, "hor_sad", "neon", 40, KVZ_STRATEGY_FPTR(hor_sad_neon_wrapper));

    success &= kvz_strategyselector_register(opaque, "pixel_var", "neon", 40, KVZ_STRATEGY_FPTR(pixel_var_neon));
  }
#endif
  return success;
}
