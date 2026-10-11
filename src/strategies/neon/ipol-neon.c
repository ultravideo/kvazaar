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
 * AArch64 NEON implementation of the subpixel interpolation sample functions.
 *
 * Bit-exact with the generic implementation for 8-bit video. For 8-bit input
 * the horizontal filter result always fits in int16 (|sum| < 2^15), so it is
 * computed with 16-bit multiply-accumulate. The vertical filter needs 32-bit
 * accumulators because it operates on those 16-bit intermediates.
 */

#include <string.h>

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "search_inter.h"
#include "strategies/neon/ipol-neon.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-ipol.h"
#include "strategyselector.h"

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

extern const int8_t kvz_g_luma_filter[4][8];
extern const int8_t kvz_g_chroma_filter[8][4];

static inline kvz_pixel clip_to_pixel_8bit(int32_t v)
{
  return (kvz_pixel)(v < 0 ? 0 : (v > 255 ? 255 : v));
}

/*
 * Horizontal pass. `base` points to the sample at (-pad, -pad) relative to the
 * block origin, so output (y, x) uses the taps base[y * stride + x + k].
 * Writes `rows` rows of `width` intermediate values into dst.
 */
/* 32-bit vertical accumulation for 8 consecutive columns starting at p = &hor[y][x]. */
static inline void ver_sum8_s32(const int16_t *p, int16_t stride, const int8_t *f, int taps,
                                int32x4_t *lo, int32x4_t *hi)
{
  int32x4_t acc_lo = vdupq_n_s32(0);
  int32x4_t acc_hi = vdupq_n_s32(0);
  for (int k = 0; k < taps; ++k) {
    int16x8_t v = vld1q_s16(p + (intptr_t)k * stride);
    acc_lo = vmlal_n_s16(acc_lo, vget_low_s16(v), f[k]);
    acc_hi = vmlal_n_s16(acc_hi, vget_high_s16(v), f[k]);
  }
  *lo = acc_lo;
  *hi = acc_hi;
}

/* Full-sample output for 8-bit: clip(((sum >> 6) + 32) >> 6). */
static inline uint8x8_t ver_finish_pixel(int32x4_t lo, int32x4_t hi)
{
  const int32x4_t off = vdupq_n_s32(32); /* wp_offset1 = 1 << (wp_shift1 - 1), wp_shift1 = 6 */
  lo = vshrq_n_s32(vaddq_s32(vshrq_n_s32(lo, 6), off), 6);
  hi = vshrq_n_s32(vaddq_s32(vshrq_n_s32(hi, 6), off), 6);
  return vqmovun_s16(vcombine_s16(vqmovn_s32(lo), vqmovn_s32(hi)));
}

/* Like ver_finish_pixel, but stores sum >> 6 in int16 first, as the generic
 * block filters do (int16_t sample = ... >> shift2). That store wraps when the
 * half-pel diagonal sum exceeds int16, so the vector path must wrap too.
 * vrshrq_n_s16 computes (x + 32) >> 6 without an intermediate overflow. */
static inline uint8x8_t ver_finish_pixel_wrap16(int32x4_t lo, int32x4_t hi)
{
  int16x8_t s = vcombine_s16(vmovn_s32(vshrq_n_s32(lo, 6)), vmovn_s32(vshrq_n_s32(hi, 6)));
  return vqmovun_s16(vrshrq_n_s16(s, 6));
}

#if defined(_MSC_VER)
#define IPOL_INLINE __forceinline
#else
#define IPOL_INLINE inline __attribute__((always_inline))
#endif

static inline void hshift_clip_row(const int16_t *h, int n, kvz_pixel *out);

/*
 * Motion-compensation sample passes. The filters are always rows of
 * kvz_g_luma_filter or kvz_g_chroma_filter, so each pass dispatches on the row
 * and runs with constant taps (zero taps skipped, u8 widening MACs for the
 * horizontal pass). Row 0 of either table is a single tap of 64, which makes
 * the horizontal pass a shift and the vertical pass a rounding shift or copy;
 * this is the integer motion vector case. All variants are exact.
 */
static IPOL_INLINE uint16x8_t tapc_u8(uint16x8_t acc, const kvz_pixel *p, int c)
{
  if (c > 0) return vmlal_u8(acc, vld1_u8(p), vdup_n_u8((uint8_t)c));
  if (c < 0) return vmlsl_u8(acc, vld1_u8(p), vdup_n_u8((uint8_t)-c));
  return acc;
}

static IPOL_INLINE void tapc_s16(int32x4_t *lo, int32x4_t *hi, const int16_t *p, int c)
{
  if (c == 0) return;
  const int16x8_t v = vld1q_s16(p);
  *lo = vmlal_n_s16(*lo, vget_low_s16(v), (int16_t)c);
  *hi = vmlal_high_n_s16(*hi, v, (int16_t)c);
}

static IPOL_INLINE void hor_rows_c(const kvz_pixel *base, int16_t src_stride, int width, int rows,
                                   int16_t dst[][LCU_WIDTH],
                                   int c0, int c1, int c2, int c3, int c4, int c5, int c6, int c7)
{
  for (int y = 0; y < rows; ++y) {
    const kvz_pixel *s = base + (intptr_t)y * src_stride;
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      uint16x8_t acc = vdupq_n_u16(0);
      acc = tapc_u8(acc, s + x + 0, c0); acc = tapc_u8(acc, s + x + 1, c1);
      acc = tapc_u8(acc, s + x + 2, c2); acc = tapc_u8(acc, s + x + 3, c3);
      acc = tapc_u8(acc, s + x + 4, c4); acc = tapc_u8(acc, s + x + 5, c5);
      acc = tapc_u8(acc, s + x + 6, c6); acc = tapc_u8(acc, s + x + 7, c7);
      vst1q_s16(&dst[y][x], vreinterpretq_s16_u16(acc));
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      if (c0) sum += c0 * s[x + 0];
      if (c1) sum += c1 * s[x + 1];
      if (c2) sum += c2 * s[x + 2];
      if (c3) sum += c3 * s[x + 3];
      if (c4) sum += c4 * s[x + 4];
      if (c5) sum += c5 * s[x + 5];
      if (c6) sum += c6 * s[x + 6];
      if (c7) sum += c7 * s[x + 7];
      dst[y][x] = (int16_t)sum;
    }
  }
}

// A filter row that is a single tap of 64 at position t: dst = sample << 6.
static IPOL_INLINE void hor_rows_shift(const kvz_pixel *base, int16_t src_stride, int width, int rows,
                                       int16_t dst[][LCU_WIDTH], int t)
{
  for (int y = 0; y < rows; ++y) {
    const kvz_pixel *s = base + (intptr_t)y * src_stride + t;
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      vst1q_s16(&dst[y][x], vreinterpretq_s16_u16(vshll_n_u8(vld1_u8(s + x), 6)));
    }
    for (; x < width; ++x) dst[y][x] = (int16_t)(s[x] << 6);
  }
}

static void hor_pass_any(const kvz_pixel *base, int16_t src_stride,
                         int width, int rows, const int8_t *f, int taps,
                         int16_t dst[][LCU_WIDTH])
{
  for (int y = 0; y < rows; ++y) {
    const kvz_pixel *s = base + (intptr_t)y * src_stride;
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int16x8_t acc = vdupq_n_s16(0);
      for (int k = 0; k < taps; ++k) {
        int16x8_t v = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s + x + k)));
        acc = vmlaq_n_s16(acc, v, f[k]);
      }
      vst1q_s16(&dst[y][x], acc);
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += f[k] * s[x + k];
      }
      dst[y][x] = (int16_t)sum;
    }
  }
}

/*
 * Horizontal pass. `base` points to the sample at (-pad, -pad) relative to the
 * block origin, so output (y, x) uses the taps base[y * stride + x + k].
 * Writes `rows` rows of `width` intermediate values into dst.
 */
static void hor_pass_s16(const kvz_pixel *base, int16_t src_stride,
                         int width, int rows, const int8_t *f, int taps,
                         int16_t dst[][LCU_WIDTH])
{
#define HOR(...) hor_rows_c(base, src_stride, width, rows, dst, __VA_ARGS__)
  if (taps == 8) {
    if (f == kvz_g_luma_filter[0]) { hor_rows_shift(base, src_stride, width, rows, dst, 3); return; }
    if (f == kvz_g_luma_filter[1]) { HOR(-1, 4, -10, 58, 17, -5, 1, 0); return; }
    if (f == kvz_g_luma_filter[2]) { HOR(-1, 4, -11, 40, 40, -11, 4, -1); return; }
    if (f == kvz_g_luma_filter[3]) { HOR(0, 1, -5, 17, 58, -10, 4, -1); return; }
  } else if (taps == 4) {
    if (f == kvz_g_chroma_filter[0]) { hor_rows_shift(base, src_stride, width, rows, dst, 1); return; }
    if (f == kvz_g_chroma_filter[1]) { HOR(-2, 58, 10, -2, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[2]) { HOR(-4, 54, 16, -2, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[3]) { HOR(-6, 46, 28, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[4]) { HOR(-4, 36, 36, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[5]) { HOR(-4, 28, 46, -6, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[6]) { HOR(-2, 16, 54, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[7]) { HOR(-2, 10, 58, -2, 0, 0, 0, 0); return; }
  }
#undef HOR
  hor_pass_any(base, src_stride, width, rows, f, taps, dst);
}

// Vertical pass with constant taps; hi_out selects the int16 (hi-precision) output.
static IPOL_INLINE void ver_rows_c(const int16_t hor[][LCU_WIDTH], int width, int height,
                                   kvz_pixel *dst, int16_t *dst_hi, int16_t dst_stride, int hi_out,
                                   int c0, int c1, int c2, int c3, int c4, int c5, int c6, int c7)
{
  for (int y = 0; y < height; ++y) {
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int32x4_t lo = vdupq_n_s32(0), hi = vdupq_n_s32(0);
      tapc_s16(&lo, &hi, &hor[y + 0][x], c0); tapc_s16(&lo, &hi, &hor[y + 1][x], c1);
      tapc_s16(&lo, &hi, &hor[y + 2][x], c2); tapc_s16(&lo, &hi, &hor[y + 3][x], c3);
      tapc_s16(&lo, &hi, &hor[y + 4][x], c4); tapc_s16(&lo, &hi, &hor[y + 5][x], c5);
      tapc_s16(&lo, &hi, &hor[y + 6][x], c6); tapc_s16(&lo, &hi, &hor[y + 7][x], c7);
      if (hi_out) {
        // (sum >> 6) stored as int16, same as the generic truncating store.
        vst1q_s16(&dst_hi[y * dst_stride + x], vcombine_s16(vshrn_n_s32(lo, 6), vshrn_n_s32(hi, 6)));
      } else {
        vst1_u8(&dst[y * dst_stride + x], ver_finish_pixel(lo, hi));
      }
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      if (c0) sum += c0 * hor[y + 0][x];
      if (c1) sum += c1 * hor[y + 1][x];
      if (c2) sum += c2 * hor[y + 2][x];
      if (c3) sum += c3 * hor[y + 3][x];
      if (c4) sum += c4 * hor[y + 4][x];
      if (c5) sum += c5 * hor[y + 5][x];
      if (c6) sum += c6 * hor[y + 6][x];
      if (c7) sum += c7 * hor[y + 7][x];
      if (hi_out) dst_hi[y * dst_stride + x] = (int16_t)(sum >> 6);
      else dst[y * dst_stride + x] = clip_to_pixel_8bit(((sum >> 6) + 32) >> 6);
    }
  }
}

static void ver_pass_any(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                         int width, int height, kvz_pixel *dst, int16_t *dst_hi,
                         int16_t dst_stride, int hi_out)
{
  for (int y = 0; y < height; ++y) {
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int32x4_t lo, hi;
      ver_sum8_s32(&hor[y][x], LCU_WIDTH, f, taps, &lo, &hi);
      if (hi_out) {
        vst1q_s16(&dst_hi[y * dst_stride + x], vcombine_s16(vshrn_n_s32(lo, 6), vshrn_n_s32(hi, 6)));
      } else {
        vst1_u8(&dst[y * dst_stride + x], ver_finish_pixel(lo, hi));
      }
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += f[k] * hor[y + k][x];
      }
      if (hi_out) dst_hi[y * dst_stride + x] = (int16_t)(sum >> 6);
      else dst[y * dst_stride + x] = clip_to_pixel_8bit(((sum >> 6) + 32) >> 6);
    }
  }
}

static IPOL_INLINE void ver_pass(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                                 int width, int height, kvz_pixel *dst, int16_t *dst_hi,
                                 int16_t dst_stride, int hi_out)
{
#define VER(...) ver_rows_c(hor, width, height, dst, dst_hi, dst_stride, hi_out, __VA_ARGS__)
  // Row 0: sum = 64 * h[t], so sum >> 6 == h[t]: a copy (hi) or a rounding shift (pixel).
  int single = -1;
  if (taps == 8 && f == kvz_g_luma_filter[0]) single = 3;
  if (taps == 4 && f == kvz_g_chroma_filter[0]) single = 1;
  if (single >= 0) {
    for (int y = 0; y < height; ++y) {
      if (hi_out) memcpy(&dst_hi[y * dst_stride], hor[y + single], (size_t)width * sizeof(int16_t));
      else hshift_clip_row(hor[y + single], width, &dst[y * dst_stride]);
    }
    return;
  }
  if (taps == 8) {
    if (f == kvz_g_luma_filter[1]) { VER(-1, 4, -10, 58, 17, -5, 1, 0); return; }
    if (f == kvz_g_luma_filter[2]) { VER(-1, 4, -11, 40, 40, -11, 4, -1); return; }
    if (f == kvz_g_luma_filter[3]) { VER(0, 1, -5, 17, 58, -10, 4, -1); return; }
  } else if (taps == 4) {
    if (f == kvz_g_chroma_filter[1]) { VER(-2, 58, 10, -2, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[2]) { VER(-4, 54, 16, -2, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[3]) { VER(-6, 46, 28, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[4]) { VER(-4, 36, 36, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[5]) { VER(-4, 28, 46, -6, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[6]) { VER(-2, 16, 54, -4, 0, 0, 0, 0); return; }
    if (f == kvz_g_chroma_filter[7]) { VER(-2, 10, 58, -2, 0, 0, 0, 0); return; }
  }
#undef VER
  ver_pass_any(hor, f, taps, width, height, dst, dst_hi, dst_stride, hi_out);
}

static void ver_pass_pixel(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                           int width, int height,
                           kvz_pixel *dst, int16_t dst_stride)
{
  ver_pass(hor, f, taps, width, height, dst, NULL, dst_stride, 0);
}

static void ver_pass_hi(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                        int width, int height,
                        int16_t *dst, int16_t dst_stride)
{
  ver_pass(hor, f, taps, width, height, NULL, dst, dst_stride, 1);
}

static void sample_quarterpel_luma_neon(const encoder_control_t * const encoder,
                                        kvz_pixel *src, int16_t src_stride,
                                        int width, int height,
                                        kvz_pixel *dst, int16_t dst_stride,
                                        int8_t hor_flag, int8_t ver_flag,
                                        const int16_t mv[2])
{
  (void)encoder; (void)hor_flag; (void)ver_flag;
  int16_t hor[KVZ_EXT_BLOCK_W_LUMA][LCU_WIDTH];
  const kvz_pixel *base = src - KVZ_LUMA_FILTER_OFFSET * src_stride - KVZ_LUMA_FILTER_OFFSET;

  hor_pass_s16(base, src_stride, width, height + KVZ_EXT_PADDING_LUMA,
               kvz_g_luma_filter[mv[0] & 3], 8, hor);
  ver_pass_pixel(hor, kvz_g_luma_filter[mv[1] & 3], 8, width, height, dst, dst_stride);
}

static void sample_quarterpel_luma_hi_neon(const encoder_control_t * const encoder,
                                           kvz_pixel *src, int16_t src_stride,
                                           int width, int height,
                                           int16_t *dst, int16_t dst_stride,
                                           int8_t hor_flag, int8_t ver_flag,
                                           const int16_t mv[2])
{
  (void)encoder; (void)hor_flag; (void)ver_flag;
  int16_t hor[KVZ_EXT_BLOCK_W_LUMA][LCU_WIDTH];
  const kvz_pixel *base = src - KVZ_LUMA_FILTER_OFFSET * src_stride - KVZ_LUMA_FILTER_OFFSET;

  hor_pass_s16(base, src_stride, width, height + KVZ_EXT_PADDING_LUMA,
               kvz_g_luma_filter[mv[0] & 3], 8, hor);
  ver_pass_hi(hor, kvz_g_luma_filter[mv[1] & 3], 8, width, height, dst, dst_stride);
}

static void sample_octpel_chroma_neon(const encoder_control_t * const encoder,
                                      kvz_pixel *src, int16_t src_stride,
                                      int width, int height,
                                      kvz_pixel *dst, int16_t dst_stride,
                                      int8_t hor_flag, int8_t ver_flag,
                                      const int16_t mv[2])
{
  (void)hor_flag; (void)ver_flag;
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  const int8_t *hor_filter = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  const int8_t *ver_filter = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  int16_t hor[KVZ_EXT_BLOCK_W_LUMA][LCU_WIDTH];
  const kvz_pixel *base = src - KVZ_CHROMA_FILTER_OFFSET * src_stride - KVZ_CHROMA_FILTER_OFFSET;

  hor_pass_s16(base, src_stride, width, height + KVZ_EXT_PADDING_CHROMA, hor_filter, 4, hor);
  ver_pass_pixel(hor, ver_filter, 4, width, height, dst, dst_stride);
}

static void sample_octpel_chroma_hi_neon(const encoder_control_t * const encoder,
                                         kvz_pixel *src, int16_t src_stride,
                                         int width, int height,
                                         int16_t *dst, int16_t dst_stride,
                                         int8_t hor_flag, int8_t ver_flag,
                                         const int16_t mv[2])
{
  (void)hor_flag; (void)ver_flag;
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  const int8_t *hor_filter = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  const int8_t *ver_filter = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  int16_t hor[KVZ_EXT_BLOCK_W_LUMA][LCU_WIDTH];
  const kvz_pixel *base = src - KVZ_CHROMA_FILTER_OFFSET * src_stride - KVZ_CHROMA_FILTER_OFFSET;

  hor_pass_s16(base, src_stride, width, height + KVZ_EXT_PADDING_CHROMA, hor_filter, 4, hor);
  ver_pass_hi(hor, ver_filter, 4, width, height, dst, dst_stride);
}

/* ---- Block filters (hpel/qpel, hor_ver and diag) ----
 * The control flow below mirrors ipol-generic.c. Only the per-row inner loops
 * are vectorized. Single-column samples that the generic code computes
 * outside the vector loops are kept scalar. */

// Scalar 8-tap horizontal sample from 8-bit samples (shift1 is 0 for 8-bit).
static inline int16_t hor8_px_scalar(const kvz_pixel *p, const int8_t *f)
{
  int32_t s = 0;
  for (int k = 0; k < 8; ++k) {
    s += f[k] * p[k];
  }
  return (int16_t)s;
}

// Scalar 8-tap vertical sample over int16 intermediates, then the two-stage
// rounding and clip used by the qpel column samples.
static inline kvz_pixel qpel_col_sample(const int16_t *col, const int8_t *f)
{
  int32_t acc = 0;
  for (int k = 0; k < 8; ++k) {
    acc += f[k] * col[k];
  }
  int16_t sample = (int16_t)(acc >> 6);
  return clip_to_pixel_8bit(((int32_t)sample + 32) >> 6);
}

// out[i] = sum_k f[k] * p[k + i] for i < n, exact in int16 for 8-bit input.
static inline void hor8_row_any(const kvz_pixel *p, const int8_t *f, int n, int16_t *out)
{
  int i = 0;
  for (; i + 8 <= n; i += 8) {
    int16x8_t acc = vdupq_n_s16(0);
    for (int k = 0; k < 8; ++k) {
      int16x8_t v = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(p + i + k)));
      acc = vmlaq_n_s16(acc, v, f[k]);
    }
    vst1q_s16(out + i, acc);
  }
  for (; i < n; ++i) {
    int32_t s = 0;
    for (int k = 0; k < 8; ++k) {
      s += f[k] * p[i + k];
    }
    out[i] = (int16_t)s;
  }
}

/*
 * Constant-tap versions of the three row helpers. The block filters only pass
 * rows of kvz_g_luma_filter, so dispatching on the row makes every tap a
 * compile-time constant: zero taps disappear, and the 8-bit passes use
 * widening u8 multiply-accumulates (vmlal_u8/vmlsl_u8) without separate
 * widening moves. The u16 accumulation is modular and the true sum fits int16
 * for 8-bit input, so reinterpreting it as int16 is exact.
 */
static inline uint16x8_t tap_u8(uint16x8_t acc, const kvz_pixel *p, int c)
{
  if (c > 0) return vmlal_u8(acc, vld1_u8(p), vdup_n_u8((uint8_t)c));
  if (c < 0) return vmlsl_u8(acc, vld1_u8(p), vdup_n_u8((uint8_t)-c));
  return acc;
}

static inline void tap_s16(int32x4_t *lo, int32x4_t *hi, const int16_t *p, int c)
{
  if (c == 0) return;
  const int16x8_t v = vld1q_s16(p);
  *lo = vmlal_n_s16(*lo, vget_low_s16(v), (int16_t)c);
  *hi = vmlal_high_n_s16(*hi, v, (int16_t)c);
}

#define DEF_ROW_HELPERS(ID, C0, C1, C2, C3, C4, C5, C6, C7) \
static inline void hor8_row_##ID(const kvz_pixel *p, int n, int16_t *out) \
{ \
  static const int8_t c[8] = { C0, C1, C2, C3, C4, C5, C6, C7 }; \
  int i = 0; \
  for (; i + 8 <= n; i += 8) { \
    uint16x8_t acc = vdupq_n_u16(0); \
    acc = tap_u8(acc, p + i + 0, C0); acc = tap_u8(acc, p + i + 1, C1); \
    acc = tap_u8(acc, p + i + 2, C2); acc = tap_u8(acc, p + i + 3, C3); \
    acc = tap_u8(acc, p + i + 4, C4); acc = tap_u8(acc, p + i + 5, C5); \
    acc = tap_u8(acc, p + i + 6, C6); acc = tap_u8(acc, p + i + 7, C7); \
    vst1q_s16(out + i, vreinterpretq_s16_u16(acc)); \
  } \
  for (; i < n; ++i) { \
    int32_t sum = 0; \
    for (int k = 0; k < 8; ++k) sum += c[k] * p[i + k]; \
    out[i] = (int16_t)sum; \
  } \
} \
static inline void ver8_px_row_##ID(const kvz_pixel *p, int16_t stride, int n, kvz_pixel *out) \
{ \
  static const int8_t c[8] = { C0, C1, C2, C3, C4, C5, C6, C7 }; \
  int i = 0; \
  for (; i + 8 <= n; i += 8) { \
    uint16x8_t acc = vdupq_n_u16(0); \
    acc = tap_u8(acc, p + 0 * (intptr_t)stride + i, C0); acc = tap_u8(acc, p + 1 * (intptr_t)stride + i, C1); \
    acc = tap_u8(acc, p + 2 * (intptr_t)stride + i, C2); acc = tap_u8(acc, p + 3 * (intptr_t)stride + i, C3); \
    acc = tap_u8(acc, p + 4 * (intptr_t)stride + i, C4); acc = tap_u8(acc, p + 5 * (intptr_t)stride + i, C5); \
    acc = tap_u8(acc, p + 6 * (intptr_t)stride + i, C6); acc = tap_u8(acc, p + 7 * (intptr_t)stride + i, C7); \
    /* (sum + 32) >> 6 with clipping; sum fits int16, so vqrshrun is exact. */ \
    vst1_u8(out + i, vqrshrun_n_s16(vreinterpretq_s16_u16(acc), 6)); \
  } \
  for (; i < n; ++i) { \
    int32_t sum = 0; \
    for (int k = 0; k < 8; ++k) sum += c[k] * p[(intptr_t)k * stride + i]; \
    out[i] = clip_to_pixel_8bit((sum + 32) >> 6); \
  } \
} \
static inline void ver8_s16_row_##ID(const int16_t *p, int16_t stride, int n, kvz_pixel *out) \
{ \
  static const int8_t c[8] = { C0, C1, C2, C3, C4, C5, C6, C7 }; \
  int i = 0; \
  for (; i + 8 <= n; i += 8) { \
    int32x4_t lo = vdupq_n_s32(0), hi = vdupq_n_s32(0); \
    tap_s16(&lo, &hi, p + 0 * (intptr_t)stride + i, C0); tap_s16(&lo, &hi, p + 1 * (intptr_t)stride + i, C1); \
    tap_s16(&lo, &hi, p + 2 * (intptr_t)stride + i, C2); tap_s16(&lo, &hi, p + 3 * (intptr_t)stride + i, C3); \
    tap_s16(&lo, &hi, p + 4 * (intptr_t)stride + i, C4); tap_s16(&lo, &hi, p + 5 * (intptr_t)stride + i, C5); \
    tap_s16(&lo, &hi, p + 6 * (intptr_t)stride + i, C6); tap_s16(&lo, &hi, p + 7 * (intptr_t)stride + i, C7); \
    vst1_u8(out + i, ver_finish_pixel_wrap16(lo, hi)); \
  } \
  for (; i < n; ++i) { \
    int32_t sum = 0; \
    for (int k = 0; k < 8; ++k) sum += c[k] * p[(intptr_t)k * stride + i]; \
    int16_t sample = (int16_t)(sum >> 6); \
    out[i] = clip_to_pixel_8bit(((int32_t)sample + 32) >> 6); \
  } \
}

DEF_ROW_HELPERS(1, -1, 4, -10, 58, 17, -5, 1, 0)
DEF_ROW_HELPERS(2, -1, 4, -11, 40, 40, -11, 4, -1)
DEF_ROW_HELPERS(3, 0, 1, -5, 17, 58, -10, 4, -1)

/* Filter row 0 is { 0, 0, 0, 64, 0, 0, 0, 0 }: the horizontal result is the
 * pixel times 64, the vertical pixel result is the pixel itself, and the
 * vertical result over intermediates is ((64 * h) >> 6) == h, an int16 value,
 * so the output is clip((h + 32) >> 6). */
static inline void hor8_row(const kvz_pixel *p, const int8_t *f, int n, int16_t *out)
{
  if (f == kvz_g_luma_filter[0]) {
    int i = 0;
    for (; i + 8 <= n; i += 8) {
      vst1q_s16(out + i, vreinterpretq_s16_u16(vshll_n_u8(vld1_u8(p + i + 3), 6)));
    }
    for (; i < n; ++i) out[i] = (int16_t)(p[i + 3] << 6);
  }
  else if (f == kvz_g_luma_filter[1]) hor8_row_1(p, n, out);
  else if (f == kvz_g_luma_filter[2]) hor8_row_2(p, n, out);
  else if (f == kvz_g_luma_filter[3]) hor8_row_3(p, n, out);
  else hor8_row_any(p, f, n, out);
}

static inline void hshift_clip_row(const int16_t *h, int n, kvz_pixel *out);
static inline void ver8_px_row_any(const kvz_pixel *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out);
static inline void ver8_s16_row_any(const int16_t *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out);

static inline void ver8_px_row(const kvz_pixel *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
{
  if (f == kvz_g_luma_filter[1]) ver8_px_row_1(p, stride, n, out);
  else if (f == kvz_g_luma_filter[2]) ver8_px_row_2(p, stride, n, out);
  else if (f == kvz_g_luma_filter[3]) ver8_px_row_3(p, stride, n, out);
  else ver8_px_row_any(p, stride, f, n, out);
}


static inline void ver8_s16_row(const int16_t *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
{
  if (f == kvz_g_luma_filter[0]) hshift_clip_row(p + 3 * (intptr_t)stride, n, out);
  else if (f == kvz_g_luma_filter[1]) ver8_s16_row_1(p, stride, n, out);
  else if (f == kvz_g_luma_filter[2]) ver8_s16_row_2(p, stride, n, out);
  else if (f == kvz_g_luma_filter[3]) ver8_s16_row_3(p, stride, n, out);
  else ver8_s16_row_any(p, stride, f, n, out);
}

// Fills `hor` rows [y0, rows) with horizontal 8-tap results (xpos = x - 2)
// and `col` rows with the first column (xpos = -3).
static inline void build_luma_hor(const kvz_pixel *src, int16_t src_stride, int y0, int rows,
                                  int width, const int8_t *f, int16_t *hor, int16_t hor_stride,
                                  int16_t *col)
{
  for (int y = y0; y < rows; ++y) {
    const kvz_pixel *row = src + (intptr_t)src_stride * (y - KVZ_LUMA_FILTER_OFFSET);
    hor8_row(row + (1 - KVZ_LUMA_FILTER_OFFSET), f, width, hor + y * hor_stride);
    col[y] = hor8_px_scalar(row - KVZ_LUMA_FILTER_OFFSET, f);
  }
}

// out[i] = clip(((sum_k f[k] * p[k * stride + i]) + 32) >> 6), 8-bit samples.
static inline void ver8_px_row_any(const kvz_pixel *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
{
  int i = 0;
  for (; i + 8 <= n; i += 8) {
    int16x8_t acc = vdupq_n_s16(0);
    for (int k = 0; k < 8; ++k) {
      int16x8_t v = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(p + (intptr_t)k * stride + i)));
      acc = vmlaq_n_s16(acc, v, f[k]);
    }
    acc = vshrq_n_s16(vaddq_s16(acc, vdupq_n_s16(32)), 6);
    vst1_u8(out + i, vqmovun_s16(acc));
  }
  for (; i < n; ++i) {
    int32_t s = 0;
    for (int k = 0; k < 8; ++k) {
      s += f[k] * p[(intptr_t)k * stride + i];
    }
    out[i] = clip_to_pixel_8bit((s + 32) >> 6);
  }
}

// out[i] = clip(((sum_k f[k] * p[k * stride + i]) >> 6) + 32) >> 6) from int16 intermediates.
static inline void ver8_s16_row_any(const int16_t *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
{
  int i = 0;
  for (; i + 8 <= n; i += 8) {
    int32x4_t lo, hi;
    ver_sum8_s32(p + i, stride, f, 8, &lo, &hi);
    vst1_u8(out + i, ver_finish_pixel_wrap16(lo, hi));
  }
  for (; i < n; ++i) {
    int32_t s = 0;
    for (int k = 0; k < 8; ++k) {
      s += f[k] * p[(intptr_t)k * stride + i];
    }
    int16_t sample = (int16_t)(s >> 6);
    out[i] = clip_to_pixel_8bit(((int32_t)sample + 32) >> 6);
  }
}

// out[i] = clip((h[i] + 32) >> 6) for a horizontal-only intermediate row.
static inline void hshift_clip_row(const int16_t *h, int n, kvz_pixel *out)
{
  int i = 0;
  for (; i + 8 <= n; i += 8) {
    int16x8_t v = vshrq_n_s16(vaddq_s16(vld1q_s16(h + i), vdupq_n_s16(32)), 6);
    vst1_u8(out + i, vqmovun_s16(v));
  }
  for (; i < n; ++i) {
    out[i] = clip_to_pixel_8bit(((int32_t)h[i] + 32) >> 6);
  }
}

// One quarter-pel row. Output x in [s, width) uses hor_row[k * stride + (x - s)].
// When s is 1, output x = 0 comes from the scalar first column `col`.
static inline void qpel_row(const int16_t *hor_row, int16_t hor_stride, const int8_t *f,
                            const int16_t *col, int s, int width, kvz_pixel *out)
{
  if (s) {
    out[0] = qpel_col_sample(col, f);
  }
  ver8_s16_row(hor_row, hor_stride, f, width - s, out + s);
}

static void filter_hpel_blocks_hor_ver_luma_neon(const encoder_control_t * encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  kvz_pixel filtered[4][LCU_LUMA_SIZE],
  int16_t hor_intermediate[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD],
  int8_t fme_level,
  int16_t hor_first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1],
  int8_t hpel_off_x, int8_t hpel_off_y)
{
  (void)encoder; (void)hpel_off_x; (void)hpel_off_y;

  const int8_t *fir0 = kvz_g_luma_filter[0];
  const int8_t *fir2 = kvz_g_luma_filter[2];

  const int16_t dst_stride = LCU_WIDTH;
  const int16_t hor_stride = LCU_WIDTH;
  const int first_row_offset = (KVZ_LUMA_FILTER_OFFSET + 1) * hor_stride;
  const int first_y = fme_level > 1 ? 0 : 1;
  const int rows = height + KVZ_EXT_PADDING_LUMA + 1;

  int16_t *col_pos0 = hor_first_cols[0];
  int16_t *col_pos2 = hor_first_cols[2];

  build_luma_hor(src, src_stride, 0, rows, width, fir0, hor_intermediate[0], hor_stride, col_pos0);
  build_luma_hor(src, src_stride, first_y, rows, width, fir2, hor_intermediate[1], hor_stride, col_pos2);

  // Right: horizontal only
  for (int y = 0; y < height; ++y) {
    hshift_clip_row(hor_intermediate[1] + first_row_offset + y * hor_stride, width, filtered[1] + y * dst_stride);
  }

  // Left: first column from the column buffer, rest copied from Right
  for (int y = 0; y < height; ++y) {
    filtered[0][y * dst_stride] = clip_to_pixel_8bit(((int32_t)col_pos2[y + KVZ_LUMA_FILTER_OFFSET + 1] + 32) >> 6);
    memcpy(&filtered[0][y * dst_stride + 1], &filtered[1][y * dst_stride], (size_t)(width - 1));
  }

  // Top: vertical only
  for (int y = 0; y < height; ++y) {
    const kvz_pixel *row = src + (intptr_t)src_stride * (y - KVZ_LUMA_FILTER_OFFSET);
    ver8_px_row(row + 1, src_stride, fir2, width, filtered[2] + y * dst_stride);
  }

  // Bottom: copy from Top, then filter the last row
  for (int y = 0; y < height - 1; ++y) {
    memcpy(&filtered[3][y * dst_stride], &filtered[2][(y + 1) * dst_stride], (size_t)width);
  }
  {
    const int y = height - 1;
    const int ypos = y - KVZ_LUMA_FILTER_OFFSET;
    ver8_px_row(src + (intptr_t)src_stride * (ypos + 1) + 1, src_stride, fir2, width, filtered[3] + y * dst_stride);
  }
}

static void filter_hpel_blocks_diag_luma_neon(const encoder_control_t * encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  kvz_pixel filtered[4][LCU_LUMA_SIZE],
  int16_t hor_intermediate[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD],
  int8_t fme_level,
  int16_t hor_first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1],
  int8_t hpel_off_x, int8_t hpel_off_y)
{
  (void)encoder; (void)src; (void)src_stride; (void)fme_level; (void)hpel_off_x; (void)hpel_off_y;

  const int8_t *fir2 = kvz_g_luma_filter[2];
  const int16_t dst_stride = LCU_WIDTH;
  const int16_t hor_stride = LCU_WIDTH;
  const int16_t *col_pos2 = hor_first_cols[2];

  // Top-right
  for (int y = 0; y < height; ++y) {
    ver8_s16_row(hor_intermediate[1] + y * hor_stride, hor_stride, fir2, width, filtered[1] + y * dst_stride);
  }

  // Top-left: first column filtered from the column buffer, rest copied from Top-right
  for (int y = 0; y < height; ++y) {
    filtered[0][y * dst_stride] = qpel_col_sample(col_pos2 + y, fir2);
    memcpy(&filtered[0][y * dst_stride + 1], &filtered[1][y * dst_stride], (size_t)(width - 1));
  }

  // Bottom-right: copy from Top-right, then filter the last row
  for (int y = 0; y < height - 1; ++y) {
    memcpy(&filtered[3][y * dst_stride], &filtered[1][(y + 1) * dst_stride], (size_t)width);
  }
  {
    const int y = height - 1;
    ver8_s16_row(hor_intermediate[1] + (y + 1) * hor_stride, hor_stride, fir2, width, filtered[3] + y * dst_stride);
  }

  // Bottom-left: copy from Top-left and Bottom-right, then the last pixel from the column buffer
  for (int y = 0; y < height - 1; ++y) {
    memcpy(&filtered[2][y * dst_stride], &filtered[0][(y + 1) * dst_stride], (size_t)width);
  }
  {
    const int y = height - 1;
    memcpy(&filtered[2][y * dst_stride + 1], &filtered[3][y * dst_stride], (size_t)(width - 1));
    filtered[2][y * dst_stride] = qpel_col_sample(col_pos2 + (y + 1), fir2);
  }
}

static void filter_qpel_blocks_hor_ver_luma_neon(const encoder_control_t * encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  kvz_pixel filtered[4][LCU_LUMA_SIZE],
  int16_t hor_intermediate[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD],
  int8_t fme_level,
  int16_t hor_first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1],
  int8_t hpel_off_x, int8_t hpel_off_y)
{
  (void)encoder; (void)fme_level;

  const int8_t *fir0 = kvz_g_luma_filter[0];
  const int8_t *fir1 = kvz_g_luma_filter[1];
  const int8_t *fir2 = kvz_g_luma_filter[2];
  const int8_t *fir3 = kvz_g_luma_filter[3];

  // Horizontal positions. Positions 0 and 2 were computed by the hpel functions.
  int16_t *hor_pos0 = hor_intermediate[0];
  int16_t *hor_pos2 = hor_intermediate[1];
  int16_t *hor_pos_l = hor_intermediate[3];
  int16_t *hor_pos_r = hor_intermediate[4];
  const int8_t *hor_fir_l = hpel_off_x != 0 ? fir1 : fir3;
  const int8_t *hor_fir_r = hpel_off_x != 0 ? fir3 : fir1;
  int16_t *col_pos_l = hor_first_cols[1];
  int16_t *col_pos_r = hor_first_cols[3];

  const int16_t dst_stride = LCU_WIDTH;
  const int16_t hor_stride = LCU_WIDTH;

  const int16_t *hor_hpel_pos = hpel_off_x != 0 ? hor_pos2 : hor_pos0;
  const int16_t *col_pos_hor = hpel_off_x != 0 ? hor_first_cols[2] : hor_first_cols[0];

  const int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  const int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  const int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  const int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;
  const int sample_off_y = hpel_off_y < 0 ? 0 : 1;
  const int sample_off_x = hpel_off_x > -1 ? 1 : 0;

  const int rows = height + KVZ_EXT_PADDING_LUMA + 1;

  // HORIZONTAL STEP
  build_luma_hor(src, src_stride, 0, rows, width, hor_fir_l, hor_pos_l, hor_stride, col_pos_l);
  build_luma_hor(src, src_stride, 0, rows, width, hor_fir_r, hor_pos_r, hor_stride, col_pos_r);

  const int8_t *ver_fir_l = hpel_off_y != 0 ? fir2 : fir0;
  const int8_t *ver_fir_r = hpel_off_y != 0 ? fir2 : fir0;
  const int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  const int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  // Left QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_l + (y + sample_off_y) * hor_stride, hor_stride, ver_fir_l,
             col_pos_l + y + sample_off_y, !off_x_fir_l, width, filtered[0] + y * dst_stride);
  }

  // Right QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_r + (y + sample_off_y) * hor_stride, hor_stride, ver_fir_r,
             col_pos_r + y + sample_off_y, !off_x_fir_r, width, filtered[1] + y * dst_stride);
  }

  // Top QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_hpel_pos + (y + off_y_fir_t) * hor_stride, hor_stride, ver_fir_t,
             col_pos_hor + y + off_y_fir_t, !sample_off_x, width, filtered[2] + y * dst_stride);
  }

  // Bottom QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_hpel_pos + (y + off_y_fir_b) * hor_stride, hor_stride, ver_fir_b,
             col_pos_hor + y + off_y_fir_b, !sample_off_x, width, filtered[3] + y * dst_stride);
  }
}

static void filter_qpel_blocks_diag_luma_neon(const encoder_control_t * encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  kvz_pixel filtered[4][LCU_LUMA_SIZE],
  int16_t hor_intermediate[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD],
  int8_t fme_level,
  int16_t hor_first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1],
  int8_t hpel_off_x, int8_t hpel_off_y)
{
  (void)encoder; (void)src; (void)src_stride; (void)fme_level;

  const int8_t *fir1 = kvz_g_luma_filter[1];
  const int8_t *fir3 = kvz_g_luma_filter[3];

  // Horizontal positions were computed by the qpel hor_ver call.
  const int16_t *hor_pos_l = hor_intermediate[3];
  const int16_t *hor_pos_r = hor_intermediate[4];
  const int16_t *col_pos_l = hor_first_cols[1];
  const int16_t *col_pos_r = hor_first_cols[3];

  const int16_t dst_stride = LCU_WIDTH;
  const int16_t hor_stride = LCU_WIDTH;

  const int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  const int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  const int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  const int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  const int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  const int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;

  // Top-left QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_l + (y + off_y_fir_t) * hor_stride, hor_stride, ver_fir_t,
             col_pos_l + y + off_y_fir_t, !off_x_fir_l, width, filtered[0] + y * dst_stride);
  }

  // Top-right QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_r + (y + off_y_fir_t) * hor_stride, hor_stride, ver_fir_t,
             col_pos_r + y + off_y_fir_t, !off_x_fir_r, width, filtered[1] + y * dst_stride);
  }

  // Bottom-left QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_l + (y + off_y_fir_b) * hor_stride, hor_stride, ver_fir_b,
             col_pos_l + y + off_y_fir_b, !off_x_fir_l, width, filtered[2] + y * dst_stride);
  }

  // Bottom-right QPEL
  for (int y = 0; y < height; ++y) {
    qpel_row(hor_pos_r + (y + off_y_fir_b) * hor_stride, hor_stride, ver_fir_b,
             col_pos_r + y + off_y_fir_b, !off_x_fir_r, width, filtered[3] + y * dst_stride);
  }
}

int kvz_strategy_register_ipol_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma", "neon", 40, &sample_quarterpel_luma_neon);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma", "neon", 40, &sample_octpel_chroma_neon);
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma_hi", "neon", 40, &sample_quarterpel_luma_hi_neon);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma_hi", "neon", 40, &sample_octpel_chroma_hi_neon);
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_hor_ver_luma", "neon", 40, &filter_hpel_blocks_hor_ver_luma_neon);
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_diag_luma", "neon", 40, &filter_hpel_blocks_diag_luma_neon);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_hor_ver_luma", "neon", 40, &filter_qpel_blocks_hor_ver_luma_neon);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_diag_luma", "neon", 40, &filter_qpel_blocks_diag_luma_neon);
  }

  return success;
}

#else

int kvz_strategy_register_ipol_neon(void* opaque, uint8_t bitdepth)
{
  (void)opaque;
  (void)bitdepth;
  return 1;
}

#endif
