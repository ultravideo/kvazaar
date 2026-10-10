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

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "search_inter.h"
#include "strategies/neon/ipol-neon.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-ipol.h"
#include "strategyselector.h"

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

extern int8_t kvz_g_luma_filter[4][8];
extern int8_t kvz_g_chroma_filter[8][4];

static inline kvz_pixel clip_to_pixel_8bit(int32_t v)
{
  return (kvz_pixel)(v < 0 ? 0 : (v > 255 ? 255 : v));
}

/*
 * Horizontal pass. `base` points to the sample at (-pad, -pad) relative to the
 * block origin, so output (y, x) uses the taps base[y * stride + x + k].
 * Writes `rows` rows of `width` intermediate values into dst.
 */
static inline void hor_pass_s16(const kvz_pixel *base, int16_t src_stride,
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

static void ver_pass_pixel(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                           int width, int height,
                           kvz_pixel *dst, int16_t dst_stride)
{
  for (int y = 0; y < height; ++y) {
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int32x4_t lo, hi;
      ver_sum8_s32(&hor[y][x], LCU_WIDTH, f, taps, &lo, &hi);
      vst1_u8(&dst[y * dst_stride + x], ver_finish_pixel(lo, hi));
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += f[k] * hor[y + k][x];
      }
      dst[y * dst_stride + x] = clip_to_pixel_8bit(((sum >> 6) + 32) >> 6);
    }
  }
}

static void ver_pass_hi(const int16_t hor[][LCU_WIDTH], const int8_t *f, int taps,
                        int width, int height,
                        int16_t *dst, int16_t dst_stride)
{
  for (int y = 0; y < height; ++y) {
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      int32x4_t lo, hi;
      ver_sum8_s32(&hor[y][x], LCU_WIDTH, f, taps, &lo, &hi);
      /* (sum >> 6) stored as int16, same as the generic truncating store. */
      vst1q_s16(&dst[y * dst_stride + x],
                vcombine_s16(vshrn_n_s32(lo, 6), vshrn_n_s32(hi, 6)));
    }
    for (; x < width; ++x) {
      int32_t sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += f[k] * hor[y + k][x];
      }
      dst[y * dst_stride + x] = (int16_t)(sum >> 6);
    }
  }
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
static inline void hor8_row(const kvz_pixel *p, const int8_t *f, int n, int16_t *out)
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
static inline void ver8_px_row(const kvz_pixel *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
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
static inline void ver8_s16_row(const int16_t *p, int16_t stride, const int8_t *f, int n, kvz_pixel *out)
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
