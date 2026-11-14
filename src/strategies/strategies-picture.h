#ifndef STRATEGIES_PICTURE_H_
#define STRATEGIES_PICTURE_H_
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
 * Interface for distortion metric functions.
 */

#include "global.h" // IWYU pragma: keep
#include "inter.h"
#include "kvazaar.h"
#include "encoderstate.h"
#include "strategies/optimized_sad_func_ptr_t.h"


typedef kvz_pixel (*pred_buffer)[32 * 32];

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NxN(suffix, n) \
/* Declare the function in advance, hopefully reducing the probability that the
 * macro expands to something unexpected and silently breaks things. */ \
static cost_pixel_nxn_func satd_ ## n ## x ## n ## _ ## suffix;\
static unsigned satd_ ## n ## x ## n ## _ ## suffix ( \
    const kvz_pixel * const block1, \
    const kvz_pixel * const block2) \
{ \
  unsigned sum = 0; \
  for (unsigned y = 0; y < (n); y += 8) { \
    unsigned row = y * (n); \
    for (unsigned x = 0; x < (n); x += 8) { \
      sum += satd_8x8_subblock_ ## suffix(&block1[row + x], (n), &block2[row + x], (n)); \
    } \
  } \
  return sum >> (KVZ_BIT_DEPTH - 8); \
}


// Function macro for defining hadamard calculating functions for dynamic size
// blocks. They calculate hadamard for integer multiples of 8x8 with the 8x8
// hadamard function.
#define SATD_ANY_SIZE(suffix) \
  static cost_pixel_any_size_func satd_any_size_ ## suffix; \
  static unsigned satd_any_size_ ## suffix ( \
      int width, int height, \
      const kvz_pixel *block1, int stride1, \
      const kvz_pixel *block2, int stride2) \
  { \
    unsigned sum = 0; \
    if (width % 8 != 0) { \
      /* Process the first column using 4x4 blocks. */ \
      for (int y = 0; y < height; y += 4) { \
        sum += kvz_satd_4x4_subblock_ ## suffix(&block1[y * stride1], stride1, \
                                                &block2[y * stride2], stride2); \
      } \
      block1 += 4; \
      block2 += 4; \
      width -= 4; \
    } \
    if (height % 8 != 0) { \
      /* Process the first row using 4x4 blocks. */ \
      for (int x = 0; x < width; x += 4) { \
        sum += kvz_satd_4x4_subblock_ ## suffix(&block1[x], stride1, \
                                                &block2[x], stride2); \
      } \
      block1 += 4 * stride1; \
      block2 += 4 * stride2; \
      height -= 4; \
    } \
    /* The rest can now be processed with 8x8 blocks. */ \
    for (int y = 0; y < height; y += 8) { \
      const kvz_pixel *row1 = &block1[y * stride1]; \
      const kvz_pixel *row2 = &block2[y * stride2]; \
      for (int x = 0; x < width; x += 8) { \
        sum += satd_8x8_subblock_ ## suffix(&row1[x], stride1, \
                                            &row2[x], stride2); \
      } \
    } \
    return sum >> (KVZ_BIT_DEPTH - 8); \
  }

typedef unsigned(reg_sad_func)(const kvz_pixel *const data1, const kvz_pixel *const data2,
  const int width, const int height,
  const unsigned stride1, const unsigned stride2);
typedef unsigned (cost_pixel_nxn_func)(const kvz_pixel *block1, const kvz_pixel *block2);
typedef unsigned (cost_pixel_any_size_func)(
    int width, int height,
    const kvz_pixel *block1, int stride1,
    const kvz_pixel *block2, int stride2
);
typedef void (cost_pixel_nxn_multi_func)(const pred_buffer preds, const kvz_pixel *orig, unsigned num_modes, unsigned *costs_out);
typedef void (cost_pixel_any_size_multi_func)(int width, int height, const kvz_pixel **preds, const int stride, const kvz_pixel *orig, const int orig_stride, unsigned num_modes, unsigned *costs_out, int8_t *valid);

typedef unsigned (pixels_calc_ssd_func)(const kvz_pixel *const ref, const kvz_pixel *const rec, const int ref_stride, const int rec_stride, const int width);
typedef optimized_sad_func_ptr_t (get_optimized_sad_func)(int32_t);
typedef uint32_t (ver_sad_func)(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                int32_t block_width, int32_t block_height,
                                uint32_t pic_stride);
typedef uint32_t (hor_sad_func)(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                int32_t width, int32_t height, uint32_t pic_stride,
                                uint32_t ref_stride, uint32_t left, uint32_t right);

typedef void (inter_recon_bipred_func)(lcu_t * const lcu,
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
  const uint8_t chroma_shift_h);

typedef double (pixel_var_func)(const kvz_pixel *buf, const uint32_t len);

/**
 * \brief SSD between an 8x8 or 16x16 block and a quarter-pel bilinear
 *        interpolation of the reference (used by the temporal filter).
 *
 * The prediction is ((4-fx)(4-fy)r00 + fx(4-fy)r10 + (4-fx)fy r01 + fx fy r11 + 8) >> 4.
 * The reference must have one readable column and row beyond the block.
 */
typedef unsigned (bilinear_ssd_func)(const kvz_pixel *org, int org_stride,
                                     const kvz_pixel *ref, int ref_stride,
                                     int size, int fx, int fy);

/**
 * \brief kvz_bilinear_ssd for several quarter-pel motion vectors.
 *
 * ref points to the colocated block. For vector i, the block at integer
 * offset (mv_x[i] >> 2, mv_y[i] >> 2) is interpolated with the fraction
 * (mv_x[i] & 3, mv_y[i] & 3). Integer positions give the plain SSD (scaled
 * to 8 bits like kvz_pixels_calc_ssd).
 */
typedef void (bilinear_ssd_multi_func)(const kvz_pixel *org, int org_stride,
                                       const kvz_pixel *ref, int ref_stride,
                                       int size, const int *mv_x, const int *mv_y,
                                       int num, unsigned *out);

/**
 * \brief Bilinear interpolation of a block (used by the temporal filter).
 *
 * dst[x] = (w00 r0[x] + w10 r0[x + 1] + w01 r1[x] + w11 r1[x + 1] +
 *           (1 << (shift - 1))) >> shift
 * where r0 and r1 are consecutive rows of ref and the weights sum to
 * 1 << shift (shift 1..6). The reference must have one readable column and
 * row beyond the block.
 */
typedef void (tf_bilinear_pred_func)(kvz_pixel *dst, int dst_stride,
                                     const kvz_pixel *ref, int ref_stride,
                                     int width, int height,
                                     int w00, int w10, int w01, int w11, int shift);

/**
 * \brief Statistics of the prediction error of a block for the temporal
 *        filter.
 *
 * With d = org - pred:
 * out[0] = sum of org, out[1] = sum of org^2, out[2] = sum of d^2,
 * and over the pixels that have a right and a lower neighbour in the block:
 * out[3] = sum of d^2, out[4] = sum of (d_right - d)^2 + (d_below - d)^2.
 */
typedef void (tf_block_stats_func)(const kvz_pixel *org, int org_stride,
                                   const kvz_pixel *pred, int pred_stride,
                                   int width, int height, int64_t out[5]);

/**
 * \brief Bilateral temporal filtering of a block.
 *
 * For each pixel, with o = src[x] and p_r = preds[r][x]:
 *   sum = o, wsum = 1, for each r: w = weights[r] * luts[r][p_r - o],
 *   sum += w * p_r, wsum += w;
 *   dst[x] = CLIP(0, max_val, (int)(sum / wsum + 0.5f))
 * evaluated in this order in single precision. luts[r] is indexed with the
 * difference of the pixel values, from -max_val to max_val.
 */
typedef void (tf_filter_block_func)(kvz_pixel *dst, int dst_stride,
                                    const kvz_pixel *src, int src_stride,
                                    const kvz_pixel *const *preds, int pred_stride,
                                    int num_refs, const float *weights,
                                    const float *const *luts,
                                    int width, int height, int max_val);

// Declare function pointers.
extern reg_sad_func * kvz_reg_sad;

extern cost_pixel_nxn_func * kvz_sad_4x4;
extern cost_pixel_nxn_func * kvz_sad_8x8;
extern cost_pixel_nxn_func * kvz_sad_16x16;
extern cost_pixel_nxn_func * kvz_sad_32x32;
extern cost_pixel_nxn_func * kvz_sad_64x64;

extern cost_pixel_nxn_func * kvz_satd_4x4;
extern cost_pixel_nxn_func * kvz_satd_8x8;
extern cost_pixel_nxn_func * kvz_satd_16x16;
extern cost_pixel_nxn_func * kvz_satd_32x32;
extern cost_pixel_nxn_func * kvz_satd_64x64;
extern cost_pixel_any_size_func *kvz_satd_any_size;

extern cost_pixel_nxn_multi_func * kvz_sad_4x4_dual;
extern cost_pixel_nxn_multi_func * kvz_sad_8x8_dual;
extern cost_pixel_nxn_multi_func * kvz_sad_16x16_dual;
extern cost_pixel_nxn_multi_func * kvz_sad_32x32_dual;
extern cost_pixel_nxn_multi_func * kvz_sad_64x64_dual;

extern cost_pixel_nxn_multi_func * kvz_satd_4x4_dual;
extern cost_pixel_nxn_multi_func * kvz_satd_8x8_dual;
extern cost_pixel_nxn_multi_func * kvz_satd_16x16_dual;
extern cost_pixel_nxn_multi_func * kvz_satd_32x32_dual;
extern cost_pixel_nxn_multi_func * kvz_satd_64x64_dual;

extern cost_pixel_any_size_multi_func *kvz_satd_any_size_quad;

extern pixels_calc_ssd_func *kvz_pixels_calc_ssd;

extern inter_recon_bipred_func * kvz_bipred_average;

extern get_optimized_sad_func *kvz_get_optimized_sad;
extern ver_sad_func *kvz_ver_sad;
extern hor_sad_func *kvz_hor_sad;

extern pixel_var_func *kvz_pixel_var;
extern bilinear_ssd_func *kvz_bilinear_ssd;
extern bilinear_ssd_multi_func *kvz_bilinear_ssd_multi;
extern tf_bilinear_pred_func *kvz_tf_bilinear_pred;
extern tf_block_stats_func *kvz_tf_block_stats;
extern tf_filter_block_func *kvz_tf_filter_block;

int kvz_strategy_register_picture(void* opaque, uint8_t bitdepth);
cost_pixel_nxn_func * kvz_pixels_get_satd_func(unsigned n);
cost_pixel_nxn_func * kvz_pixels_get_sad_func(unsigned n);
cost_pixel_nxn_multi_func * kvz_pixels_get_satd_dual_func(unsigned n);
cost_pixel_nxn_multi_func * kvz_pixels_get_sad_dual_func(unsigned n);

#define STRATEGIES_PICTURE_EXPORTS \
  {"reg_sad", (void**) &kvz_reg_sad}, \
  {"sad_4x4", (void**) &kvz_sad_4x4}, \
  {"sad_8x8", (void**) &kvz_sad_8x8}, \
  {"sad_16x16", (void**) &kvz_sad_16x16}, \
  {"sad_32x32", (void**) &kvz_sad_32x32}, \
  {"sad_64x64", (void**) &kvz_sad_64x64}, \
  {"satd_4x4", (void**) &kvz_satd_4x4}, \
  {"satd_8x8", (void**) &kvz_satd_8x8}, \
  {"satd_16x16", (void**) &kvz_satd_16x16}, \
  {"satd_32x32", (void**) &kvz_satd_32x32}, \
  {"satd_64x64", (void**) &kvz_satd_64x64}, \
  {"satd_any_size", (void**) &kvz_satd_any_size}, \
  {"sad_4x4_dual", (void**) &kvz_sad_4x4_dual}, \
  {"sad_8x8_dual", (void**) &kvz_sad_8x8_dual}, \
  {"sad_16x16_dual", (void**) &kvz_sad_16x16_dual}, \
  {"sad_32x32_dual", (void**) &kvz_sad_32x32_dual}, \
  {"sad_64x64_dual", (void**) &kvz_sad_64x64_dual}, \
  {"satd_4x4_dual", (void**) &kvz_satd_4x4_dual}, \
  {"satd_8x8_dual", (void**) &kvz_satd_8x8_dual}, \
  {"satd_16x16_dual", (void**) &kvz_satd_16x16_dual}, \
  {"satd_32x32_dual", (void**) &kvz_satd_32x32_dual}, \
  {"satd_64x64_dual", (void**) &kvz_satd_64x64_dual}, \
  {"satd_any_size_quad", (void**) &kvz_satd_any_size_quad}, \
  {"pixels_calc_ssd", (void**) &kvz_pixels_calc_ssd}, \
  {"bipred_average", (void**) &kvz_bipred_average}, \
  {"get_optimized_sad", (void**) &kvz_get_optimized_sad}, \
  {"ver_sad", (void**) &kvz_ver_sad}, \
  {"hor_sad", (void**) &kvz_hor_sad}, \
  {"pixel_var", (void**) &kvz_pixel_var}, \
  {"bilinear_ssd", (void**) &kvz_bilinear_ssd}, \
  {"bilinear_ssd_multi", (void**) &kvz_bilinear_ssd_multi}, \
  {"tf_bilinear_pred", (void**) &kvz_tf_bilinear_pred}, \
  {"tf_block_stats", (void**) &kvz_tf_block_stats}, \
  {"tf_filter_block", (void**) &kvz_tf_filter_block}, \



#endif //STRATEGIES_PICTURE_H_
