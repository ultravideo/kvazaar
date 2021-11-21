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

#include "strategies/generic/picture-generic.h"

#include <stdlib.h>

#include "strategies/strategies-picture.h"
#include "strategyselector.h"

// Function to clip int16_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel kvz_fast_clip_16bit_to_pixel(int16_t value)
{
  // Ensure that compiler generates arithmetic shift from ">>" 
#if defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)

  if (value & ~PIXEL_MAX) {
    int16_t temp = (-value) >> 15;
#if KVZ_BIT_DEPTH == 10
    temp &= PIXEL_MAX;
#endif
    return temp;
  }
  else {
    return value;
  }
#else
  return CLIP(PIXEL_MIN, PIXEL_MAX, value);
#endif
}

// Function to clip int32_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel kvz_fast_clip_32bit_to_pixel(int32_t value)
{
  // Ensure that compiler generates arithmetic shift from ">>" 
#if defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)

  if (value & ~PIXEL_MAX) {
    int32_t temp = (-value) >> 31;
#if KVZ_BIT_DEPTH == 10
    temp &= PIXEL_MAX;
#endif
    return temp;
  }
  else {
    return value;
  }
#else
  return CLIP(PIXEL_MIN, PIXEL_MAX, value);
#endif
}

/**
 * \brief Calculate Sum of Absolute Differences (SAD)
 *
 * Calculate Sum of Absolute Differences (SAD) between two rectangular regions
 * located in arbitrary points in the picture.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param stride  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned reg_sad_generic(const kvz_pixel * const data1, const kvz_pixel * const data2,
                         const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  int y, x;
  unsigned sad = 0;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }

  return sad;
}

/**
 * \brief  Transform differences between two 4x4 blocks.
 * From HM 13.0
 */
static int32_t hadamard_4x4_generic(int32_t diff[4*4])
{
  int32_t m[4 * 4];
  m[0] = diff[0] + diff[12];
  m[1] = diff[1] + diff[13];
  m[2] = diff[2] + diff[14];
  m[3] = diff[3] + diff[15];
  m[4] = diff[4] + diff[8];
  m[5] = diff[5] + diff[9];
  m[6] = diff[6] + diff[10];
  m[7] = diff[7] + diff[11];
  m[8] = diff[4] - diff[8];
  m[9] = diff[5] - diff[9];
  m[10] = diff[6] - diff[10];
  m[11] = diff[7] - diff[11];
  m[12] = diff[0] - diff[12];
  m[13] = diff[1] - diff[13];
  m[14] = diff[2] - diff[14];
  m[15] = diff[3] - diff[15];

  int32_t d[4 * 4];
  d[0] = m[0] + m[4];
  d[1] = m[1] + m[5];
  d[2] = m[2] + m[6];
  d[3] = m[3] + m[7];
  d[4] = m[8] + m[12];
  d[5] = m[9] + m[13];
  d[6] = m[10] + m[14];
  d[7] = m[11] + m[15];
  d[8] = m[0] - m[4];
  d[9] = m[1] - m[5];
  d[10] = m[2] - m[6];
  d[11] = m[3] - m[7];
  d[12] = m[12] - m[8];
  d[13] = m[13] - m[9];
  d[14] = m[14] - m[10];
  d[15] = m[15] - m[11];

  m[0] = d[0] + d[3];
  m[1] = d[1] + d[2];
  m[2] = d[1] - d[2];
  m[3] = d[0] - d[3];
  m[4] = d[4] + d[7];
  m[5] = d[5] + d[6];
  m[6] = d[5] - d[6];
  m[7] = d[4] - d[7];
  m[8] = d[8] + d[11];
  m[9] = d[9] + d[10];
  m[10] = d[9] - d[10];
  m[11] = d[8] - d[11];
  m[12] = d[12] + d[15];
  m[13] = d[13] + d[14];
  m[14] = d[13] - d[14];
  m[15] = d[12] - d[15];

  d[0] = m[0] + m[1];
  d[1] = m[0] - m[1];
  d[2] = m[2] + m[3];
  d[3] = m[3] - m[2];
  d[4] = m[4] + m[5];
  d[5] = m[4] - m[5];
  d[6] = m[6] + m[7];
  d[7] = m[7] - m[6];
  d[8] = m[8] + m[9];
  d[9] = m[8] - m[9];
  d[10] = m[10] + m[11];
  d[11] = m[11] - m[10];
  d[12] = m[12] + m[13];
  d[13] = m[12] - m[13];
  d[14] = m[14] + m[15];
  d[15] = m[15] - m[14];

  int32_t satd = 0;
  for (int i = 0; i < 16; i++) {
    satd += abs(d[i]);
  }
  satd = ((satd + 1) >> 1);

  return satd;
}

/**
 * \brief  Calculate SATD between two 4x4 blocks.
 */
static unsigned satd_4x4_generic(const kvz_pixel *piOrg, const kvz_pixel *piCur)
{
  int32_t diff[4 * 4];
  for (int i = 0; i < 4 * 4; i++) {
    diff[i] = piOrg[i] - piCur[i];
  }
  return hadamard_4x4_generic(diff);
}

/**
* \brief  Calculate SATD between two 4x4 blocks inside bigger arrays.
*/
unsigned kvz_satd_4x4_subblock_generic(const kvz_pixel * buf1,
                                       const int32_t     stride1,
                                       const kvz_pixel * buf2,
                                       const int32_t     stride2)
{
  int32_t diff[4 * 4];
  for (int y = 0; y < 4; y++) {
    for (int x = 0; x < 4; x++) {
      diff[x + y * 4] = buf1[x + y * stride1] - buf2[x + y * stride2];
    }
  }
  return hadamard_4x4_generic(diff);
}

void kvz_satd_4x4_subblock_quad_generic(const kvz_pixel *preds[4],
                                       const int stride,
                                       const kvz_pixel *orig,
                                       const int orig_stride,
                                       unsigned costs[4])
{
  int32_t diff[4][4 * 4];
  for (int y = 0; y < 4; y++) {
    for (int x = 0; x < 4; x++) {
      diff[0][x + y * 4] = orig[x + y * orig_stride] - preds[0][x + y * stride];
      diff[1][x + y * 4] = orig[x + y * orig_stride] - preds[1][x + y * stride];
      diff[2][x + y * 4] = orig[x + y * orig_stride] - preds[2][x + y * stride];
      diff[3][x + y * 4] = orig[x + y * orig_stride] - preds[3][x + y * stride];
    }
  }

  costs[0] = hadamard_4x4_generic(diff[0]);
  costs[1] = hadamard_4x4_generic(diff[1]);
  costs[2] = hadamard_4x4_generic(diff[2]);
  costs[3] = hadamard_4x4_generic(diff[3]);
}

/**
* \brief  Calculate SATD between two 8x8 blocks inside bigger arrays.
*/
static unsigned satd_8x8_subblock_generic(const kvz_pixel * piOrg, const int32_t iStrideOrg,
  const kvz_pixel * piCur, const int32_t iStrideCur)
{
  int32_t k, i, j, jj, sad = 0;
  int32_t diff[64], m1[8][8], m2[8][8], m3[8][8];

  for (k = 0; k < 64; k += 8) {
    diff[k + 0] = piOrg[0] - piCur[0];
    diff[k + 1] = piOrg[1] - piCur[1];
    diff[k + 2] = piOrg[2] - piCur[2];
    diff[k + 3] = piOrg[3] - piCur[3];
    diff[k + 4] = piOrg[4] - piCur[4];
    diff[k + 5] = piOrg[5] - piCur[5];
    diff[k + 6] = piOrg[6] - piCur[6];
    diff[k + 7] = piOrg[7] - piCur[7];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  // horizontal
  for (j = 0; j < 8; ++j) {
    jj = j << 3;
    m2[j][0] = diff[jj] + diff[jj + 4];
    m2[j][1] = diff[jj + 1] + diff[jj + 5];
    m2[j][2] = diff[jj + 2] + diff[jj + 6];
    m2[j][3] = diff[jj + 3] + diff[jj + 7];
    m2[j][4] = diff[jj] - diff[jj + 4];
    m2[j][5] = diff[jj + 1] - diff[jj + 5];
    m2[j][6] = diff[jj + 2] - diff[jj + 6];
    m2[j][7] = diff[jj + 3] - diff[jj + 7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  // vertical
  for (i = 0; i < 8; ++i) {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 64; ++i) {
    sad += abs(((int*)m2)[i]);
  }

  sad = (sad + 2) >> 2;

  return sad;
}

static void satd_8x8_subblock_quad_generic(const kvz_pixel **preds,
                                       const int stride,
                                       const kvz_pixel *orig,
                                       const int orig_stride,
                                       unsigned *costs)
{
  costs[0] = satd_8x8_subblock_generic(orig, orig_stride, preds[0], stride);
  costs[1] = satd_8x8_subblock_generic(orig, orig_stride, preds[1], stride);
  costs[2] = satd_8x8_subblock_generic(orig, orig_stride, preds[2], stride);
  costs[3] = satd_8x8_subblock_generic(orig, orig_stride, preds[3], stride);
}

// These macros define sadt_16bit_NxN for N = 8, 16, 32, 64
SATD_NxN(generic,  8)
SATD_NxN(generic, 16)
SATD_NxN(generic, 32)
SATD_NxN(generic, 64)
SATD_ANY_SIZE(generic)


// Declare these functions to make sure the signature of the macro matches.
static cost_pixel_nxn_multi_func satd_4x4_dual_generic;
static cost_pixel_nxn_multi_func satd_8x8_dual_generic;
static cost_pixel_nxn_multi_func satd_16x16_dual_generic;
static cost_pixel_nxn_multi_func satd_32x32_dual_generic;
static cost_pixel_nxn_multi_func satd_64x64_dual_generic;

#define SATD_DUAL_NXN(n, pixel_type) \
static void satd_ ## n ## x ## n ## _dual_generic( \
  const pred_buffer preds, const pixel_type * const orig, unsigned num_modes, unsigned *costs_out) \
{ \
  unsigned x, y; \
  unsigned sum = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  sum += satd_8x8_subblock_generic(&preds[0][row + x], (n), &orig[row + x], (n)); \
  } \
  } \
  costs_out[0] = sum>>(KVZ_BIT_DEPTH-8); \
  \
  sum = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  sum += satd_8x8_subblock_generic(&preds[1][row + x], (n), &orig[row + x], (n)); \
  } \
  } \
  costs_out[1] = sum>>(KVZ_BIT_DEPTH-8); \
}

static void satd_4x4_dual_generic(const pred_buffer preds, const kvz_pixel * const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = satd_4x4_generic(orig, preds[0]);
  costs_out[1] = satd_4x4_generic(orig, preds[1]);
}

SATD_DUAL_NXN(8, kvz_pixel)
SATD_DUAL_NXN(16, kvz_pixel)
SATD_DUAL_NXN(32, kvz_pixel)
SATD_DUAL_NXN(64, kvz_pixel)

#define SATD_ANY_SIZE_MULTI_GENERIC(suffix, num_parallel_blocks) \
  static cost_pixel_any_size_multi_func satd_any_size_## suffix; \
  static void satd_any_size_ ## suffix ( \
      int width, int height, \
      const kvz_pixel **preds, \
      const int stride, \
      const kvz_pixel *orig, \
      const int orig_stride, \
      unsigned num_modes, \
      unsigned *costs_out, \
      int8_t *valid) \
  { \
    unsigned sums[num_parallel_blocks] = { 0 }; \
    const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };\
    const kvz_pixel *orig_ptr = orig; \
    costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0; \
    if (width % 8 != 0) { \
      /* Process the first column using 4x4 blocks. */ \
      for (int y = 0; y < height; y += 4) { \
        kvz_satd_4x4_subblock_ ## suffix(preds, stride, orig, orig_stride, sums); \
            } \
      orig_ptr += 4; \
      for(int blk = 0; blk < num_parallel_blocks; ++blk){\
        pred_ptrs[blk] += 4; \
            }\
      width -= 4; \
            } \
    if (height % 8 != 0) { \
      /* Process the first row using 4x4 blocks. */ \
      for (int x = 0; x < width; x += 4 ) { \
        kvz_satd_4x4_subblock_ ## suffix(pred_ptrs, stride, orig_ptr, orig_stride, sums); \
            } \
      orig_ptr += 4 * orig_stride; \
      for(int blk = 0; blk < num_parallel_blocks; ++blk){\
        pred_ptrs[blk] += 4 * stride; \
            }\
      height -= 4; \
        } \
    /* The rest can now be processed with 8x8 blocks. */ \
    for (int y = 0; y < height; y += 8) { \
      orig_ptr = &orig[y * orig_stride]; \
      pred_ptrs[0] = &preds[0][y * stride]; \
      pred_ptrs[1] = &preds[1][y * stride]; \
      pred_ptrs[2] = &preds[2][y * stride]; \
      pred_ptrs[3] = &preds[3][y * stride]; \
      for (int x = 0; x < width; x += 8) { \
        satd_8x8_subblock_ ## suffix(pred_ptrs, stride, orig_ptr, orig_stride, sums); \
        orig_ptr += 8; \
        pred_ptrs[0] += 8; \
        pred_ptrs[1] += 8; \
        pred_ptrs[2] += 8; \
        pred_ptrs[3] += 8; \
        costs_out[0] += sums[0]; \
        costs_out[1] += sums[1]; \
        costs_out[2] += sums[2]; \
        costs_out[3] += sums[3]; \
      } \
    } \
    for(int i = 0; i < num_parallel_blocks; ++i){\
      costs_out[i] = costs_out[i] >> (KVZ_BIT_DEPTH - 8);\
    } \
    return; \
  }

SATD_ANY_SIZE_MULTI_GENERIC(quad_generic, 4)

// Function macro for defining SAD calculating functions
// for fixed size blocks.
#define SAD_NXN(n, pixel_type) \
static unsigned sad_ ##  n ## x ## n ## _generic( \
  const pixel_type * const block1, const pixel_type * const block2) \
{ \
  unsigned i; \
  unsigned sum = 0; \
  for (i = 0; i < (n)*(n); ++i) { \
  sum += abs(block1[i] - block2[i]); \
  } \
  return sum>>(KVZ_BIT_DEPTH-8); \
}

// Declare these functions to make sure the signature of the macro matches.
static cost_pixel_nxn_func sad_4x4_generic;
static cost_pixel_nxn_func sad_8x8_generic;
static cost_pixel_nxn_func sad_16x16_generic;
static cost_pixel_nxn_func sad_32x32_generic;
static cost_pixel_nxn_func sad_64x64_generic;

// These macros define sad_16bit_nxn functions for n = 4, 8, 16, 32, 64
// with function signatures of cost_16bit_nxn_func.
// They are used through get_pixel_sad_func.
SAD_NXN(4, kvz_pixel)
SAD_NXN(8, kvz_pixel)
SAD_NXN(16, kvz_pixel)
SAD_NXN(32, kvz_pixel)
SAD_NXN(64, kvz_pixel)

// Declare these functions to make sure the signature of the macro matches.
static cost_pixel_nxn_multi_func sad_4x4_dual_generic;
static cost_pixel_nxn_multi_func sad_8x8_dual_generic;
static cost_pixel_nxn_multi_func sad_16x16_dual_generic;
static cost_pixel_nxn_multi_func sad_32x32_dual_generic;
static cost_pixel_nxn_multi_func sad_64x64_dual_generic;

// Function macro for defining SAD calculating functions
// for fixed size blocks.
#define SAD_DUAL_NXN(n, pixel_type) \
static void sad_ ##  n ## x ## n ## _dual_generic( \
  const pred_buffer preds, const pixel_type * const orig, unsigned num_modes, unsigned *costs_out) \
{ \
  unsigned i; \
  unsigned sum = 0; \
  for (i = 0; i < (n)*(n); ++i) { \
  sum += abs(preds[0][i] - orig[i]); \
  } \
  costs_out[0] = sum>>(KVZ_BIT_DEPTH-8); \
  \
  sum = 0; \
  for (i = 0; i < (n)*(n); ++i) { \
  sum += abs(preds[1][i] - orig[i]); \
  } \
  costs_out[1] = sum>>(KVZ_BIT_DEPTH-8); \
}

SAD_DUAL_NXN(4, kvz_pixel)
SAD_DUAL_NXN(8, kvz_pixel)
SAD_DUAL_NXN(16, kvz_pixel)
SAD_DUAL_NXN(32, kvz_pixel)
SAD_DUAL_NXN(64, kvz_pixel)

static unsigned pixels_calc_ssd_generic(const kvz_pixel *const ref, const kvz_pixel *const rec,
                 const int ref_stride, const int rec_stride,
                 const int width)
{
  int ssd = 0;
  int y, x;

  for (y = 0; y < width; ++y) {
    for (x = 0; x < width; ++x) {
      int diff = ref[x + y * ref_stride] - rec[x + y * rec_stride];
      ssd += diff * diff;
    }
  }

  return ssd >> (2*(KVZ_BIT_DEPTH-8));
}

static void bipred_average_px_px(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t offset = 1 << (shift - 1);

  for (int i = 0; i < pu_w * pu_h; ++i)
  {
    int y = i / pu_w;
    int x = i % pu_w;
    int16_t sample_L0 = px_L0[i] << (14 - KVZ_BIT_DEPTH);
    int16_t sample_L1 = px_L1[i] << (14 - KVZ_BIT_DEPTH);
    int32_t rounded = (sample_L0 + sample_L1 + offset) >> shift;
    dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
  }
}

static void bipred_average_im_im(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t offset = 1 << (shift - 1);

  for (int i = 0; i < pu_w * pu_h; ++i)
  {
    int y = i / pu_w;
    int x = i % pu_w;
    int16_t sample_L0 = im_L0[i];
    int16_t sample_L1 = im_L1[i];
    int32_t rounded = (sample_L0 + sample_L1 + offset) >> shift;
    dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
  }
}

static void bipred_average_px_im(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t offset = 1 << (shift - 1);

  for (int i = 0; i < pu_w * pu_h; ++i)
  {
    int y = i / pu_w;
    int x = i % pu_w;
    int16_t sample_px = px[i] << (14 - KVZ_BIT_DEPTH);
    int16_t sample_im = im[i];
    int32_t rounded = (sample_px + sample_im + offset) >> shift;
    dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
  }
}

static void bipred_average_generic(lcu_t *const lcu,
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
  const bool predict_chroma) {

  //After reconstruction, merge the predictors by taking an average of each pixel
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
    unsigned pb_offset = SUB_SCU(pu_y) / 2 * LCU_WIDTH_C + SUB_SCU(pu_x) / 2;
    unsigned pb_w = pu_w / 2;
    unsigned pb_h = pu_h / 2;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_px_px(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_px_px(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, LCU_WIDTH_C);

    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_im_im(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_im_im(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, LCU_WIDTH_C);

    } else {
      kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_px_im(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_px_im(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, LCU_WIDTH_C);
    }
  }
}


static optimized_sad_func_ptr_t get_optimized_sad_generic(int32_t width)
{
  return NULL;
}

/**
 * \brief Vertically interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static uint32_t ver_sad_generic(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                int block_width, int block_height, unsigned pic_stride)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref_data[x]);
    }
  }

  return sad;
}

/**
 * \brief Horizontally interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width   Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned hor_sad(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride, unsigned ref_stride)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref_data[y * ref_stride]);
    }
  }

  return sad;
}


static uint32_t hor_sad_generic(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                int32_t width, int32_t height, uint32_t pic_stride,
                                uint32_t ref_stride, uint32_t left, uint32_t right)
{
  uint32_t result = 0;
  if (left) {
    result += hor_sad    (pic_data, ref_data + left, left,
                          height, pic_stride, ref_stride);

    result += kvz_reg_sad(pic_data + left, ref_data + left, width - left,
                          height, pic_stride, ref_stride);
  } else if (right) {
    result += kvz_reg_sad(pic_data, ref_data, width - right,
                          height, pic_stride, ref_stride);

    result += hor_sad    (pic_data + width - right,
                          ref_data + width - right - 1,
                          right, height, pic_stride, ref_stride);
  } else {
    result += kvz_reg_sad(pic_data, ref_data, width,
                          height, pic_stride, ref_stride);
  }
  return result;
}

// Calculate pixel value variance. Takes in arrays of kvz_pixel
static double pixel_var_generic(const kvz_pixel *arr, const uint32_t len)
{
  double var = 0;
  double arr_mean = 0;

  // Calculate array mean
  int i = 0;
  double sum = 0;

  for (; i < len; ++i) {
    sum += arr[i];
  }
  arr_mean = sum / (double)len;

  // Calculate array variance
  for (i = 0; i < len; ++i) {
    double tmp = (double)arr[i] - arr_mean;
    var += tmp*tmp;
  }

  var /= len;

  return var;
}

int kvz_strategy_register_picture_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "reg_sad", "generic", 0, &reg_sad_generic);

  success &= kvz_strategyselector_register(opaque, "sad_4x4", "generic", 0, &sad_4x4_generic);
  success &= kvz_strategyselector_register(opaque, "sad_8x8", "generic", 0, &sad_8x8_generic);
  success &= kvz_strategyselector_register(opaque, "sad_16x16", "generic", 0, &sad_16x16_generic);
  success &= kvz_strategyselector_register(opaque, "sad_32x32", "generic", 0, &sad_32x32_generic);
  success &= kvz_strategyselector_register(opaque, "sad_64x64", "generic", 0, &sad_64x64_generic);

  success &= kvz_strategyselector_register(opaque, "satd_4x4", "generic", 0, &satd_4x4_generic);
  success &= kvz_strategyselector_register(opaque, "satd_8x8", "generic", 0, &satd_8x8_generic);
  success &= kvz_strategyselector_register(opaque, "satd_16x16", "generic", 0, &satd_16x16_generic);
  success &= kvz_strategyselector_register(opaque, "satd_32x32", "generic", 0, &satd_32x32_generic);
  success &= kvz_strategyselector_register(opaque, "satd_64x64", "generic", 0, &satd_64x64_generic);

  success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "generic", 0, &sad_4x4_dual_generic);
  success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "generic", 0, &sad_8x8_dual_generic);
  success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "generic", 0, &sad_16x16_dual_generic);
  success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "generic", 0, &sad_32x32_dual_generic);
  success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "generic", 0, &sad_64x64_dual_generic);

  success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "generic", 0, &satd_4x4_dual_generic);
  success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "generic", 0, &satd_8x8_dual_generic);
  success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "generic", 0, &satd_16x16_dual_generic);
  success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "generic", 0, &satd_32x32_dual_generic);
  success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "generic", 0, &satd_64x64_dual_generic);
  success &= kvz_strategyselector_register(opaque, "satd_any_size", "generic", 0, &satd_any_size_generic);
  success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "generic", 0, &satd_any_size_quad_generic);

  success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "generic", 0, &pixels_calc_ssd_generic);
  success &= kvz_strategyselector_register(opaque, "bipred_average", "generic", 0, &bipred_average_generic);

  success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "generic", 0, &get_optimized_sad_generic);
  success &= kvz_strategyselector_register(opaque, "ver_sad", "generic", 0, &ver_sad_generic);
  success &= kvz_strategyselector_register(opaque, "hor_sad", "generic", 0, &hor_sad_generic);

  success &= kvz_strategyselector_register(opaque, "pixel_var", "generic", 0, &pixel_var_generic);

  return success;
}
