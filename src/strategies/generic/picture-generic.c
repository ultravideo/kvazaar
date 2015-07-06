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

/*
 * \file
 */

#include <stdlib.h>

#include "strategyselector.h"

// Function to clip int16_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel fast_clip_16bit_to_pixel(int16_t value)
{
  // Ensure that compiler generates arithmetic shift from ">>" 
#if defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)

  if (value & ~PIXEL_MAX) {
    int16_t temp = (-value) >> 15;
#if BITDEPTH == 10
    temp &= PIXEL_MAX;
#endif
    return temp;
  }
  else {
    return value;
  }
#else
  CLIP(PIXEL_MIN, PIXEL_MAX, value);
#endif
}

// Function to clip int32_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel fast_clip_32bit_to_pixel(int32_t value)
{
  // Ensure that compiler generates arithmetic shift from ">>" 
#if defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)

  if (value & ~PIXEL_MAX) {
    int32_t temp = (-value) >> 31;
#if BITDEPTH == 10
    temp &= PIXEL_MAX;
#endif
    return temp;
  }
  else {
    return value;
  }
#else
  CLIP(PIXEL_MIN, PIXEL_MAX, value);
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
 * \brief  Calculate SATD between two 4x4 blocks inside bigger arrays.
 * From HM 13.0
 */
static unsigned satd_8bit_4x4_generic(const kvz_pixel *piOrg, const kvz_pixel *piCur)
{
  int32_t k, satd = 0, diff[16], m[16], d[16];
  for (k = 0; k < 16; ++k) {
    diff[k] = piOrg[k] - piCur[k];
  }

  /*===== hadamard transform =====*/
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

  for (k = 0; k<16; ++k) {
    satd += abs(d[k]);
  }
  satd = ((satd + 1) >> 1);

  return satd;
}

/**
* \brief  Calculate SATD between two 8x8 blocks inside bigger arrays.
*/
unsigned satd_16bit_8x8_general(const kvz_pixel * piOrg, const int32_t iStrideOrg,
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

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NXN(n, pixel_type, suffix) \
unsigned satd_ ## suffix ## _ ## n ## x ## n ## _generic( \
  const pixel_type * const block1, const pixel_type * const block2) \
{ \
  unsigned x, y; \
  unsigned sum = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  sum += satd_16bit_8x8_general(&block1[row + x], (n), &block2[row + x], (n)); \
  } \
  } \
  return sum; \
}

// Declare these functions to make sure the signature of the macro matches.
cost_pixel_nxn_func satd_8bit_4x4_generic;
cost_pixel_nxn_func satd_8bit_8x8_generic;
cost_pixel_nxn_func satd_8bit_16x16_generic;
cost_pixel_nxn_func satd_8bit_32x32_generic;
cost_pixel_nxn_func satd_8bit_64x64_generic;

// These macros define sadt_16bit_NxN for N = 8, 16, 32, 64
SATD_NXN(8, kvz_pixel, 8bit)
SATD_NXN(16, kvz_pixel, 8bit)
SATD_NXN(32, kvz_pixel, 8bit)
SATD_NXN(64, kvz_pixel, 8bit)

// Function macro for defining SAD calculating functions
// for fixed size blocks.
#define SAD_NXN(n, pixel_type, suffix) \
static unsigned sad_ ## suffix ## _ ##  n ## x ## n ## _generic( \
  const pixel_type * const block1, const pixel_type * const block2) \
{ \
  unsigned i; \
  unsigned sum = 0; \
  for (i = 0; i < (n)*(n); ++i) { \
  sum += abs(block1[i] - block2[i]); \
  } \
  return sum; \
}

// Declare these functions to make sure the signature of the macro matches.
static cost_pixel_nxn_func sad_8bit_4x4_generic;
static cost_pixel_nxn_func sad_8bit_8x8_generic;
static cost_pixel_nxn_func sad_8bit_16x16_generic;
static cost_pixel_nxn_func sad_8bit_32x32_generic;
static cost_pixel_nxn_func sad_8bit_64x64_generic;

// These macros define sad_16bit_nxn functions for n = 4, 8, 16, 32, 64
// with function signatures of cost_16bit_nxn_func.
// They are used through get_pixel_sad_func.
SAD_NXN(4, kvz_pixel, 8bit)
SAD_NXN(8, kvz_pixel, 8bit)
SAD_NXN(16, kvz_pixel, 8bit)
SAD_NXN(32, kvz_pixel, 8bit)
SAD_NXN(64, kvz_pixel, 8bit)


int strategy_register_picture_generic(void* opaque)
{
  bool success = true;

  success &= strategyselector_register(opaque, "reg_sad", "generic", 0, &reg_sad_generic);

  success &= strategyselector_register(opaque, "sad_8bit_4x4", "generic", 0, &sad_8bit_4x4_generic);
  success &= strategyselector_register(opaque, "sad_8bit_8x8", "generic", 0, &sad_8bit_8x8_generic);
  success &= strategyselector_register(opaque, "sad_8bit_16x16", "generic", 0, &sad_8bit_16x16_generic);
  success &= strategyselector_register(opaque, "sad_8bit_32x32", "generic", 0, &sad_8bit_32x32_generic);
  success &= strategyselector_register(opaque, "sad_8bit_64x64", "generic", 0, &sad_8bit_64x64_generic);

  success &= strategyselector_register(opaque, "satd_8bit_4x4", "generic", 0, &satd_8bit_4x4_generic);
  success &= strategyselector_register(opaque, "satd_8bit_8x8", "generic", 0, &satd_8bit_8x8_generic);
  success &= strategyselector_register(opaque, "satd_8bit_16x16", "generic", 0, &satd_8bit_16x16_generic);
  success &= strategyselector_register(opaque, "satd_8bit_32x32", "generic", 0, &satd_8bit_32x32_generic);
  success &= strategyselector_register(opaque, "satd_8bit_64x64", "generic", 0, &satd_8bit_64x64_generic);

  return success;
}
