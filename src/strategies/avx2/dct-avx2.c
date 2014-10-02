/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 */

#include <stdlib.h>

#include "dct-avx2.h"
#include "strategyselector.h"
#include "tables.h"

#if COMPILE_INTEL_AVX2
#include <immintrin.h>

extern const int16_t g_t4[4][4];
extern const int16_t g_t8[8][8];
extern const int16_t g_t16[16][16];
extern const int16_t g_t32[32][32];

static const int16_t dst_4x4[4][4] = 
{
  { 29, 55, 74, 84 },
  { 74, 74, 0, -74 },
  { 84, -29, -74, 55 },
  { 55, -84, 74, -29 }
};

static const int16_t dst_4x4_t[4][4] =
{
  { 29, 74, 84, 55 },
  { 55, 74, -29, -84 },
  { 74, 0, -74, 74 },
  { 84, -74, 55, -29 }
};

static const int16_t g_t4_t[4][4] =
{
  { 64, 83, 64, 36, },
  { 64, 36, -64, -83, },
  { 64, -36, -64, 83, },
  { 64, -83, 64, -36 }
};

static const int16_t g_t8_t[8][8] =
{
  { 64, 89, 83, 75, 64, 50, 36, 18, },
  { 64, 75, 36, -18, -64, -89, -83, -50, },
  { 64, 50, -36, -89, -64, 18, 83, 75, },
  { 64, 18, -83, -50, 64, 75, -36, -89, },
  { 64, -18, -83, 50, 64, -75, -36, 89, },
  { 64, -50, -36, 89, -64, -18, 83, -75, },
  { 64, -75, 36, 18, -64, 89, -83, 50, },
  { 64, -89, 83, -75, 64, -50, 36, -18 }
};

static const int16_t g_t16_t[16][16] =
{
  { 64, 90, 89, 87, 83, 80, 75, 70, 64, 57, 50, 43, 36, 25, 18, 9, },
  { 64, 87, 75, 57, 36, 9, -18, -43, -64, -80, -89, -90, -83, -70, -50, -25, },
  { 64, 80, 50, 9, -36, -70, -89, -87, -64, -25, 18, 57, 83, 90, 75, 43, },
  { 64, 70, 18, -43, -83, -87, -50, 9, 64, 90, 75, 25, -36, -80, -89, -57, },
  { 64, 57, -18, -80, -83, -25, 50, 90, 64, -9, -75, -87, -36, 43, 89, 70, },
  { 64, 43, -50, -90, -36, 57, 89, 25, -64, -87, -18, 70, 83, 9, -75, -80, },
  { 64, 25, -75, -70, 36, 90, 18, -80, -64, 43, 89, 9, -83, -57, 50, 87, },
  { 64, 9, -89, -25, 83, 43, -75, -57, 64, 70, -50, -80, 36, 87, -18, -90, },
  { 64, -9, -89, 25, 83, -43, -75, 57, 64, -70, -50, 80, 36, -87, -18, 90, },
  { 64, -25, -75, 70, 36, -90, 18, 80, -64, -43, 89, -9, -83, 57, 50, -87, },
  { 64, -43, -50, 90, -36, -57, 89, -25, -64, 87, -18, -70, 83, -9, -75, 80, },
  { 64, -57, -18, 80, -83, 25, 50, -90, 64, 9, -75, 87, -36, -43, 89, -70, },
  { 64, -70, 18, 43, -83, 87, -50, -9, 64, -90, 75, -25, -36, 80, -89, 57, },
  { 64, -80, 50, -9, -36, 70, -89, 87, -64, 25, 18, -57, 83, -90, 75, -43, },
  { 64, -87, 75, -57, 36, -9, -18, 43, -64, 80, -89, 90, -83, 70, -50, 25, },
  { 64, -90, 89, -87, 83, -80, 75, -70, 64, -57, 50, -43, 36, -25, 18, -9 }
};

static const int16_t g_t32_t[32][32] =
{
  { 64, 90, 90, 90, 89, 88, 87, 85, 83, 82, 80, 78, 75, 73, 70, 67, 64, 61, 57, 54, 50, 46, 43, 38, 36, 31, 25, 22, 18, 13, 9, 4, },
  { 64, 90, 87, 82, 75, 67, 57, 46, 36, 22, 9, -4, -18, -31, -43, -54, -64, -73, -80, -85, -89, -90, -90, -88, -83, -78, -70, -61, -50, -38, -25, -13, },
  { 64, 88, 80, 67, 50, 31, 9, -13, -36, -54, -70, -82, -89, -90, -87, -78, -64, -46, -25, -4, 18, 38, 57, 73, 83, 90, 90, 85, 75, 61, 43, 22, },
  { 64, 85, 70, 46, 18, -13, -43, -67, -83, -90, -87, -73, -50, -22, 9, 38, 64, 82, 90, 88, 75, 54, 25, -4, -36, -61, -80, -90, -89, -78, -57, -31, },
  { 64, 82, 57, 22, -18, -54, -80, -90, -83, -61, -25, 13, 50, 78, 90, 85, 64, 31, -9, -46, -75, -90, -87, -67, -36, 4, 43, 73, 89, 88, 70, 38, },
  { 64, 78, 43, -4, -50, -82, -90, -73, -36, 13, 57, 85, 89, 67, 25, -22, -64, -88, -87, -61, -18, 31, 70, 90, 83, 54, 9, -38, -75, -90, -80, -46, },
  { 64, 73, 25, -31, -75, -90, -70, -22, 36, 78, 90, 67, 18, -38, -80, -90, -64, -13, 43, 82, 89, 61, 9, -46, -83, -88, -57, -4, 50, 85, 87, 54, },
  { 64, 67, 9, -54, -89, -78, -25, 38, 83, 85, 43, -22, -75, -90, -57, 4, 64, 90, 70, 13, -50, -88, -80, -31, 36, 82, 87, 46, -18, -73, -90, -61, },
  { 64, 61, -9, -73, -89, -46, 25, 82, 83, 31, -43, -88, -75, -13, 57, 90, 64, -4, -70, -90, -50, 22, 80, 85, 36, -38, -87, -78, -18, 54, 90, 67, },
  { 64, 54, -25, -85, -75, -4, 70, 88, 36, -46, -90, -61, 18, 82, 80, 13, -64, -90, -43, 38, 89, 67, -9, -78, -83, -22, 57, 90, 50, -31, -87, -73, },
  { 64, 46, -43, -90, -50, 38, 90, 54, -36, -90, -57, 31, 89, 61, -25, -88, -64, 22, 87, 67, -18, -85, -70, 13, 83, 73, -9, -82, -75, 4, 80, 78, },
  { 64, 38, -57, -88, -18, 73, 80, -4, -83, -67, 25, 90, 50, -46, -90, -31, 64, 85, 9, -78, -75, 13, 87, 61, -36, -90, -43, 54, 89, 22, -70, -82, },
  { 64, 31, -70, -78, 18, 90, 43, -61, -83, 4, 87, 54, -50, -88, -9, 82, 64, -38, -90, -22, 75, 73, -25, -90, -36, 67, 80, -13, -89, -46, 57, 85, },
  { 64, 22, -80, -61, 50, 85, -9, -90, -36, 73, 70, -38, -89, -4, 87, 46, -64, -78, 25, 90, 18, -82, -57, 54, 83, -13, -90, -31, 75, 67, -43, -88, },
  { 64, 13, -87, -38, 75, 61, -57, -78, 36, 88, -9, -90, -18, 85, 43, -73, -64, 54, 80, -31, -89, 4, 90, 22, -83, -46, 70, 67, -50, -82, 25, 90, },
  { 64, 4, -90, -13, 89, 22, -87, -31, 83, 38, -80, -46, 75, 54, -70, -61, 64, 67, -57, -73, 50, 78, -43, -82, 36, 85, -25, -88, 18, 90, -9, -90, },
  { 64, -4, -90, 13, 89, -22, -87, 31, 83, -38, -80, 46, 75, -54, -70, 61, 64, -67, -57, 73, 50, -78, -43, 82, 36, -85, -25, 88, 18, -90, -9, 90, },
  { 64, -13, -87, 38, 75, -61, -57, 78, 36, -88, -9, 90, -18, -85, 43, 73, -64, -54, 80, 31, -89, -4, 90, -22, -83, 46, 70, -67, -50, 82, 25, -90, },
  { 64, -22, -80, 61, 50, -85, -9, 90, -36, -73, 70, 38, -89, 4, 87, -46, -64, 78, 25, -90, 18, 82, -57, -54, 83, 13, -90, 31, 75, -67, -43, 88, },
  { 64, -31, -70, 78, 18, -90, 43, 61, -83, -4, 87, -54, -50, 88, -9, -82, 64, 38, -90, 22, 75, -73, -25, 90, -36, -67, 80, 13, -89, 46, 57, -85, },
  { 64, -38, -57, 88, -18, -73, 80, 4, -83, 67, 25, -90, 50, 46, -90, 31, 64, -85, 9, 78, -75, -13, 87, -61, -36, 90, -43, -54, 89, -22, -70, 82, },
  { 64, -46, -43, 90, -50, -38, 90, -54, -36, 90, -57, -31, 89, -61, -25, 88, -64, -22, 87, -67, -18, 85, -70, -13, 83, -73, -9, 82, -75, -4, 80, -78, },
  { 64, -54, -25, 85, -75, 4, 70, -88, 36, 46, -90, 61, 18, -82, 80, -13, -64, 90, -43, -38, 89, -67, -9, 78, -83, 22, 57, -90, 50, 31, -87, 73, },
  { 64, -61, -9, 73, -89, 46, 25, -82, 83, -31, -43, 88, -75, 13, 57, -90, 64, 4, -70, 90, -50, -22, 80, -85, 36, 38, -87, 78, -18, -54, 90, -67, },
  { 64, -67, 9, 54, -89, 78, -25, -38, 83, -85, 43, 22, -75, 90, -57, -4, 64, -90, 70, -13, -50, 88, -80, 31, 36, -82, 87, -46, -18, 73, -90, 61, },
  { 64, -73, 25, 31, -75, 90, -70, 22, 36, -78, 90, -67, 18, 38, -80, 90, -64, 13, 43, -82, 89, -61, 9, 46, -83, 88, -57, 4, 50, -85, 87, -54, },
  { 64, -78, 43, 4, -50, 82, -90, 73, -36, -13, 57, -85, 89, -67, 25, 22, -64, 88, -87, 61, -18, -31, 70, -90, 83, -54, 9, 38, -75, 90, -80, 46, },
  { 64, -82, 57, -22, -18, 54, -80, 90, -83, 61, -25, -13, 50, -78, 90, -85, 64, -31, -9, 46, -75, 90, -87, 67, -36, -4, 43, -73, 89, -88, 70, -38, },
  { 64, -85, 70, -46, 18, 13, -43, 67, -83, 90, -87, 73, -50, 22, 9, -38, 64, -82, 90, -88, 75, -54, 25, 4, -36, 61, -80, 90, -89, 78, -57, 31, },
  { 64, -88, 80, -67, 50, -31, 9, 13, -36, 54, -70, 82, -89, 90, -87, 78, -64, 46, -25, 4, 18, -38, 57, -73, 83, -90, 90, -85, 75, -61, 43, -22, },
  { 64, -90, 87, -82, 75, -67, 57, -46, 36, -22, 9, 4, -18, 31, -43, 54, -64, 73, -80, 85, -89, 90, -90, 88, -83, 78, -70, 61, -50, 38, -25, 13, },
  { 64, -90, 90, -90, 89, -88, 87, -85, 83, -82, 80, -78, 75, -73, 70, -67, 64, -61, 57, -54, 50, -46, 43, -38, 36, -31, 25, -22, 18, -13, 9, -4 }
};

/**
* \brief AVX2 transform functions
*
* TODO: description
*
* \param TODO
*
* \returns TODO
*/
static void transpose_4x4_16bit(const int16_t *src, int16_t *dst)
{
  __m256i original;
  __m128i upper, lower, tmp0, tmp1;
  original = _mm256_loadu_si256((__m256i*)src);
  upper = _mm256_castsi256_si128(original);
  lower = _mm256_extracti128_si256(original, 1);

  tmp0 = _mm_unpacklo_epi16(upper, lower);
  tmp1 = _mm_unpackhi_epi16(upper, lower);

  upper = _mm_unpacklo_epi16(tmp0, tmp1);
  lower = _mm_unpackhi_epi16(tmp0, tmp1);

  _mm256_storeu_si256((__m256i*)dst, _mm256_inserti128_si256(_mm256_castsi128_si256(upper), lower, 1));
}


static void transpose_8x8_16bit(const int16_t *src, int16_t *dst)
{
  __m256i row_pair[4], tmp[4];
  row_pair[0] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)src + 0)), _mm_loadu_si128((__m128i*)src + 4), 1);
  row_pair[1] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)src + 1)), _mm_loadu_si128((__m128i*)src + 5), 1);
  row_pair[2] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)src + 2)), _mm_loadu_si128((__m128i*)src + 6), 1);
  row_pair[3] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)src + 3)), _mm_loadu_si128((__m128i*)src + 7), 1);

  tmp[0] = _mm256_unpacklo_epi16(row_pair[0], row_pair[1]);
  tmp[1] = _mm256_unpackhi_epi16(row_pair[0], row_pair[1]);
  tmp[2] = _mm256_unpacklo_epi16(row_pair[2], row_pair[3]);
  tmp[3] = _mm256_unpackhi_epi16(row_pair[2], row_pair[3]);

  row_pair[0] = _mm256_unpacklo_epi32(tmp[0], tmp[2]);
  row_pair[1] = _mm256_unpackhi_epi32(tmp[0], tmp[2]);
  row_pair[2] = _mm256_unpacklo_epi32(tmp[1], tmp[3]);
  row_pair[3] = _mm256_unpackhi_epi32(tmp[1], tmp[3]);

  tmp[0] = _mm256_unpacklo_epi64(row_pair[0], row_pair[2]);
  tmp[1] = _mm256_unpackhi_epi64(row_pair[0], row_pair[2]);
  tmp[2] = _mm256_unpacklo_epi64(row_pair[1], row_pair[3]);
  tmp[3] = _mm256_unpackhi_epi64(row_pair[1], row_pair[3]);

  tmp[0] = _mm256_permute4x64_epi64(tmp[0], 0 + 8 + 16 + 192);
  tmp[1] = _mm256_permute4x64_epi64(tmp[1], 0 + 8 + 16 + 192);
  tmp[2] = _mm256_permute4x64_epi64(tmp[2], 0 + 8 + 16 + 192);
  tmp[3] = _mm256_permute4x64_epi64(tmp[3], 0 + 8 + 16 + 192);

  _mm_storeu_si128((__m128i*)dst + 0, _mm256_castsi256_si128(tmp[0]));
  _mm_storeu_si128((__m128i*)dst + 1, _mm256_castsi256_si128(tmp[1]));
  _mm_storeu_si128((__m128i*)dst + 2, _mm256_castsi256_si128(tmp[2]));
  _mm_storeu_si128((__m128i*)dst + 3, _mm256_castsi256_si128(tmp[3]));

  _mm_storeu_si128((__m128i*)dst + 4, _mm256_extracti128_si256(tmp[0], 1));
  _mm_storeu_si128((__m128i*)dst + 5, _mm256_extracti128_si256(tmp[1], 1));
  _mm_storeu_si128((__m128i*)dst + 6, _mm256_extracti128_si256(tmp[2], 1));
  _mm_storeu_si128((__m128i*)dst + 7, _mm256_extracti128_si256(tmp[3], 1));
}

static void transpose_16x16_16bit(const int16_t *src, int16_t *dst)
{
  int i;
  __m256i row[16], tmp[16];
  for (i = 0; i < 16; ++i) {
    row[i] = _mm256_loadu_si256((__m256i*) src + i);
  }

  for (i = 0; i < 16; i += 4) {
    tmp[i + 0] = _mm256_unpacklo_epi16(row[i + 0], row[i + 1]);
    tmp[i + 1] = _mm256_unpackhi_epi16(row[i + 0], row[i + 1]);
    tmp[i + 2] = _mm256_unpacklo_epi16(row[i + 2], row[i + 3]);
    tmp[i + 3] = _mm256_unpackhi_epi16(row[i + 2], row[i + 3]);
  }
  for (i = 0; i < 16; i += 4) {
    row[i + 0] = _mm256_unpacklo_epi32(tmp[i + 0], tmp[i + 2]);
    row[i + 1] = _mm256_unpackhi_epi32(tmp[i + 0], tmp[i + 2]);
    row[i + 2] = _mm256_unpacklo_epi32(tmp[i + 1], tmp[i + 3]);
    row[i + 3] = _mm256_unpackhi_epi32(tmp[i + 1], tmp[i + 3]);
  }

  for (i = 0; i < 8; i += 2) {
    tmp[i + 0] = _mm256_unpacklo_epi64(row[i / 2 + 0], row[i / 2 + 4]);
    tmp[i + 1] = _mm256_unpackhi_epi64(row[i / 2 + 0], row[i / 2 + 4]);
  }
  for (i = 8; i < 16; i += 2) {
    tmp[i + 0] = _mm256_unpacklo_epi64(row[i / 2 + 4], row[i / 2 + 8]);
    tmp[i + 1] = _mm256_unpackhi_epi64(row[i / 2 + 4], row[i / 2 + 8]);
  }


  for (i = 0; i < 8; ++i) {
    _mm_storeu_si128((__m128i*)dst + 2 * i, _mm256_extracti128_si256(tmp[i], 0));
    _mm_storeu_si128(((__m128i*)dst) + 2 * i + 1, _mm256_extracti128_si256(tmp[i + 8], 0));
  }
  for (i = 8; i < 16; ++i) {
    _mm_storeu_si128((__m128i*)dst + 2 * i, _mm256_extracti128_si256(tmp[(i + 8) % 16], 1));
    _mm_storeu_si128((__m128i*)dst + 2 * i + 1, _mm256_extracti128_si256(tmp[i], 1));
  }
}


static void transpose_32x32_16bit(const int16_t *src, int16_t *dst)
{
  int i;
  __m256i row[32][2], tmp[32][2];
  for (i = 0; i < 32; ++i) {
    row[i][0] = _mm256_loadu_si256((__m256i*) src + 2 * i);
    row[i][1] = _mm256_loadu_si256((__m256i*) src + 2 * i + 1);
  }

  for (i = 0; i < 32; i += 4) {
    tmp[i + 0][0] = _mm256_unpacklo_epi16(row[i + 0][0], row[i + 1][0]);
    tmp[i + 1][0] = _mm256_unpackhi_epi16(row[i + 0][0], row[i + 1][0]);
    tmp[i + 2][0] = _mm256_unpacklo_epi16(row[i + 2][0], row[i + 3][0]);
    tmp[i + 3][0] = _mm256_unpackhi_epi16(row[i + 2][0], row[i + 3][0]);

    tmp[i + 0][1] = _mm256_unpacklo_epi16(row[i + 0][1], row[i + 1][1]);
    tmp[i + 1][1] = _mm256_unpackhi_epi16(row[i + 0][1], row[i + 1][1]);
    tmp[i + 2][1] = _mm256_unpacklo_epi16(row[i + 2][1], row[i + 3][1]);
    tmp[i + 3][1] = _mm256_unpackhi_epi16(row[i + 2][1], row[i + 3][1]);
  }
  for (i = 0; i < 32; i += 4) {
    row[i + 0][0] = _mm256_unpacklo_epi32(tmp[i + 0][0], tmp[i + 2][0]);
    row[i + 1][0] = _mm256_unpackhi_epi32(tmp[i + 0][0], tmp[i + 2][0]);
    row[i + 2][0] = _mm256_unpacklo_epi32(tmp[i + 1][0], tmp[i + 3][0]);
    row[i + 3][0] = _mm256_unpackhi_epi32(tmp[i + 1][0], tmp[i + 3][0]);

    row[i + 0][1] = _mm256_unpacklo_epi32(tmp[i + 0][1], tmp[i + 2][1]);
    row[i + 1][1] = _mm256_unpackhi_epi32(tmp[i + 0][1], tmp[i + 2][1]);
    row[i + 2][1] = _mm256_unpacklo_epi32(tmp[i + 1][1], tmp[i + 3][1]);
    row[i + 3][1] = _mm256_unpackhi_epi32(tmp[i + 1][1], tmp[i + 3][1]);
  }

  for (i = 0; i < 8; i += 2) {
    tmp[i + 0][0] = _mm256_unpacklo_epi64(row[i / 2 + 0][0], row[i / 2 + 4][0]);
    tmp[i + 1][0] = _mm256_unpackhi_epi64(row[i / 2 + 0][0], row[i / 2 + 4][0]);

    tmp[i + 0][1] = _mm256_unpacklo_epi64(row[i / 2 + 0][1], row[i / 2 + 4][1]);
    tmp[i + 1][1] = _mm256_unpackhi_epi64(row[i / 2 + 0][1], row[i / 2 + 4][1]);
  }

  for (i = 8; i < 16; i += 2) {
    tmp[i + 0][0] = _mm256_unpacklo_epi64(row[i / 2 + 4][0], row[i / 2 + 8][0]);
    tmp[i + 1][0] = _mm256_unpackhi_epi64(row[i / 2 + 4][0], row[i / 2 + 8][0]);

    tmp[i + 0][1] = _mm256_unpacklo_epi64(row[i / 2 + 4][1], row[i / 2 + 8][1]);
    tmp[i + 1][1] = _mm256_unpackhi_epi64(row[i / 2 + 4][1], row[i / 2 + 8][1]);
  }

  for (i = 16; i < 24; i += 2) {
    tmp[i + 0][0] = _mm256_unpacklo_epi64(row[i / 2 + 8][0], row[i / 2 + 12][0]);
    tmp[i + 1][0] = _mm256_unpackhi_epi64(row[i / 2 + 8][0], row[i / 2 + 12][0]);

    tmp[i + 0][1] = _mm256_unpacklo_epi64(row[i / 2 + 8][1], row[i / 2 + 12][1]);
    tmp[i + 1][1] = _mm256_unpackhi_epi64(row[i / 2 + 8][1], row[i / 2 + 12][1]);
  }

  for (i = 24; i < 32; i += 2) {
    tmp[i + 0][0] = _mm256_unpacklo_epi64(row[i / 2 + 12][0], row[i / 2 + 16][0]);
    tmp[i + 1][0] = _mm256_unpackhi_epi64(row[i / 2 + 12][0], row[i / 2 + 16][0]);

    tmp[i + 0][1] = _mm256_unpacklo_epi64(row[i / 2 + 12][1], row[i / 2 + 16][1]);
    tmp[i + 1][1] = _mm256_unpackhi_epi64(row[i / 2 + 12][1], row[i / 2 + 16][1]);
  }


  for (i = 0; i < 8; ++i) {
    _mm_storeu_si128((__m128i*)dst + 4 * i, _mm256_extracti128_si256(tmp[i][0], 0));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 1, _mm256_extracti128_si256(tmp[i + 8][0], 0));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 2, _mm256_extracti128_si256(tmp[i + 16][0], 0));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 3, _mm256_extracti128_si256(tmp[i + 24][0], 0));
  }
  for (i = 8; i < 16; ++i) {
    _mm_storeu_si128((__m128i*)dst + 4 * i, _mm256_extracti128_si256(tmp[i - 8][0], 1));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 1, _mm256_extracti128_si256(tmp[i + 8 - 8][0], 1));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 2, _mm256_extracti128_si256(tmp[i + 16 - 8][0], 1));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 3, _mm256_extracti128_si256(tmp[i + 24 - 8][0], 1));
  }
  for (i = 16; i < 24; ++i) {
    _mm_storeu_si128((__m128i*)dst + 4 * i, _mm256_extracti128_si256(tmp[i - 16][1], 0));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 1, _mm256_extracti128_si256(tmp[i + 8 - 16][1], 0));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 2, _mm256_extracti128_si256(tmp[i + 16 - 16][1], 0));
    _mm_storeu_si128(((__m128i*)dst) + 4 * i + 3, _mm256_extracti128_si256(tmp[(i + 24 - 16)][1], 0));
  }
  for (i = 24; i < 32; ++i) {
    _mm_storeu_si128((__m128i*)dst + 4 * i, _mm256_extracti128_si256(tmp[(i - 24) % 32][1], 1));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 1, _mm256_extracti128_si256(tmp[i + 8 - 24][1], 1));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 2, _mm256_extracti128_si256(tmp[(i + 16 - 24) % 32][1], 1));
    _mm_storeu_si128((__m128i*)dst + 4 * i + 3, _mm256_extracti128_si256(tmp[i + 24 - 24][1], 1));
  }
}


static void mul_clip_matrix_4x4_avx2(const int16_t *first, const int16_t *second, int16_t *dst, int32_t shift)
{
  __m256i b[2], a, result,  even[2], odd[2];

  const int32_t add = 1 << (shift - 1);

  a = _mm256_loadu_si256((__m256i*) first);
  b[0] = _mm256_loadu_si256((__m256i*) second);

  b[0] = _mm256_unpacklo_epi16(b[0], _mm256_srli_si256(b[0], 8));
  b[1] = _mm256_permute2x128_si256(b[0], b[0], 1 + 16);
  b[0] = _mm256_permute2x128_si256(b[0], b[0], 0);

  even[0] = _mm256_shuffle_epi32(a, 0);
  odd[0] = _mm256_shuffle_epi32(a, 1 + 4 + 16 + 64);

  even[0] = _mm256_madd_epi16(even[0], b[0]);
  odd[0] = _mm256_madd_epi16(odd[0], b[1]);

  result = _mm256_add_epi32(even[0], odd[0]);
  result = _mm256_add_epi32(result, _mm256_set1_epi32(add));
  result = _mm256_srai_epi32(result, shift);

  even[1] = _mm256_shuffle_epi32(a, 2 + 8 + 32 + 128);
  odd[1] = _mm256_shuffle_epi32(a, 3 + 12 + 48 + 192);

  even[1] = _mm256_madd_epi16(even[1], b[0]);
  odd[1] = _mm256_madd_epi16(odd[1], b[1]);

  odd[1] = _mm256_add_epi32(even[1], odd[1]);
  odd[1] = _mm256_add_epi32(odd[1], _mm256_set1_epi32(add));
  odd[1] = _mm256_srai_epi32(odd[1], shift);

  result = _mm256_packs_epi32(result, odd[1]);

  _mm256_storeu_si256((__m256i*)dst, result);
}

static void mul_clip_matrix_8x8_avx2(const int16_t *first, const int16_t *second, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i b[2], accu[8], even[2], odd[2];

  const int32_t add = 1 << (shift - 1);

  b[0] = _mm256_loadu_si256((__m256i*) second);

  b[1] = _mm256_unpackhi_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
  b[0] = _mm256_unpacklo_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
  b[0] = _mm256_inserti128_si256(b[0], _mm256_castsi256_si128(b[1]), 1);

  for (i = 0; i < 8; i += 2) {

    even[0] = _mm256_set1_epi32(((int32_t*)first)[4 * i]);
    even[0] = _mm256_madd_epi16(even[0], b[0]);
    accu[i] = even[0];
    
    odd[0] = _mm256_set1_epi32(((int32_t*)first)[4 * (i + 1)]);
    odd[0] = _mm256_madd_epi16(odd[0], b[0]);
    accu[i + 1] = odd[0];
  }

  for (j = 1; j < 4; ++j) {

    b[0] = _mm256_loadu_si256((__m256i*)second + j);

    b[1] = _mm256_unpackhi_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
    b[0] = _mm256_unpacklo_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
    b[0] = _mm256_inserti128_si256(b[0], _mm256_castsi256_si128(b[1]), 1);

    for (i = 0; i < 8; i += 2) {
    
      even[0] = _mm256_set1_epi32(((int32_t*)first)[4 * i + j]);
      even[0] = _mm256_madd_epi16(even[0], b[0]);
      accu[i] = _mm256_add_epi32(accu[i], even[0]);

      odd[0] = _mm256_set1_epi32(((int32_t*)first)[4 * (i + 1) + j]);
      odd[0] = _mm256_madd_epi16(odd[0], b[0]);
      accu[i + 1] = _mm256_add_epi32(accu[i + 1], odd[0]);
    }
  }

  for (i = 0; i < 8; i += 2) {
     __m256i result, first_half, second_half;
  
    first_half = _mm256_srai_epi32(_mm256_add_epi32(accu[i], _mm256_set1_epi32(add)), shift);
    second_half = _mm256_srai_epi32(_mm256_add_epi32(accu[i + 1], _mm256_set1_epi32(add)), shift);
    result = _mm256_permute4x64_epi64(_mm256_packs_epi32(first_half, second_half), 0 + 8 + 16 + 192);
    _mm256_storeu_si256((__m256i*)dst + i / 2, result);

  }
}

static void mul_clip_matrix_16x16_avx2(const int16_t *first, const int16_t *second, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i row[4], accu[16][2], even, odd;

  const int32_t stride = 8;

  const int32_t add = 1 << (shift - 1);

  row[0] = _mm256_loadu_si256((__m256i*) second);
  row[1] = _mm256_loadu_si256((__m256i*) second + 1);
  row[2] = _mm256_unpacklo_epi16(row[0], row[1]);
  row[3] = _mm256_unpackhi_epi16(row[0], row[1]);
  row[0] = _mm256_permute2x128_si256(row[2], row[3], 0 + 32);
  row[1] = _mm256_permute2x128_si256(row[2], row[3], 1 + 48);

  for (i = 0; i < 16; i += 2) {

    even = _mm256_set1_epi32(((int32_t*)first)[stride * i]);
    accu[i][0] = _mm256_madd_epi16(even, row[0]);
    accu[i][1] = _mm256_madd_epi16(even, row[1]);

    odd = _mm256_set1_epi32(((int32_t*)first)[stride * (i + 1)]);
    accu[i+1][0] = _mm256_madd_epi16(odd, row[0]);
    accu[i+1][1] = _mm256_madd_epi16(odd, row[1]);
  }

  for (j = 2; j < 16; j+=2) {

    row[0] = _mm256_loadu_si256((__m256i*)second + j);
    row[1] = _mm256_loadu_si256((__m256i*)second + j + 1);
    row[2] = _mm256_unpacklo_epi16(row[0], row[1]);
    row[3] = _mm256_unpackhi_epi16(row[0], row[1]);
    row[0] = _mm256_permute2x128_si256(row[2], row[3], 0 + 32);
    row[1] = _mm256_permute2x128_si256(row[2], row[3], 1 + 48);

    for (i = 0; i < 16; i += 2) {

      even = _mm256_set1_epi32(((int32_t*)first)[stride * i + j/2]);
      accu[i][0] = _mm256_add_epi32(accu[i][0], _mm256_madd_epi16(even, row[0]));
      accu[i][1] = _mm256_add_epi32(accu[i][1], _mm256_madd_epi16(even, row[1]));

      odd = _mm256_set1_epi32(((int32_t*)first)[stride * (i + 1) + j/2]);
      accu[i + 1][0] = _mm256_add_epi32(accu[i + 1][0], _mm256_madd_epi16(odd, row[0]));
      accu[i + 1][1] = _mm256_add_epi32(accu[i + 1][1], _mm256_madd_epi16(odd, row[1]));

    }
  }

  for (i = 0; i < 16; ++i) {
    __m256i result, first_half, second_half;

    first_half = _mm256_srai_epi32(_mm256_add_epi32(accu[i][0], _mm256_set1_epi32(add)), shift);
    second_half = _mm256_srai_epi32(_mm256_add_epi32(accu[i][1], _mm256_set1_epi32(add)), shift);
    result = _mm256_permute4x64_epi64(_mm256_packs_epi32(first_half, second_half), 0 + 8 + 16 + 192);
    _mm256_storeu_si256((__m256i*)dst + i, result);

  }
}

static void mul_clip_matrix_32x32_avx2(const int16_t *first, const int16_t *second, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i row[4], tmp[2], accu[32][4], even, odd;

  const int32_t stride = 16;

  const int32_t add = 1 << (shift - 1);

  row[0] = _mm256_loadu_si256((__m256i*) second);
  row[1] = _mm256_loadu_si256((__m256i*) second + 2);
  tmp[0] = _mm256_unpacklo_epi16(row[0], row[1]);
  tmp[1] = _mm256_unpackhi_epi16(row[0], row[1]);
  row[0] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
  row[1] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

  row[2] = _mm256_loadu_si256((__m256i*) second + 1);
  row[3] = _mm256_loadu_si256((__m256i*) second + 3);
  tmp[0] = _mm256_unpacklo_epi16(row[2], row[3]);
  tmp[1] = _mm256_unpackhi_epi16(row[2], row[3]);
  row[2] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
  row[3] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

  for (i = 0; i < 32; i += 2) {

    even = _mm256_set1_epi32(((int32_t*)first)[stride * i]);
    accu[i][0] = _mm256_madd_epi16(even, row[0]);
    accu[i][1] = _mm256_madd_epi16(even, row[1]);
    accu[i][2] = _mm256_madd_epi16(even, row[2]);
    accu[i][3] = _mm256_madd_epi16(even, row[3]);

    odd = _mm256_set1_epi32(((int32_t*)first)[stride * (i + 1)]);
    accu[i + 1][0] = _mm256_madd_epi16(odd, row[0]);
    accu[i + 1][1] = _mm256_madd_epi16(odd, row[1]);
    accu[i + 1][2] = _mm256_madd_epi16(odd, row[2]);
    accu[i + 1][3] = _mm256_madd_epi16(odd, row[3]);
  }

  for (j = 4; j < 64; j += 4) {

    row[0] = _mm256_loadu_si256((__m256i*)second + j);
    row[1] = _mm256_loadu_si256((__m256i*)second + j + 2);
    tmp[0] = _mm256_unpacklo_epi16(row[0], row[1]);
    tmp[1] = _mm256_unpackhi_epi16(row[0], row[1]);
    row[0] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
    row[1] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

    row[2] = _mm256_loadu_si256((__m256i*) second + j + 1);
    row[3] = _mm256_loadu_si256((__m256i*) second + j + 3);
    tmp[0] = _mm256_unpacklo_epi16(row[2], row[3]);
    tmp[1] = _mm256_unpackhi_epi16(row[2], row[3]);
    row[2] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
    row[3] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

    for (i = 0; i < 32; i += 2) {

      even = _mm256_set1_epi32(((int32_t*)first)[stride * i + j / 4]);
      accu[i][0] = _mm256_add_epi32(accu[i][0], _mm256_madd_epi16(even, row[0]));
      accu[i][1] = _mm256_add_epi32(accu[i][1], _mm256_madd_epi16(even, row[1]));
      accu[i][2] = _mm256_add_epi32(accu[i][2], _mm256_madd_epi16(even, row[2]));
      accu[i][3] = _mm256_add_epi32(accu[i][3], _mm256_madd_epi16(even, row[3]));

      odd = _mm256_set1_epi32(((int32_t*)first)[stride * (i + 1) + j / 4]);
      accu[i + 1][0] = _mm256_add_epi32(accu[i + 1][0], _mm256_madd_epi16(odd, row[0]));
      accu[i + 1][1] = _mm256_add_epi32(accu[i + 1][1], _mm256_madd_epi16(odd, row[1]));
      accu[i + 1][2] = _mm256_add_epi32(accu[i + 1][2], _mm256_madd_epi16(odd, row[2]));
      accu[i + 1][3] = _mm256_add_epi32(accu[i + 1][3], _mm256_madd_epi16(odd, row[3]));

    }
  }

  for (i = 0; i < 32; ++i) {
    __m256i result, first_quarter, second_quarter, third_quarter, fourth_quarter;

    first_quarter = _mm256_srai_epi32(_mm256_add_epi32(accu[i][0], _mm256_set1_epi32(add)), shift);
    second_quarter = _mm256_srai_epi32(_mm256_add_epi32(accu[i][1], _mm256_set1_epi32(add)), shift);
    third_quarter = _mm256_srai_epi32(_mm256_add_epi32(accu[i][2], _mm256_set1_epi32(add)), shift);
    fourth_quarter = _mm256_srai_epi32(_mm256_add_epi32(accu[i][3], _mm256_set1_epi32(add)), shift);
    result = _mm256_permute4x64_epi64(_mm256_packs_epi32(first_quarter, second_quarter), 0 + 8 + 16 + 192);
    _mm256_storeu_si256((__m256i*)dst + 2 * i, result);
    result = _mm256_permute4x64_epi64(_mm256_packs_epi32(third_quarter, fourth_quarter), 0 + 8 + 16 + 192);
    _mm256_storeu_si256((__m256i*)dst + 2 * i + 1, result);

  }
}

static void matrix_dst_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4 * 4];

  mul_clip_matrix_4x4_avx2(src, (int16_t*)dst_4x4_t, tmp, shift0);
  mul_clip_matrix_4x4_avx2((int16_t*)dst_4x4, tmp, dst, shift1);
}

static void matrix_idst_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4 * 4];

  mul_clip_matrix_4x4_avx2((int16_t*)dst_4x4_t, src, tmp, shift0);
  mul_clip_matrix_4x4_avx2(tmp, (int16_t*)dst_4x4, dst, shift1);
}

static void matrix_transform_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4 * 4];

  mul_clip_matrix_4x4_avx2(src, (int16_t*)g_t4_t, tmp, shift0);
  mul_clip_matrix_4x4_avx2((int16_t*)g_t4, tmp, dst, shift1);
}

static void matrix_itransform_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4*4];

  mul_clip_matrix_4x4_avx2((int16_t*)g_t4_t, src, tmp, shift0);
  mul_clip_matrix_4x4_avx2(tmp, (int16_t*)g_t4, dst, shift1);
}

static void matrix_transform_2d_8x8_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[8 * 8];

  mul_clip_matrix_8x8_avx2(src, (int16_t*)g_t8_t, tmp, shift0);
  mul_clip_matrix_8x8_avx2((int16_t*)g_t8, tmp, dst, shift1);
}

static void matrix_itransform_2d_8x8_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[8 * 8];

  mul_clip_matrix_8x8_avx2((int16_t*)g_t8_t, src, tmp, shift0);
  mul_clip_matrix_8x8_avx2(tmp, (int16_t*)g_t8, dst, shift1);
}

static void matrix_transform_2d_16x16_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[16 * 16];

  mul_clip_matrix_16x16_avx2(src, (int16_t*)g_t16_t, tmp, shift0);
  mul_clip_matrix_16x16_avx2((int16_t*)g_t16, tmp, dst, shift1);
}

static void matrix_itransform_2d_16x16_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[16 * 16];

  mul_clip_matrix_16x16_avx2((int16_t*)g_t16_t, src, tmp, shift0);
  mul_clip_matrix_16x16_avx2(tmp, (int16_t*)g_t16, dst, shift1);
}

static void matrix_transform_2d_32x32_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[32 * 32];

  mul_clip_matrix_32x32_avx2(src, (int16_t*)g_t32_t, tmp, shift0);
  mul_clip_matrix_32x32_avx2((int16_t*)g_t32, tmp, dst, shift1);
}

static void matrix_itransform_2d_32x32_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[32 * 32];

  mul_clip_matrix_32x32_avx2((int16_t*) g_t32_t, src, tmp, shift0);
  mul_clip_matrix_32x32_avx2(tmp, (int16_t*)g_t32, dst, shift1);
}

static void matrix_dst_4x4_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[4] + 8;
  matrix_dst_2d_4x4_avx2(src, dst, (const int16_t*)dst_4x4, shift_1st, shift_2nd);
}

static void matrix_idst_4x4_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_idst_2d_4x4_avx2(src, dst, (const int16_t*)dst_4x4, shift_1st, shift_2nd);
}

static void matrix_dct_4x4_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[4] + 8;
  matrix_transform_2d_4x4_avx2(src, dst, (const int16_t*)g_t4, shift_1st, shift_2nd);
}

static void matrix_idct_4x4_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_itransform_2d_4x4_avx2(src, dst, (const int16_t*)g_t4, shift_1st, shift_2nd);
}

static void matrix_dct_8x8_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[8] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[8] + 8;
  matrix_transform_2d_8x8_avx2(src, dst, (const int16_t*)g_t8, shift_1st, shift_2nd);
}

static void matrix_idct_8x8_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_itransform_2d_8x8_avx2(src, dst, (const int16_t*)g_t8_t, shift_1st, shift_2nd);
}

static void matrix_dct_16x16_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[16] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[16] + 8;
  matrix_transform_2d_16x16_avx2(src, dst, (const int16_t*)g_t16, shift_1st, shift_2nd);
}

static void matrix_idct_16x16_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_itransform_2d_16x16_avx2(src, dst, (const int16_t*)g_t16, shift_1st, shift_2nd);
}

static void matrix_dct_32x32_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[32] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[32] + 8;
  matrix_transform_2d_32x32_avx2(src, dst, (const int16_t*)g_t32, shift_1st, shift_2nd);
}

static void matrix_idct_32x32_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_itransform_2d_32x32_avx2(src, dst, (const int16_t*)g_t32, shift_1st, shift_2nd);
}
#endif //COMPILE_INTEL_AVX2

int strategy_register_dct_avx2(void* opaque)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  success &= strategyselector_register(opaque, "fast_forward_dst_4x4", "avx2", 40, &matrix_dst_4x4_avx2);

  success &= strategyselector_register(opaque, "dct_4x4", "avx2", 40, &matrix_dct_4x4_avx2);
  success &= strategyselector_register(opaque, "dct_8x8", "avx2", 40, &matrix_dct_8x8_avx2);
  success &= strategyselector_register(opaque, "dct_16x16", "avx2", 40, &matrix_dct_16x16_avx2);
  success &= strategyselector_register(opaque, "dct_32x32", "avx2", 40, &matrix_dct_32x32_avx2);

  success &= strategyselector_register(opaque, "fast_inverse_dst_4x4", "avx2", 40, &matrix_idst_4x4_avx2);

  success &= strategyselector_register(opaque, "idct_4x4", "avx2", 40, &matrix_idct_4x4_avx2);
  success &= strategyselector_register(opaque, "idct_8x8", "avx2", 40, &matrix_idct_8x8_avx2);
  success &= strategyselector_register(opaque, "idct_16x16", "avx2", 40, &matrix_idct_16x16_avx2);
  success &= strategyselector_register(opaque, "idct_32x32", "avx2", 40, &matrix_idct_32x32_avx2);
#endif //COMPILE_INTEL_AVX2  
  return success;
}
