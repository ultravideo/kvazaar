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


static void mul_matrix_4x4_avx2(const int16_t *first, const int16_t *second, int16_t *dst, int32_t shift)
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

static void mul_matrix_8x8_avx2(const int16_t *first, const int16_t *second, int16_t *dst, const int32_t shift)
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

static void mul_matrix_16x16_avx2(const int16_t *first, const int16_t *second, int16_t *dst, const int32_t shift)
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

static void matrix_transform_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4 * 4];
  int16_t transposed[16];

  transpose_4x4_16bit(transform, transposed);
  mul_matrix_4x4_avx2(src, transposed, tmp, shift0);
  mul_matrix_4x4_avx2(transform, tmp, dst, shift1);
}

static void matrix_itransform_2d_4x4_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[4*4];
  int16_t transposed[16];

  transpose_4x4_16bit(transform, transposed);
  mul_matrix_4x4_avx2(transposed, src, tmp, shift0);
  mul_matrix_4x4_avx2(tmp, transform, dst, shift1);
}

static void matrix_transform_2d_8x8_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[8 * 8];
  int16_t transposed[64];

  transpose_8x8_16bit(transform, transposed);
  mul_matrix_8x8_avx2(src, transposed, tmp, shift0);
  mul_matrix_8x8_avx2(transform, tmp, dst, shift1);
}

static void matrix_itransform_2d_8x8_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[8 * 8];
  int16_t transposed[64];

  transpose_8x8_16bit(transform, transposed);
  mul_matrix_8x8_avx2(transposed, src, tmp, shift0);
  mul_matrix_8x8_avx2(tmp, transform, dst, shift1);
}

static void matrix_transform_2d_16x16_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[16 * 16];
  int16_t transposed[16 * 16];

  transpose_16x16_16bit(transform, transposed);
  mul_matrix_16x16_avx2(src, transposed, tmp, shift0);
  mul_matrix_16x16_avx2(transform, tmp, dst, shift1);
}

static void matrix_itransform_2d_16x16_avx2(const int16_t *src, int16_t *dst, const int16_t *transform, const int16_t shift0, const int16_t shift1)
{
  int16_t tmp[16 * 16];
  int16_t transposed[16 * 16];

  transpose_16x16_16bit(transform, transposed);
  mul_matrix_16x16_avx2(transposed, src, tmp, shift0);
  mul_matrix_16x16_avx2(tmp, transform, dst, shift1);
}

static void partial_butterfly_32_avx2(short *src, short *dst,
  int32_t shift)
{
  int32_t j, k;
  int32_t e[16], o[16];
  int32_t ee[8], eo[8];
  int32_t eee[4], eeo[4];
  int32_t eeee[2], eeeo[2];
  int32_t add = 1 << (shift - 1);
  const int32_t line = 32;

  for (j = 0; j < line; j++) {
    // E and O
    for (k = 0; k < 16; k++) {
      e[k] = src[k] + src[31 - k];
      o[k] = src[k] - src[31 - k];
    }
    // EE and EO
    for (k = 0; k < 8; k++) {
      ee[k] = e[k] + e[15 - k];
      eo[k] = e[k] - e[15 - k];
    }
    // EEE and EEO
    for (k = 0; k < 4; k++) {
      eee[k] = ee[k] + ee[7 - k];
      eeo[k] = ee[k] - ee[7 - k];
    }
    // EEEE and EEEO
    eeee[0] = eee[0] + eee[3];
    eeeo[0] = eee[0] - eee[3];
    eeee[1] = eee[1] + eee[2];
    eeeo[1] = eee[1] - eee[2];

    dst[0] = (short)((g_t32[0][0] * eeee[0] + g_t32[0][1] * eeee[1] + add) >> shift);
    dst[16 * line] = (short)((g_t32[16][0] * eeee[0] + g_t32[16][1] * eeee[1] + add) >> shift);
    dst[8 * line] = (short)((g_t32[8][0] * eeeo[0] + g_t32[8][1] * eeeo[1] + add) >> shift);
    dst[24 * line] = (short)((g_t32[24][0] * eeeo[0] + g_t32[24][1] * eeeo[1] + add) >> shift);
    for (k = 4; k < 32; k += 8) {
      dst[k*line] = (short)((g_t32[k][0] * eeo[0] + g_t32[k][1] * eeo[1] + g_t32[k][2] * eeo[2] + g_t32[k][3] * eeo[3] + add) >> shift);
    }
    for (k = 2; k < 32; k += 4) {
      dst[k*line] = (short)((g_t32[k][0] * eo[0] + g_t32[k][1] * eo[1] + g_t32[k][2] * eo[2] + g_t32[k][3] * eo[3] +
        g_t32[k][4] * eo[4] + g_t32[k][5] * eo[5] + g_t32[k][6] * eo[6] + g_t32[k][7] * eo[7] + add) >> shift);
    }
    for (k = 1; k < 32; k += 2) {
      dst[k*line] = (short)((g_t32[k][0] * o[0] + g_t32[k][1] * o[1] + g_t32[k][2] * o[2] + g_t32[k][3] * o[3] +
        g_t32[k][4] * o[4] + g_t32[k][5] * o[5] + g_t32[k][6] * o[6] + g_t32[k][7] * o[7] +
        g_t32[k][8] * o[8] + g_t32[k][9] * o[9] + g_t32[k][10] * o[10] + g_t32[k][11] * o[11] +
        g_t32[k][12] * o[12] + g_t32[k][13] * o[13] + g_t32[k][14] * o[14] + g_t32[k][15] * o[15] + add) >> shift);
    }
    src += 32;
    dst++;
  }
}


static void partial_butterfly_inverse_32_avx2(int16_t *src, int16_t *dst,
  int32_t shift)
{
  int32_t j, k;
  int32_t e[16], o[16];
  int32_t ee[8], eo[8];
  int32_t eee[4], eeo[4];
  int32_t eeee[2], eeeo[2];
  int32_t add = 1 << (shift - 1);
  const int32_t line = 32;

  for (j = 0; j<line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 16; k++) {
      o[k] = g_t32[1][k] * src[line] + g_t32[3][k] * src[3 * line] + g_t32[5][k] * src[5 * line] + g_t32[7][k] * src[7 * line] +
        g_t32[9][k] * src[9 * line] + g_t32[11][k] * src[11 * line] + g_t32[13][k] * src[13 * line] + g_t32[15][k] * src[15 * line] +
        g_t32[17][k] * src[17 * line] + g_t32[19][k] * src[19 * line] + g_t32[21][k] * src[21 * line] + g_t32[23][k] * src[23 * line] +
        g_t32[25][k] * src[25 * line] + g_t32[27][k] * src[27 * line] + g_t32[29][k] * src[29 * line] + g_t32[31][k] * src[31 * line];
    }
    for (k = 0; k < 8; k++) {
      eo[k] = g_t32[2][k] * src[2 * line] + g_t32[6][k] * src[6 * line] + g_t32[10][k] * src[10 * line] + g_t32[14][k] * src[14 * line] +
        g_t32[18][k] * src[18 * line] + g_t32[22][k] * src[22 * line] + g_t32[26][k] * src[26 * line] + g_t32[30][k] * src[30 * line];
    }
    for (k = 0; k < 4; k++) {
      eeo[k] = g_t32[4][k] * src[4 * line] + g_t32[12][k] * src[12 * line] + g_t32[20][k] * src[20 * line] + g_t32[28][k] * src[28 * line];
    }
    eeeo[0] = g_t32[8][0] * src[8 * line] + g_t32[24][0] * src[24 * line];
    eeeo[1] = g_t32[8][1] * src[8 * line] + g_t32[24][1] * src[24 * line];
    eeee[0] = g_t32[0][0] * src[0] + g_t32[16][0] * src[16 * line];
    eeee[1] = g_t32[0][1] * src[0] + g_t32[16][1] * src[16 * line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    eee[0] = eeee[0] + eeeo[0];
    eee[3] = eeee[0] - eeeo[0];
    eee[1] = eeee[1] + eeeo[1];
    eee[2] = eeee[1] - eeeo[1];
    for (k = 0; k < 4; k++) {
      ee[k] = eee[k] + eeo[k];
      ee[k + 4] = eee[3 - k] - eeo[3 - k];
    }
    for (k = 0; k < 8; k++) {
      e[k] = ee[k] + eo[k];
      e[k + 8] = ee[7 - k] - eo[7 - k];
    }
    for (k = 0; k<16; k++) {
      dst[k] = (short)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k + 16] = (short)MAX(-32768, MIN(32767, (e[15 - k] - o[15 - k] + add) >> shift));
    }
    src++;
    dst += 32;
  }
}

#define DCT_NXN_AVX2(n) \
static void dct_ ## n ## x ## n ## _avx2(int8_t bitdepth, int16_t *block, int16_t *coeff) { \
  \
  int16_t tmp[n*n]; \
  int32_t shift_1st = g_convert_to_bit[n] + 1 + (bitdepth - 8); \
  int32_t shift_2nd = g_convert_to_bit[n] + 8; \
  \
  partial_butterfly_ ## n ## _avx2(block, tmp, shift_1st); \
  partial_butterfly_ ## n ## _avx2(tmp, coeff, shift_2nd); \
}

#define IDCT_NXN_AVX2(n) \
static void idct_ ## n ## x ## n ## _avx2(int8_t bitdepth, int16_t *block, int16_t *coeff) { \
\
  int16_t tmp[n*n]; \
  int32_t shift_1st = 7; \
  int32_t shift_2nd = 12 - (bitdepth - 8); \
\
  partial_butterfly_inverse_ ## n ## _avx2(coeff, tmp, shift_1st); \
  partial_butterfly_inverse_ ## n ## _avx2(tmp, block, shift_2nd); \
}

DCT_NXN_AVX2(32);

IDCT_NXN_AVX2(32);

static void matrix_dst_4x4_avx2(int8_t bitdepth, int16_t *src, int16_t *dst)
{
  int32_t shift_1st = g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[4] + 8;
  matrix_transform_2d_4x4_avx2(src, dst, (const int16_t*)dst_4x4, shift_1st, shift_2nd);
}

static void matrix_idst_4x4_avx2(int8_t bitdepth, int16_t *dst, int16_t *src)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  matrix_itransform_2d_4x4_avx2(src, dst, (const int16_t*)dst_4x4, shift_1st, shift_2nd);
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
  matrix_itransform_2d_8x8_avx2(src, dst, (const int16_t*)g_t8, shift_1st, shift_2nd);
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
#endif //COMPILE_INTEL_AVX2

int strategy_register_dct_avx2(void* opaque)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  success &= strategyselector_register(opaque, "fast_forward_dst_4x4", "avx2", 40, &matrix_dst_4x4_avx2);

  success &= strategyselector_register(opaque, "dct_4x4", "avx2", 40, &matrix_dct_4x4_avx2);
  success &= strategyselector_register(opaque, "dct_8x8", "avx2", 40, &matrix_dct_8x8_avx2);
  success &= strategyselector_register(opaque, "dct_16x16", "avx2", 40, &matrix_dct_16x16_avx2);
  success &= strategyselector_register(opaque, "dct_32x32", "avx2", 40, &dct_32x32_avx2);

  success &= strategyselector_register(opaque, "fast_inverse_dst_4x4", "avx2", 40, &matrix_idst_4x4_avx2);

  success &= strategyselector_register(opaque, "idct_4x4", "avx2", 40, &matrix_idct_4x4_avx2);
  success &= strategyselector_register(opaque, "idct_8x8", "avx2", 40, &matrix_idct_8x8_avx2);
  success &= strategyselector_register(opaque, "idct_16x16", "avx2", 40, &matrix_idct_16x16_avx2);
  success &= strategyselector_register(opaque, "idct_32x32", "avx2", 40, &idct_32x32_avx2);
#endif //COMPILE_INTEL_AVX2  
  return success;
}
