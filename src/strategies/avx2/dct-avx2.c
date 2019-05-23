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

#include "strategies/avx2/dct-avx2.h"

#if COMPILE_INTEL_AVX2
#include <immintrin.h>

#include "strategyselector.h"
#include "tables.h"

extern const int16_t kvz_g_dst_4[4][4];
extern const int16_t kvz_g_dct_4[4][4];
extern const int16_t kvz_g_dct_8[8][8];
extern const int16_t kvz_g_dct_16[16][16];
extern const int16_t kvz_g_dct_32[32][32];

extern const int16_t kvz_g_dst_4_t[4][4];
extern const int16_t kvz_g_dct_4_t[4][4];
extern const int16_t kvz_g_dct_8_t[8][8];
extern const int16_t kvz_g_dct_16_t[16][16];
extern const int16_t kvz_g_dct_32_t[32][32];

/*
* \file
* \brief AVX2 transformations.
*/

// 4x4 matrix multiplication with value clipping.
// Parameters: Two 4x4 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_4x4_avx2(const int16_t *left, const int16_t *right, int16_t *dst, int32_t shift)
{
  __m256i b[2], a, result, even[2], odd[2];

  const int32_t add = 1 << (shift - 1);

  a = _mm256_loadu_si256((__m256i*) left);
  b[0] = _mm256_loadu_si256((__m256i*) right);

  // Interleave values in both 128-bit lanes
  b[0] = _mm256_unpacklo_epi16(b[0], _mm256_srli_si256(b[0], 8));
  b[1] = _mm256_permute2x128_si256(b[0], b[0], 1 + 16);
  b[0] = _mm256_permute2x128_si256(b[0], b[0], 0);

  // Fill both 128-lanes with the first pair of 16-bit factors in the lane.
  even[0] = _mm256_shuffle_epi32(a, 0);
  odd[0] = _mm256_shuffle_epi32(a, 1 + 4 + 16 + 64);

  // Multiply packed elements and sum pairs. Input 16-bit output 32-bit.
  even[0] = _mm256_madd_epi16(even[0], b[0]);
  odd[0] = _mm256_madd_epi16(odd[0], b[1]);

  // Add the halves of the dot product and
  // round.
  result = _mm256_add_epi32(even[0], odd[0]);
  result = _mm256_add_epi32(result, _mm256_set1_epi32(add));
  result = _mm256_srai_epi32(result, shift);

  //Repeat for the remaining parts
  even[1] = _mm256_shuffle_epi32(a, 2 + 8 + 32 + 128);
  odd[1] = _mm256_shuffle_epi32(a, 3 + 12 + 48 + 192);

  even[1] = _mm256_madd_epi16(even[1], b[0]);
  odd[1] = _mm256_madd_epi16(odd[1], b[1]);

  odd[1] = _mm256_add_epi32(even[1], odd[1]);
  odd[1] = _mm256_add_epi32(odd[1], _mm256_set1_epi32(add));
  odd[1] = _mm256_srai_epi32(odd[1], shift);

  // Truncate to 16-bit values
  result = _mm256_packs_epi32(result, odd[1]);

  _mm256_storeu_si256((__m256i*)dst, result);
}

// 8x8 matrix multiplication with value clipping.
// Parameters: Two 8x8 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
//
static void mul_clip_matrix_8x8_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i b[2], accu[8], even[2], odd[2];

  const int32_t add = 1 << (shift - 1);

  b[0] = _mm256_loadu_si256((__m256i*) right);

  b[1] = _mm256_unpackhi_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
  b[0] = _mm256_unpacklo_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
  b[0] = _mm256_inserti128_si256(b[0], _mm256_castsi256_si128(b[1]), 1);

  for (i = 0; i < 8; i += 2) {

    even[0] = _mm256_set1_epi32(((int32_t*)left)[4 * i]);
    even[0] = _mm256_madd_epi16(even[0], b[0]);
    accu[i] = even[0];

    odd[0] = _mm256_set1_epi32(((int32_t*)left)[4 * (i + 1)]);
    odd[0] = _mm256_madd_epi16(odd[0], b[0]);
    accu[i + 1] = odd[0];
  }

  for (j = 1; j < 4; ++j) {

    b[0] = _mm256_loadu_si256((__m256i*)right + j);

    b[1] = _mm256_unpackhi_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
    b[0] = _mm256_unpacklo_epi16(b[0], _mm256_castsi128_si256(_mm256_extracti128_si256(b[0], 1)));
    b[0] = _mm256_inserti128_si256(b[0], _mm256_castsi256_si128(b[1]), 1);

    for (i = 0; i < 8; i += 2) {

      even[0] = _mm256_set1_epi32(((int32_t*)left)[4 * i + j]);
      even[0] = _mm256_madd_epi16(even[0], b[0]);
      accu[i] = _mm256_add_epi32(accu[i], even[0]);

      odd[0] = _mm256_set1_epi32(((int32_t*)left)[4 * (i + 1) + j]);
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

static INLINE __m256i swap_lanes(__m256i v)
{
  return _mm256_permute4x64_epi64(v, _MM_SHUFFLE(1, 0, 3, 2));
}

static INLINE __m256i truncate(__m256i v, __m256i debias, int32_t shift)
{
  __m256i truncable = _mm256_add_epi32 (v,         debias);
  return              _mm256_srai_epi32(truncable, shift);
}

static void matrix_dct_8x8_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[8] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[8] + 8;

  const int32_t add1    = 1 << (shift_1st - 1);
  const __m256i debias1 = _mm256_set1_epi32(add1);

  const int32_t add2    = 1 << (shift_2nd - 1);
  const __m256i debias2 = _mm256_set1_epi32(add2);

  const __m256i *dct    = (__m256i *)&(kvz_g_dct_8[0][0]);

  // Keep upper row intact and swap neighboring 16-bit words in lower row
  const __m256i shuf_lorow_mask =
      _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                       8,  9,  10, 11, 12, 13, 14, 15,
                       18, 19, 16, 17, 22, 23, 20, 21,
                       26, 27, 24, 25, 30, 31, 28, 29);

  __m256i tmpres[4];

  // Dual Rows, because two 8x16b words fit in one YMM
  __m256i i_dr_0      = _mm256_loadu_si256((__m256i *)input + 0);
  __m256i i_dr_1      = _mm256_loadu_si256((__m256i *)input + 1);
  __m256i i_dr_2      = _mm256_loadu_si256((__m256i *)input + 2);
  __m256i i_dr_3      = _mm256_loadu_si256((__m256i *)input + 3);

  __m256i i_dr_0_swp  = swap_lanes(i_dr_0);
  __m256i i_dr_1_swp  = swap_lanes(i_dr_1);
  __m256i i_dr_2_swp  = swap_lanes(i_dr_2);
  __m256i i_dr_3_swp  = swap_lanes(i_dr_3);

  /*
   * Multiply input by the tranpose of DCT matrix into tmpres, and DCT matrix
   * by tmpres - this is then our output matrix
   *
   * It's easier to implement an AVX2 matrix multiplication if you can multiply
   * the left term with the transpose of the right term. Here things are stored
   * row-wise, not column-wise, so we can effectively read DCT_T column-wise
   * into YMM registers by reading DCT row-wise. Also because of this, the
   * first multiplication is hacked to produce the transpose of the result
   * instead, since it will be used in similar fashion as the right operand
   * in the second multiplication.
   */
  for (int dry = 0; dry < 4; dry++) {

    // Read columns of DCT matrix's transpose by reading rows of DCT matrix
    __m256i d_dr        = _mm256_loadu_si256(dct + dry);

    __m256i prod0       = _mm256_madd_epi16(d_dr,     i_dr_0);
    __m256i prod0_swp   = _mm256_madd_epi16(d_dr,     i_dr_0_swp);
    __m256i prod1       = _mm256_madd_epi16(d_dr,     i_dr_1);
    __m256i prod1_swp   = _mm256_madd_epi16(d_dr,     i_dr_1_swp);
    __m256i prod2       = _mm256_madd_epi16(d_dr,     i_dr_2);
    __m256i prod2_swp   = _mm256_madd_epi16(d_dr,     i_dr_2_swp);
    __m256i prod3       = _mm256_madd_epi16(d_dr,     i_dr_3);
    __m256i prod3_swp   = _mm256_madd_epi16(d_dr,     i_dr_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate(hsum2c_0, debias1, shift_1st);
    __m256i hsum2c_1_tr = truncate(hsum2c_1, debias1, shift_1st);

    __m256i tmp_dr      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    tmpres[dry]         = _mm256_shuffle_epi8(tmp_dr, shuf_lorow_mask);
  }

  __m256i t_dr_0      = tmpres[0];
  __m256i t_dr_1      = tmpres[1];
  __m256i t_dr_2      = tmpres[2];
  __m256i t_dr_3      = tmpres[3];

  __m256i t_dr_0_swp  = swap_lanes(t_dr_0);
  __m256i t_dr_1_swp  = swap_lanes(t_dr_1);
  __m256i t_dr_2_swp  = swap_lanes(t_dr_2);
  __m256i t_dr_3_swp  = swap_lanes(t_dr_3);

  for (int dry = 0; dry < 4; dry++) {
    __m256i d_dr        = _mm256_loadu_si256(dct + dry);

    __m256i prod0       = _mm256_madd_epi16(d_dr,     t_dr_0);
    __m256i prod0_swp   = _mm256_madd_epi16(d_dr,     t_dr_0_swp);
    __m256i prod1       = _mm256_madd_epi16(d_dr,     t_dr_1);
    __m256i prod1_swp   = _mm256_madd_epi16(d_dr,     t_dr_1_swp);
    __m256i prod2       = _mm256_madd_epi16(d_dr,     t_dr_2);
    __m256i prod2_swp   = _mm256_madd_epi16(d_dr,     t_dr_2_swp);
    __m256i prod3       = _mm256_madd_epi16(d_dr,     t_dr_3);
    __m256i prod3_swp   = _mm256_madd_epi16(d_dr,     t_dr_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate(hsum2c_0, debias2, shift_2nd);
    __m256i hsum2c_1_tr = truncate(hsum2c_1, debias2, shift_2nd);

    __m256i tmp_dr      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    __m256i final_dr    = _mm256_shuffle_epi8(tmp_dr, shuf_lorow_mask);

    _mm256_storeu_si256(((__m256i *)output) + dry, final_dr);
  }
}

// 16x16 matrix multiplication with value clipping.
// Parameters: Two 16x16 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_16x16_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i row[4], accu[16][2], even, odd;

  const int32_t stride = 8;

  const int32_t add = 1 << (shift - 1);

  row[0] = _mm256_loadu_si256((__m256i*) right);
  row[1] = _mm256_loadu_si256((__m256i*) right + 1);
  row[2] = _mm256_unpacklo_epi16(row[0], row[1]);
  row[3] = _mm256_unpackhi_epi16(row[0], row[1]);
  row[0] = _mm256_permute2x128_si256(row[2], row[3], 0 + 32);
  row[1] = _mm256_permute2x128_si256(row[2], row[3], 1 + 48);

  for (i = 0; i < 16; i += 2) {

    even = _mm256_set1_epi32(((int32_t*)left)[stride * i]);
    accu[i][0] = _mm256_madd_epi16(even, row[0]);
    accu[i][1] = _mm256_madd_epi16(even, row[1]);

    odd = _mm256_set1_epi32(((int32_t*)left)[stride * (i + 1)]);
    accu[i + 1][0] = _mm256_madd_epi16(odd, row[0]);
    accu[i + 1][1] = _mm256_madd_epi16(odd, row[1]);
  }

  for (j = 2; j < 16; j += 2) {

    row[0] = _mm256_loadu_si256((__m256i*)right + j);
    row[1] = _mm256_loadu_si256((__m256i*)right + j + 1);
    row[2] = _mm256_unpacklo_epi16(row[0], row[1]);
    row[3] = _mm256_unpackhi_epi16(row[0], row[1]);
    row[0] = _mm256_permute2x128_si256(row[2], row[3], 0 + 32);
    row[1] = _mm256_permute2x128_si256(row[2], row[3], 1 + 48);

    for (i = 0; i < 16; i += 2) {

      even = _mm256_set1_epi32(((int32_t*)left)[stride * i + j / 2]);
      accu[i][0] = _mm256_add_epi32(accu[i][0], _mm256_madd_epi16(even, row[0]));
      accu[i][1] = _mm256_add_epi32(accu[i][1], _mm256_madd_epi16(even, row[1]));

      odd = _mm256_set1_epi32(((int32_t*)left)[stride * (i + 1) + j / 2]);
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

// 32x32 matrix multiplication with value clipping.
// Parameters: Two 32x32 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_32x32_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  int i, j;
  __m256i row[4], tmp[2], accu[32][4], even, odd;

  const int32_t stride = 16;

  const int32_t add = 1 << (shift - 1);

  row[0] = _mm256_loadu_si256((__m256i*) right);
  row[1] = _mm256_loadu_si256((__m256i*) right + 2);
  tmp[0] = _mm256_unpacklo_epi16(row[0], row[1]);
  tmp[1] = _mm256_unpackhi_epi16(row[0], row[1]);
  row[0] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
  row[1] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

  row[2] = _mm256_loadu_si256((__m256i*) right + 1);
  row[3] = _mm256_loadu_si256((__m256i*) right + 3);
  tmp[0] = _mm256_unpacklo_epi16(row[2], row[3]);
  tmp[1] = _mm256_unpackhi_epi16(row[2], row[3]);
  row[2] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
  row[3] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

  for (i = 0; i < 32; i += 2) {

    even = _mm256_set1_epi32(((int32_t*)left)[stride * i]);
    accu[i][0] = _mm256_madd_epi16(even, row[0]);
    accu[i][1] = _mm256_madd_epi16(even, row[1]);
    accu[i][2] = _mm256_madd_epi16(even, row[2]);
    accu[i][3] = _mm256_madd_epi16(even, row[3]);

    odd = _mm256_set1_epi32(((int32_t*)left)[stride * (i + 1)]);
    accu[i + 1][0] = _mm256_madd_epi16(odd, row[0]);
    accu[i + 1][1] = _mm256_madd_epi16(odd, row[1]);
    accu[i + 1][2] = _mm256_madd_epi16(odd, row[2]);
    accu[i + 1][3] = _mm256_madd_epi16(odd, row[3]);
  }

  for (j = 4; j < 64; j += 4) {

    row[0] = _mm256_loadu_si256((__m256i*)right + j);
    row[1] = _mm256_loadu_si256((__m256i*)right + j + 2);
    tmp[0] = _mm256_unpacklo_epi16(row[0], row[1]);
    tmp[1] = _mm256_unpackhi_epi16(row[0], row[1]);
    row[0] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
    row[1] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

    row[2] = _mm256_loadu_si256((__m256i*) right + j + 1);
    row[3] = _mm256_loadu_si256((__m256i*) right + j + 3);
    tmp[0] = _mm256_unpacklo_epi16(row[2], row[3]);
    tmp[1] = _mm256_unpackhi_epi16(row[2], row[3]);
    row[2] = _mm256_permute2x128_si256(tmp[0], tmp[1], 0 + 32);
    row[3] = _mm256_permute2x128_si256(tmp[0], tmp[1], 1 + 48);

    for (i = 0; i < 32; i += 2) {

      even = _mm256_set1_epi32(((int32_t*)left)[stride * i + j / 4]);
      accu[i][0] = _mm256_add_epi32(accu[i][0], _mm256_madd_epi16(even, row[0]));
      accu[i][1] = _mm256_add_epi32(accu[i][1], _mm256_madd_epi16(even, row[1]));
      accu[i][2] = _mm256_add_epi32(accu[i][2], _mm256_madd_epi16(even, row[2]));
      accu[i][3] = _mm256_add_epi32(accu[i][3], _mm256_madd_epi16(even, row[3]));

      odd = _mm256_set1_epi32(((int32_t*)left)[stride * (i + 1) + j / 4]);
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

// Macro that generates 2D transform functions with clipping values.
// Sets correct shift values and matrices according to transform type and
// block size. Performs matrix multiplication horizontally and vertically.
#define TRANSFORM(type, n) static void matrix_ ## type ## _ ## n ## x ## n ## _avx2(int8_t bitdepth, const int16_t *input, int16_t *output)\
{\
  int32_t shift_1st = kvz_g_convert_to_bit[n] + 1 + (bitdepth - 8); \
  int32_t shift_2nd = kvz_g_convert_to_bit[n] + 8; \
  int16_t tmp[n * n];\
  const int16_t *tdct = &kvz_g_ ## type ## _ ## n ## _t[0][0];\
  const int16_t *dct = &kvz_g_ ## type ## _ ## n [0][0];\
\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(input, tdct, tmp, shift_1st);\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(dct, tmp, output, shift_2nd);\
}\

// Macro that generates 2D inverse transform functions with clipping values.
// Sets correct shift values and matrices according to transform type and
// block size. Performs matrix multiplication horizontally and vertically.
#define ITRANSFORM(type, n) \
static void matrix_i ## type ## _## n ## x ## n ## _avx2(int8_t bitdepth, const int16_t *input, int16_t *output)\
{\
  int32_t shift_1st = 7; \
  int32_t shift_2nd = 12 - (bitdepth - 8); \
  int16_t tmp[n * n];\
  const int16_t *tdct = &kvz_g_ ## type ## _ ## n ## _t[0][0];\
  const int16_t *dct = &kvz_g_ ## type ## _ ## n [0][0];\
\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tdct, input, tmp, shift_1st);\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tmp, dct, output, shift_2nd);\
}\

// Generate all the transform functions
TRANSFORM(dst, 4);
TRANSFORM(dct, 4);

// Ha, we've got a tailored implementation for this
// TRANSFORM(dct, 8);

TRANSFORM(dct, 16);
TRANSFORM(dct, 32);

ITRANSFORM(dst, 4);
ITRANSFORM(dct, 4);
ITRANSFORM(dct, 8);
ITRANSFORM(dct, 16);
ITRANSFORM(dct, 32);

#endif //COMPILE_INTEL_AVX2

int kvz_strategy_register_dct_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "fast_forward_dst_4x4", "avx2", 40, &matrix_dst_4x4_avx2);

    success &= kvz_strategyselector_register(opaque, "dct_4x4", "avx2", 40, &matrix_dct_4x4_avx2);
    success &= kvz_strategyselector_register(opaque, "dct_8x8", "avx2", 40, &matrix_dct_8x8_avx2);
    success &= kvz_strategyselector_register(opaque, "dct_16x16", "avx2", 40, &matrix_dct_16x16_avx2);
    success &= kvz_strategyselector_register(opaque, "dct_32x32", "avx2", 40, &matrix_dct_32x32_avx2);

    success &= kvz_strategyselector_register(opaque, "fast_inverse_dst_4x4", "avx2", 40, &matrix_idst_4x4_avx2);

    success &= kvz_strategyselector_register(opaque, "idct_4x4", "avx2", 40, &matrix_idct_4x4_avx2);
    success &= kvz_strategyselector_register(opaque, "idct_8x8", "avx2", 40, &matrix_idct_8x8_avx2);
    success &= kvz_strategyselector_register(opaque, "idct_16x16", "avx2", 40, &matrix_idct_16x16_avx2);
    success &= kvz_strategyselector_register(opaque, "idct_32x32", "avx2", 40, &matrix_idct_32x32_avx2);
  }
#endif //COMPILE_INTEL_AVX2  
  return success;
}
