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

static INLINE __m256i swap_lanes(__m256i v)
{
  return _mm256_permute4x64_epi64(v, _MM_SHUFFLE(1, 0, 3, 2));
}

static INLINE __m256i truncate(__m256i v, __m256i debias, int32_t shift)
{
  __m256i truncable = _mm256_add_epi32 (v,         debias);
  return              _mm256_srai_epi32(truncable, shift);
}

// 4x4 matrix multiplication with value clipping.
// Parameters: Two 4x4 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static __m256i mul_clip_matrix_4x4_avx2(const __m256i left, const __m256i right, int shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i right_los = _mm256_permute4x64_epi64(right, _MM_SHUFFLE(2, 0, 2, 0));
  __m256i right_his = _mm256_permute4x64_epi64(right, _MM_SHUFFLE(3, 1, 3, 1));

  __m256i right_cols_up = _mm256_unpacklo_epi16(right_los, right_his);
  __m256i right_cols_dn = _mm256_unpackhi_epi16(right_los, right_his);

  __m256i left_slice1 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i left_slice2 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(1, 1, 1, 1));
  __m256i left_slice3 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(2, 2, 2, 2));
  __m256i left_slice4 = _mm256_shuffle_epi32(left, _MM_SHUFFLE(3, 3, 3, 3));

  __m256i prod1 = _mm256_madd_epi16(left_slice1, right_cols_up);
  __m256i prod2 = _mm256_madd_epi16(left_slice2, right_cols_dn);
  __m256i prod3 = _mm256_madd_epi16(left_slice3, right_cols_up);
  __m256i prod4 = _mm256_madd_epi16(left_slice4, right_cols_dn);

  __m256i rows_up = _mm256_add_epi32(prod1, prod2);
  __m256i rows_dn = _mm256_add_epi32(prod3, prod4);

  __m256i rows_up_tr = truncate(rows_up, debias, shift);
  __m256i rows_dn_tr = truncate(rows_dn, debias, shift);

  __m256i result = _mm256_packs_epi32(rows_up_tr, rows_dn_tr);
  return result;
}

static void matrix_dst_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[4] + 8;
  const int16_t *tdst = &kvz_g_dst_4_t[0][0];
  const int16_t *dst  = &kvz_g_dst_4  [0][0];

  __m256i tdst_v = _mm256_load_si256((const __m256i *) tdst);
  __m256i  dst_v = _mm256_load_si256((const __m256i *)  dst);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(in_v,  tdst_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(dst_v, tmp,    shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_idst_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);

  const int16_t *tdst = &kvz_g_dst_4_t[0][0];
  const int16_t *dst  = &kvz_g_dst_4  [0][0];

  __m256i tdst_v = _mm256_load_si256((const __m256i *)tdst);
  __m256i  dst_v = _mm256_load_si256((const __m256i *) dst);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(tdst_v, in_v,  shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(tmp,    dst_v, shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_dct_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[4] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[4] + 8;
  const int16_t *tdct = &kvz_g_dct_4_t[0][0];
  const int16_t *dct  = &kvz_g_dct_4  [0][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i *) tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i *)  dct);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(in_v,  tdct_v, shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(dct_v, tmp,    shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void matrix_idct_4x4_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);

  const int16_t *tdct = &kvz_g_dct_4_t[0][0];
  const int16_t *dct  = &kvz_g_dct_4  [0][0];

  __m256i tdct_v = _mm256_load_si256((const __m256i *)tdct);
  __m256i  dct_v = _mm256_load_si256((const __m256i *) dct);
  __m256i   in_v = _mm256_load_si256((const __m256i *)input);

  __m256i tmp    = mul_clip_matrix_4x4_avx2(tdct_v, in_v,  shift_1st);
  __m256i result = mul_clip_matrix_4x4_avx2(tmp,    dct_v, shift_2nd);

  _mm256_store_si256((__m256i *)output, result);
}

static void mul_clip_matrix_8x8_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  const __m256i transp_mask = _mm256_broadcastsi128_si256(_mm_setr_epi8(0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15));

  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i left_dr[4] = {
    _mm256_load_si256((const __m256i *)left + 0),
    _mm256_load_si256((const __m256i *)left + 1),
    _mm256_load_si256((const __m256i *)left + 2),
    _mm256_load_si256((const __m256i *)left + 3),
  };
  __m256i right_dr[4] = {
    _mm256_load_si256((const __m256i *)right + 0),
    _mm256_load_si256((const __m256i *)right + 1),
    _mm256_load_si256((const __m256i *)right + 2),
    _mm256_load_si256((const __m256i *)right + 3),
  };

  __m256i rdrs_rearr[8];

  // Rearrange right matrix
  for (int32_t dry = 0; dry < 4; dry++) {
    __m256i rdr = right_dr[dry];
    __m256i rdr_los = _mm256_permute4x64_epi64(rdr, _MM_SHUFFLE(2, 0, 2, 0));
    __m256i rdr_his = _mm256_permute4x64_epi64(rdr, _MM_SHUFFLE(3, 1, 3, 1));

    __m256i rdr_lo_rearr = _mm256_shuffle_epi8(rdr_los, transp_mask);
    __m256i rdr_hi_rearr = _mm256_shuffle_epi8(rdr_his, transp_mask);

    rdrs_rearr[dry * 2 + 0] = rdr_lo_rearr;
    rdrs_rearr[dry * 2 + 1] = rdr_hi_rearr;
  }

  // Double-Row Y for destination matrix
  for (int32_t dry = 0; dry < 4; dry++) {
    __m256i ldr = left_dr[dry];

    __m256i ldr_slice12 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(0, 0, 0, 0));
    __m256i ldr_slice34 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(1, 1, 1, 1));
    __m256i ldr_slice56 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(2, 2, 2, 2));
    __m256i ldr_slice78 = _mm256_shuffle_epi32(ldr, _MM_SHUFFLE(3, 3, 3, 3));

    __m256i prod1 = _mm256_madd_epi16(ldr_slice12, rdrs_rearr[0]);
    __m256i prod2 = _mm256_madd_epi16(ldr_slice12, rdrs_rearr[1]);
    __m256i prod3 = _mm256_madd_epi16(ldr_slice34, rdrs_rearr[2]);
    __m256i prod4 = _mm256_madd_epi16(ldr_slice34, rdrs_rearr[3]);
    __m256i prod5 = _mm256_madd_epi16(ldr_slice56, rdrs_rearr[4]);
    __m256i prod6 = _mm256_madd_epi16(ldr_slice56, rdrs_rearr[5]);
    __m256i prod7 = _mm256_madd_epi16(ldr_slice78, rdrs_rearr[6]);
    __m256i prod8 = _mm256_madd_epi16(ldr_slice78, rdrs_rearr[7]);

    __m256i lo_1 = _mm256_add_epi32(prod1, prod3);
    __m256i hi_1 = _mm256_add_epi32(prod2, prod4);
    __m256i lo_2 = _mm256_add_epi32(prod5, prod7);
    __m256i hi_2 = _mm256_add_epi32(prod6, prod8);

    __m256i lo   = _mm256_add_epi32(lo_1,  lo_2);
    __m256i hi   = _mm256_add_epi32(hi_1,  hi_2);

    __m256i lo_tr = truncate(lo, debias, shift);
    __m256i hi_tr = truncate(hi, debias, shift);

    __m256i final_dr = _mm256_packs_epi32(lo_tr, hi_tr);

    _mm256_store_si256((__m256i *)dst + dry, final_dr);
  }
}

// Multiplies A by B_T's transpose and stores result's transpose in output,
// which should be an array of 4 __m256i's
static void matmul_8x8_a_bt_t(const int16_t *a, const int16_t *b_t,
    __m256i *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  // Keep upper row intact and swap neighboring 16-bit words in lower row
  const __m256i shuf_lorow_mask =
      _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                       8,  9,  10, 11, 12, 13, 14, 15,
                       18, 19, 16, 17, 22, 23, 20, 21,
                       26, 27, 24, 25, 30, 31, 28, 29);

  const __m256i *b_t_256 = (const __m256i *)b_t;

  // Dual Rows, because two 8x16b words fit in one YMM
  __m256i a_dr_0      = _mm256_load_si256((__m256i *)a + 0);
  __m256i a_dr_1      = _mm256_load_si256((__m256i *)a + 1);
  __m256i a_dr_2      = _mm256_load_si256((__m256i *)a + 2);
  __m256i a_dr_3      = _mm256_load_si256((__m256i *)a + 3);

  __m256i a_dr_0_swp  = swap_lanes(a_dr_0);
  __m256i a_dr_1_swp  = swap_lanes(a_dr_1);
  __m256i a_dr_2_swp  = swap_lanes(a_dr_2);
  __m256i a_dr_3_swp  = swap_lanes(a_dr_3);

  for (int dry = 0; dry < 4; dry++) {

    // Read dual columns of B matrix by reading rows of its transpose
    __m256i b_dc        = _mm256_load_si256(b_t_256 + dry);

    __m256i prod0       = _mm256_madd_epi16(b_dc,     a_dr_0);
    __m256i prod0_swp   = _mm256_madd_epi16(b_dc,     a_dr_0_swp);
    __m256i prod1       = _mm256_madd_epi16(b_dc,     a_dr_1);
    __m256i prod1_swp   = _mm256_madd_epi16(b_dc,     a_dr_1_swp);
    __m256i prod2       = _mm256_madd_epi16(b_dc,     a_dr_2);
    __m256i prod2_swp   = _mm256_madd_epi16(b_dc,     a_dr_2_swp);
    __m256i prod3       = _mm256_madd_epi16(b_dc,     a_dr_3);
    __m256i prod3_swp   = _mm256_madd_epi16(b_dc,     a_dr_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate(hsum2c_0, debias, shift);
    __m256i hsum2c_1_tr = truncate(hsum2c_1, debias, shift);

    __m256i tmp_dc      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    output[dry]         = _mm256_shuffle_epi8(tmp_dc, shuf_lorow_mask);
  }
}

// Multiplies A by B_T's transpose and stores result in output
// which should be an array of 4 __m256i's
static void matmul_8x8_a_bt(const int16_t *a, const __m256i *b_t,
    int16_t *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const __m256i shuf_lorow_mask =
      _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                       8,  9,  10, 11, 12, 13, 14, 15,
                       18, 19, 16, 17, 22, 23, 20, 21,
                       26, 27, 24, 25, 30, 31, 28, 29);

  const __m256i *a_256 = (const __m256i *)a;

  __m256i b_dc_0      = b_t[0];
  __m256i b_dc_1      = b_t[1];
  __m256i b_dc_2      = b_t[2];
  __m256i b_dc_3      = b_t[3];

  __m256i b_dc_0_swp  = swap_lanes(b_dc_0);
  __m256i b_dc_1_swp  = swap_lanes(b_dc_1);
  __m256i b_dc_2_swp  = swap_lanes(b_dc_2);
  __m256i b_dc_3_swp  = swap_lanes(b_dc_3);

  for (int dry = 0; dry < 4; dry++) {
    __m256i a_dr        = _mm256_load_si256(a_256 + dry);

    __m256i prod0       = _mm256_madd_epi16(a_dr,     b_dc_0);
    __m256i prod0_swp   = _mm256_madd_epi16(a_dr,     b_dc_0_swp);
    __m256i prod1       = _mm256_madd_epi16(a_dr,     b_dc_1);
    __m256i prod1_swp   = _mm256_madd_epi16(a_dr,     b_dc_1_swp);
    __m256i prod2       = _mm256_madd_epi16(a_dr,     b_dc_2);
    __m256i prod2_swp   = _mm256_madd_epi16(a_dr,     b_dc_2_swp);
    __m256i prod3       = _mm256_madd_epi16(a_dr,     b_dc_3);
    __m256i prod3_swp   = _mm256_madd_epi16(a_dr,     b_dc_3_swp);

    __m256i hsum0       = _mm256_hadd_epi32(prod0,    prod0_swp);
    __m256i hsum1       = _mm256_hadd_epi32(prod1,    prod1_swp);
    __m256i hsum2       = _mm256_hadd_epi32(prod2,    prod2_swp);
    __m256i hsum3       = _mm256_hadd_epi32(prod3,    prod3_swp);

    __m256i hsum2c_0    = _mm256_hadd_epi32(hsum0,    hsum1);
    __m256i hsum2c_1    = _mm256_hadd_epi32(hsum2,    hsum3);

    __m256i hsum2c_0_tr = truncate(hsum2c_0, debias, shift);
    __m256i hsum2c_1_tr = truncate(hsum2c_1, debias, shift);

    __m256i tmp_dr      = _mm256_packs_epi32(hsum2c_0_tr, hsum2c_1_tr);

    __m256i final_dr    = _mm256_shuffle_epi8(tmp_dr, shuf_lorow_mask);

    _mm256_store_si256((__m256i *)output + dry, final_dr);
  }
}

static void matrix_dct_8x8_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[8] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[8] + 8;

  const int16_t *dct  = &kvz_g_dct_8[0][0];

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

  __m256i tmpres[4];

  matmul_8x8_a_bt_t(input,  dct, tmpres, shift_1st);
  matmul_8x8_a_bt  (dct, tmpres, output, shift_2nd);
}

static void matrix_idct_8x8_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[8 * 8];

  const int16_t *tdct = &kvz_g_dct_8_t[0][0];
  const int16_t *dct  = &kvz_g_dct_8  [0][0];

  mul_clip_matrix_8x8_avx2(tdct, input, tmp,    shift_1st);
  mul_clip_matrix_8x8_avx2(tmp,  dct,   output, shift_2nd);

  /*
   * Because:
   * out = tdct * input * dct = tdct * (input * dct) = tdct * (input * transpose(tdct))
   * This could almost be done this way:
   *
   * matmul_8x8_a_bt_t(input, tdct, debias1, shift_1st, tmp);
   * matmul_8x8_a_bt  (tdct,  tmp,  debias2, shift_2nd, output);
   *
   * But not really, since it will fall victim to some very occasional
   * rounding errors. Sadly.
   */
}

static void matmul_16x16_a_bt_t(const int16_t *a, const int16_t *b_t, __m256i *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  for (int32_t x = 0; x < 16; x++) {
    __m256i bt_c = _mm256_load_si256((const __m256i *)b_t + x);

    __m256i results_32[2];

    // First Row Offset
    for (int32_t fro = 0; fro < 2; fro++) {
      // Read first rows 0, 1, 2, 3, 8, 9, 10, 11, and then next 4
      __m256i a_r0  = _mm256_load_si256((const __m256i *)a + fro * 4 + 0);
      __m256i a_r1  = _mm256_load_si256((const __m256i *)a + fro * 4 + 1);
      __m256i a_r2  = _mm256_load_si256((const __m256i *)a + fro * 4 + 2);
      __m256i a_r3  = _mm256_load_si256((const __m256i *)a + fro * 4 + 3);
      __m256i a_r8  = _mm256_load_si256((const __m256i *)a + fro * 4 + 8);
      __m256i a_r9  = _mm256_load_si256((const __m256i *)a + fro * 4 + 9);
      __m256i a_r10 = _mm256_load_si256((const __m256i *)a + fro * 4 + 10);
      __m256i a_r11 = _mm256_load_si256((const __m256i *)a + fro * 4 + 11);

      __m256i p0  = _mm256_madd_epi16(bt_c, a_r0);
      __m256i p1  = _mm256_madd_epi16(bt_c, a_r1);
      __m256i p2  = _mm256_madd_epi16(bt_c, a_r2);
      __m256i p3  = _mm256_madd_epi16(bt_c, a_r3);
      __m256i p8  = _mm256_madd_epi16(bt_c, a_r8);
      __m256i p9  = _mm256_madd_epi16(bt_c, a_r9);
      __m256i p10 = _mm256_madd_epi16(bt_c, a_r10);
      __m256i p11 = _mm256_madd_epi16(bt_c, a_r11);

      // Combine low lanes from P0 and P8, high lanes from them, and the same
      // with P1:P9 and so on
      __m256i p0l = _mm256_permute2x128_si256(p0, p8,  0x20);
      __m256i p0h = _mm256_permute2x128_si256(p0, p8,  0x31);
      __m256i p1l = _mm256_permute2x128_si256(p1, p9,  0x20);
      __m256i p1h = _mm256_permute2x128_si256(p1, p9,  0x31);
      __m256i p2l = _mm256_permute2x128_si256(p2, p10, 0x20);
      __m256i p2h = _mm256_permute2x128_si256(p2, p10, 0x31);
      __m256i p3l = _mm256_permute2x128_si256(p3, p11, 0x20);
      __m256i p3h = _mm256_permute2x128_si256(p3, p11, 0x31);

      __m256i s0  = _mm256_add_epi32(p0l, p0h);
      __m256i s1  = _mm256_add_epi32(p1l, p1h);
      __m256i s2  = _mm256_add_epi32(p2l, p2h);
      __m256i s3  = _mm256_add_epi32(p3l, p3h);

      __m256i s4  = _mm256_unpacklo_epi64(s0, s1);
      __m256i s5  = _mm256_unpackhi_epi64(s0, s1);
      __m256i s6  = _mm256_unpacklo_epi64(s2, s3);
      __m256i s7  = _mm256_unpackhi_epi64(s2, s3);

      __m256i s8  = _mm256_add_epi32(s4, s5);
      __m256i s9  = _mm256_add_epi32(s6, s7);

      __m256i res = _mm256_hadd_epi32(s8, s9);
      results_32[fro] = truncate(res, debias, shift);
    }
    __m256i final_col = _mm256_packs_epi32(results_32[0], results_32[1]);
    output[x] = final_col;
  }
}

static void matmul_16x16_a_bt(const int16_t *a, const __m256i *b_t, int16_t *output, const int8_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  for (int32_t y = 0; y < 16; y++) {
    __m256i a_r = _mm256_load_si256((const __m256i *)a + y);
    __m256i results_32[2];

    for (int32_t fco = 0; fco < 2; fco++) {
      // Read first cols 0, 1, 2, 3, 8, 9, 10, 11, and then next 4
      __m256i bt_c0  = b_t[fco * 4 + 0];
      __m256i bt_c1  = b_t[fco * 4 + 1];
      __m256i bt_c2  = b_t[fco * 4 + 2];
      __m256i bt_c3  = b_t[fco * 4 + 3];
      __m256i bt_c8  = b_t[fco * 4 + 8];
      __m256i bt_c9  = b_t[fco * 4 + 9];
      __m256i bt_c10 = b_t[fco * 4 + 10];
      __m256i bt_c11 = b_t[fco * 4 + 11];

      __m256i p0  = _mm256_madd_epi16(a_r, bt_c0);
      __m256i p1  = _mm256_madd_epi16(a_r, bt_c1);
      __m256i p2  = _mm256_madd_epi16(a_r, bt_c2);
      __m256i p3  = _mm256_madd_epi16(a_r, bt_c3);
      __m256i p8  = _mm256_madd_epi16(a_r, bt_c8);
      __m256i p9  = _mm256_madd_epi16(a_r, bt_c9);
      __m256i p10 = _mm256_madd_epi16(a_r, bt_c10);
      __m256i p11 = _mm256_madd_epi16(a_r, bt_c11);

      // Combine low lanes from P0 and P8, high lanes from them, and the same
      // with P1:P9 and so on
      __m256i p0l = _mm256_permute2x128_si256(p0, p8,  0x20);
      __m256i p0h = _mm256_permute2x128_si256(p0, p8,  0x31);
      __m256i p1l = _mm256_permute2x128_si256(p1, p9,  0x20);
      __m256i p1h = _mm256_permute2x128_si256(p1, p9,  0x31);
      __m256i p2l = _mm256_permute2x128_si256(p2, p10, 0x20);
      __m256i p2h = _mm256_permute2x128_si256(p2, p10, 0x31);
      __m256i p3l = _mm256_permute2x128_si256(p3, p11, 0x20);
      __m256i p3h = _mm256_permute2x128_si256(p3, p11, 0x31);

      __m256i s0  = _mm256_add_epi32(p0l, p0h);
      __m256i s1  = _mm256_add_epi32(p1l, p1h);
      __m256i s2  = _mm256_add_epi32(p2l, p2h);
      __m256i s3  = _mm256_add_epi32(p3l, p3h);

      __m256i s4  = _mm256_unpacklo_epi64(s0, s1);
      __m256i s5  = _mm256_unpackhi_epi64(s0, s1);
      __m256i s6  = _mm256_unpacklo_epi64(s2, s3);
      __m256i s7  = _mm256_unpackhi_epi64(s2, s3);

      __m256i s8  = _mm256_add_epi32(s4, s5);
      __m256i s9  = _mm256_add_epi32(s6, s7);

      __m256i res = _mm256_hadd_epi32(s8, s9);
      results_32[fco] = truncate(res, debias, shift);
    }
    __m256i final_col = _mm256_packs_epi32(results_32[0], results_32[1]);
    _mm256_store_si256((__m256i *)output + y, final_col);
  }
}

// NOTE: The strides measured by s_stride_log2 and d_stride_log2 are in units
// of 16 coeffs, not 1!
static void transpose_16x16_stride(const int16_t *src,
                                         int16_t *dst,
                                         uint8_t  s_stride_log2,
                                         uint8_t  d_stride_log2)
{
  __m256i tmp_128[16];
  for (uint32_t i = 0; i < 16; i += 8) {

    // After every n-bit unpack, 2n-bit units in the vectors will be in
    // correct order. Pair words first, then dwords, then qwords. After that,
    // whole lanes will be correct.
    __m256i tmp_32[8];
    __m256i tmp_64[8];

    __m256i m[8] = {
      _mm256_load_si256((const __m256i *)src + ((i + 0) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 1) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 2) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 3) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 4) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 5) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 6) << s_stride_log2)),
      _mm256_load_si256((const __m256i *)src + ((i + 7) << s_stride_log2)),
    };

    tmp_32[0]      = _mm256_unpacklo_epi16(     m[0],      m[1]);
    tmp_32[1]      = _mm256_unpacklo_epi16(     m[2],      m[3]);
    tmp_32[2]      = _mm256_unpackhi_epi16(     m[0],      m[1]);
    tmp_32[3]      = _mm256_unpackhi_epi16(     m[2],      m[3]);

    tmp_32[4]      = _mm256_unpacklo_epi16(     m[4],      m[5]);
    tmp_32[5]      = _mm256_unpacklo_epi16(     m[6],      m[7]);
    tmp_32[6]      = _mm256_unpackhi_epi16(     m[4],      m[5]);
    tmp_32[7]      = _mm256_unpackhi_epi16(     m[6],      m[7]);


    tmp_64[0]      = _mm256_unpacklo_epi32(tmp_32[0], tmp_32[1]);
    tmp_64[1]      = _mm256_unpacklo_epi32(tmp_32[2], tmp_32[3]);
    tmp_64[2]      = _mm256_unpackhi_epi32(tmp_32[0], tmp_32[1]);
    tmp_64[3]      = _mm256_unpackhi_epi32(tmp_32[2], tmp_32[3]);

    tmp_64[4]      = _mm256_unpacklo_epi32(tmp_32[4], tmp_32[5]);
    tmp_64[5]      = _mm256_unpacklo_epi32(tmp_32[6], tmp_32[7]);
    tmp_64[6]      = _mm256_unpackhi_epi32(tmp_32[4], tmp_32[5]);
    tmp_64[7]      = _mm256_unpackhi_epi32(tmp_32[6], tmp_32[7]);


    tmp_128[i + 0] = _mm256_unpacklo_epi64(tmp_64[0], tmp_64[4]);
    tmp_128[i + 1] = _mm256_unpackhi_epi64(tmp_64[0], tmp_64[4]);
    tmp_128[i + 2] = _mm256_unpacklo_epi64(tmp_64[2], tmp_64[6]);
    tmp_128[i + 3] = _mm256_unpackhi_epi64(tmp_64[2], tmp_64[6]);

    tmp_128[i + 4] = _mm256_unpacklo_epi64(tmp_64[1], tmp_64[5]);
    tmp_128[i + 5] = _mm256_unpackhi_epi64(tmp_64[1], tmp_64[5]);
    tmp_128[i + 6] = _mm256_unpacklo_epi64(tmp_64[3], tmp_64[7]);
    tmp_128[i + 7] = _mm256_unpackhi_epi64(tmp_64[3], tmp_64[7]);
  }

  for (uint32_t i = 0; i < 8; i++) {
    uint32_t loid     = i + 0;
    uint32_t hiid     = i + 8;

    uint32_t dst_loid = loid << d_stride_log2;
    uint32_t dst_hiid = hiid << d_stride_log2;

    __m256i lo       = tmp_128[loid];
    __m256i hi       = tmp_128[hiid];
    __m256i final_lo = _mm256_permute2x128_si256(lo, hi, 0x20);
    __m256i final_hi = _mm256_permute2x128_si256(lo, hi, 0x31);

    _mm256_store_si256((__m256i *)dst + dst_loid, final_lo);
    _mm256_store_si256((__m256i *)dst + dst_hiid, final_hi);
  }
}

static void transpose_16x16(const int16_t *src, int16_t *dst)
{
  transpose_16x16_stride(src, dst, 0, 0);
}

static void transpose_32x32(const int16_t *src, int16_t *dst)
{
  transpose_16x16_stride(src +       0, dst +       0, 1, 1);
  transpose_16x16_stride(src +      16, dst + 16 * 32, 1, 1);
  transpose_16x16_stride(src + 16 * 32, dst +      16, 1, 1);
  transpose_16x16_stride(src + 16 * 33, dst + 16 * 33, 1, 1);
}

static __m256i truncate_inv(__m256i v, int32_t shift)
{
  int32_t add = 1 << (shift - 1);

  __m256i debias  = _mm256_set1_epi32(add);
  __m256i v2      = _mm256_add_epi32 (v,  debias);
  __m256i trunced = _mm256_srai_epi32(v2, shift);
  return  trunced;
}

static __m256i extract_odds(__m256i v)
{
  // 0 1 2 3 4 5 6 7 | 8 9 a b c d e f => 1 3 5 7 1 3 5 7 | 9 b d f 9 b d f
  const __m256i oddmask = _mm256_setr_epi8( 2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15);

  __m256i tmp = _mm256_shuffle_epi8 (v,   oddmask);
  return _mm256_permute4x64_epi64   (tmp, _MM_SHUFFLE(3, 1, 2, 0));
}

static __m256i extract_combine_odds(__m256i v0, __m256i v1)
{
  // 0 1 2 3 4 5 6 7 | 8 9 a b c d e f => 1 3 5 7 1 3 5 7 | 9 b d f 9 b d f
  const __m256i oddmask = _mm256_setr_epi8( 2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15,
                                            2,  3,  6,  7, 10, 11, 14, 15);

  __m256i tmp0 = _mm256_shuffle_epi8(v0,   oddmask);
  __m256i tmp1 = _mm256_shuffle_epi8(v1,   oddmask);

  __m256i tmp2 = _mm256_blend_epi32 (tmp0, tmp1, 0xcc); // 1100 1100

  return _mm256_permute4x64_epi64   (tmp2, _MM_SHUFFLE(3, 1, 2, 0));
}

// Extract items 2, 6, A and E from first four columns of DCT, order them as
// follows:
// D0,2 D0,6 D1,2 D1,6 D1,a D1,e D0,a D0,e | D2,2 D2,6 D3,2 D3,6 D3,a D3,e D2,a D2,e
static __m256i extract_26ae(const __m256i *tdct)
{
  // 02 03 22 23 06 07 26 27 | 0a 0b 2a 2b 02 0f 2e 2f
  // =>
  // 02 06 22 26 02 06 22 26 | 2a 2e 0a 0e 2a 2e 0a 0e
  const __m256i evens_mask = _mm256_setr_epi8( 0,  1,  8,  9,  4,  5, 12, 13,
                                               0,  1,  8,  9,  4,  5, 12, 13,
                                               4,  5, 12, 13,  0,  1,  8,  9,
                                               4,  5, 12, 13,  0,  1,  8,  9);

  __m256i shufd_0 = _mm256_shuffle_epi32(tdct[0], _MM_SHUFFLE(2, 3, 0, 1));
  __m256i shufd_2 = _mm256_shuffle_epi32(tdct[2], _MM_SHUFFLE(2, 3, 0, 1));

  __m256i cmbd_01 = _mm256_blend_epi32(shufd_0, tdct[1], 0xaa); // 1010 1010
  __m256i cmbd_23 = _mm256_blend_epi32(shufd_2, tdct[3], 0xaa); // 1010 1010

  __m256i evens_01 = _mm256_shuffle_epi8(cmbd_01, evens_mask);
  __m256i evens_23 = _mm256_shuffle_epi8(cmbd_23, evens_mask);

  __m256i evens_0123 = _mm256_unpacklo_epi64(evens_01, evens_23);

  return _mm256_permute4x64_epi64(evens_0123, _MM_SHUFFLE(3, 1, 2, 0));
}

// 2 6 2 6 a e a e | 2 6 2 6 a e a e
static __m256i extract_26ae_vec(__m256i col)
{
  const __m256i mask_26ae = _mm256_set1_epi32(0x0d0c0504);

  // 2 6 2 6 2 6 2 6 | a e a e a e a e
  __m256i reord = _mm256_shuffle_epi8     (col,   mask_26ae);
  __m256i final = _mm256_permute4x64_epi64(reord, _MM_SHUFFLE(3, 1, 2, 0));
  return  final;
}

// D00 D80 D01 D81 D41 Dc1 D40 Dc0 | D40 Dc0 D41 Dc1 D01 D81 D00 D80
static __m256i extract_d048c(const __m256i *tdct)
{
  const __m256i final_shuf = _mm256_setr_epi8( 0,  1,  8,  9,  2,  3, 10, 11,
                                               6,  7, 14, 15,  4,  5, 12, 13,
                                               4,  5, 12, 13,  6,  7, 14, 15,
                                               2,  3, 10, 11,  0,  1,  8,  9);
  __m256i c0 = tdct[0];
  __m256i c1 = tdct[1];

  __m256i c1_2  = _mm256_slli_epi32       (c1,    16);
  __m256i cmbd  = _mm256_blend_epi16      (c0,    c1_2, 0x22); // 0010 0010
  __m256i cmbd2 = _mm256_shuffle_epi32    (cmbd,  _MM_SHUFFLE(2, 0, 2, 0));
  __m256i cmbd3 = _mm256_permute4x64_epi64(cmbd2, _MM_SHUFFLE(3, 1, 2, 0));
  __m256i final = _mm256_shuffle_epi8     (cmbd3, final_shuf);

  return final;
}

// 0 8 0 8 4 c 4 c | 4 c 4 c 0 8 0 8
static __m256i extract_d048c_vec(__m256i col)
{
  const __m256i shufmask = _mm256_setr_epi8( 0,  1,  0,  1,  8,  9,  8,  9,
                                             8,  9,  8,  9,  0,  1,  0,  1,
                                             0,  1,  0,  1,  8,  9,  8,  9,
                                             8,  9,  8,  9,  0,  1,  0,  1);

  __m256i col_db4s = _mm256_shuffle_epi8     (col, shufmask);
  __m256i col_los  = _mm256_permute4x64_epi64(col_db4s, _MM_SHUFFLE(1, 1, 0, 0));
  __m256i col_his  = _mm256_permute4x64_epi64(col_db4s, _MM_SHUFFLE(3, 3, 2, 2));

  __m256i final    = _mm256_unpacklo_epi16   (col_los,  col_his);
  return final;
}

static void partial_butterfly_inverse_16_avx2(const int16_t *src, int16_t *dst, int32_t shift)
{
  __m256i tsrc[16];

  const uint32_t width = 16;

  const int16_t *tdct = &kvz_g_dct_16_t[0][0];

  const __m256i  eo_signmask = _mm256_setr_epi32( 1,  1,  1,  1, -1, -1, -1, -1);
  const __m256i eeo_signmask = _mm256_setr_epi32( 1,  1, -1, -1, -1, -1,  1,  1);
  const __m256i   o_signmask = _mm256_set1_epi32(-1);

  const __m256i final_shufmask = _mm256_setr_epi8( 0,  1,  2,  3,  4,  5,  6,  7,
                                                   8,  9, 10, 11, 12, 13, 14, 15,
                                                   6,  7,  4,  5,  2,  3,  0,  1,
                                                  14, 15, 12, 13, 10, 11,  8,  9);
  transpose_16x16(src, (int16_t *)tsrc);

  const __m256i dct_cols[8] = {
    _mm256_load_si256((const __m256i *)tdct + 0),
    _mm256_load_si256((const __m256i *)tdct + 1),
    _mm256_load_si256((const __m256i *)tdct + 2),
    _mm256_load_si256((const __m256i *)tdct + 3),
    _mm256_load_si256((const __m256i *)tdct + 4),
    _mm256_load_si256((const __m256i *)tdct + 5),
    _mm256_load_si256((const __m256i *)tdct + 6),
    _mm256_load_si256((const __m256i *)tdct + 7),
  };

  // These contain: D1,0 D3,0 D5,0 D7,0 D9,0 Db,0 Dd,0 Df,0 | D1,4 D3,4 D5,4 D7,4 D9,4 Db,4 Dd,4 Df,4
  //                D1,1 D3,1 D5,1 D7,1 D9,1 Db,1 Dd,1 Df,1 | D1,5 D3,5 D5,5 D7,5 D9,5 Db,5 Dd,5 Df,5
  //                D1,2 D3,2 D5,2 D7,2 D9,2 Db,2 Dd,2 Df,2 | D1,6 D3,6 D5,6 D7,6 D9,6 Db,6 Dd,6 Df,6
  //                D1,3 D3,3 D5,3 D7,3 D9,3 Db,3 Dd,3 Df,3 | D1,7 D3,7 D5,7 D7,7 D9,7 Db,7 Dd,7 Df,7
  __m256i dct_col_odds[4];
  for (uint32_t j = 0; j < 4; j++) {
    dct_col_odds[j] = extract_combine_odds(dct_cols[j + 0], dct_cols[j + 4]);
  }
  for (uint32_t j = 0; j < width; j++) {
    __m256i col = tsrc[j];
    __m256i odds = extract_odds(col);

    __m256i o04   = _mm256_madd_epi16           (odds,     dct_col_odds[0]);
    __m256i o15   = _mm256_madd_epi16           (odds,     dct_col_odds[1]);
    __m256i o26   = _mm256_madd_epi16           (odds,     dct_col_odds[2]);
    __m256i o37   = _mm256_madd_epi16           (odds,     dct_col_odds[3]);

    __m256i o0145 = _mm256_hadd_epi32           (o04,      o15);
    __m256i o2367 = _mm256_hadd_epi32           (o26,      o37);

    __m256i o     = _mm256_hadd_epi32           (o0145,    o2367);

    // D0,2 D0,6 D1,2 D1,6 D1,a D1,e D0,a D0,e | D2,2 D2,6 D3,2 D3,6 D3,a D3,e D2,a D2,e
    __m256i d_db2 = extract_26ae(dct_cols);

    // 2 6 2 6 a e a e | 2 6 2 6 a e a e
    __m256i t_db2 = extract_26ae_vec            (col);

    __m256i eo_parts  = _mm256_madd_epi16       (d_db2,    t_db2);
    __m256i eo_parts2 = _mm256_shuffle_epi32    (eo_parts, _MM_SHUFFLE(0, 1, 2, 3));

    // EO0 EO1 EO1 EO0 | EO2 EO3 EO3 EO2
    __m256i eo        = _mm256_add_epi32        (eo_parts, eo_parts2);
    __m256i eo2       = _mm256_permute4x64_epi64(eo,       _MM_SHUFFLE(1, 3, 2, 0));
    __m256i eo3       = _mm256_sign_epi32       (eo2,      eo_signmask);

    __m256i d_db4     = extract_d048c           (dct_cols);
    __m256i t_db4     = extract_d048c_vec       (col);
    __m256i eee_eeo   = _mm256_madd_epi16       (d_db4,   t_db4);

    __m256i eee_eee   = _mm256_permute4x64_epi64(eee_eeo,  _MM_SHUFFLE(3, 0, 3, 0));
    __m256i eeo_eeo1  = _mm256_permute4x64_epi64(eee_eeo,  _MM_SHUFFLE(1, 2, 1, 2));

    __m256i eeo_eeo2  = _mm256_sign_epi32       (eeo_eeo1, eeo_signmask);

    // EE0 EE1 EE2 EE3 | EE3 EE2 EE1 EE0
    __m256i ee        = _mm256_add_epi32        (eee_eee,  eeo_eeo2);
    __m256i e         = _mm256_add_epi32        (ee,       eo3);

    __m256i o_neg     = _mm256_sign_epi32       (o,        o_signmask);
    __m256i o_lo      = _mm256_blend_epi32      (o,        o_neg, 0xf0); // 1111 0000
    __m256i o_hi      = _mm256_blend_epi32      (o,        o_neg, 0x0f); // 0000 1111

    __m256i res_lo    = _mm256_add_epi32        (e,        o_lo);
    __m256i res_hi    = _mm256_add_epi32        (e,        o_hi);
    __m256i res_hi2   = _mm256_permute4x64_epi64(res_hi,   _MM_SHUFFLE(1, 0, 3, 2));

    __m256i res_lo_t  = truncate_inv(res_lo,  shift);
    __m256i res_hi_t  = truncate_inv(res_hi2, shift);

    __m256i res_16_1  = _mm256_packs_epi32      (res_lo_t, res_hi_t);
    __m256i final     = _mm256_shuffle_epi8     (res_16_1, final_shufmask);

    _mm256_store_si256((__m256i *)dst + j, final);
  }
}

static void matrix_idct_16x16_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (bitdepth - 8);
  ALIGNED(64) int16_t tmp[16 * 16];

  partial_butterfly_inverse_16_avx2(input, tmp,    shift_1st);
  partial_butterfly_inverse_16_avx2(tmp,   output, shift_2nd);
}

static void matrix_dct_16x16_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[16] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[16] + 8;

  const int16_t *dct  = &kvz_g_dct_16[0][0];

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

  __m256i tmpres[16];
  matmul_16x16_a_bt_t(input,  dct, tmpres, shift_1st);
  matmul_16x16_a_bt  (dct, tmpres, output, shift_2nd);
}

static __m256i get_overflows(const __m256i a, const __m256i b, const __m256i res, const __m256i of_adjust_mask)
{
  const __m256i ones = _mm256_set1_epi16(1);

  __m256i src_signdiff = _mm256_xor_si256  (a,            b);
  __m256i a_r_signdiff = _mm256_xor_si256  (a,            res);

  __m256i of_possible  = _mm256_xor_si256  (src_signdiff, of_adjust_mask);

  __m256i overflows    = _mm256_and_si256  (of_possible,  a_r_signdiff);
          overflows    = _mm256_srai_epi16 (overflows,    15);

  __m256i of_signs     = _mm256_srai_epi16 (a,            15);
          of_signs     = _mm256_or_si256   (of_signs,     ones);
  return                 _mm256_and_si256  (overflows,    of_signs);
}

/*
 * You need more than 16 bits to store the result of signed 16b-16b operation,
 * this one stores the low 16b in lo and high 16b in hi. The high 16 bits can
 * only be -1, 0 or 1.
 *
 * of_possible_mask is either all zero bits for subtraction, or all ones for
 * addition
 */
static void sub_16_16_hilo(const __m256i a, const __m256i b, __m256i *lo, __m256i *hi)
{
  const __m256i zero = _mm256_setzero_si256();

  *lo = _mm256_sub_epi16(a, b);
  *hi = get_overflows(a, b, *lo, zero);
}

static void add_16_16_hilo(const __m256i a, const __m256i b, __m256i *lo, __m256i *hi)
{
  const __m256i ff     = _mm256_set1_epi8(-1);

  *lo = _mm256_add_epi16(a, b);
  *hi = get_overflows(a, b, *lo, ff);
}

static __m256i reverse_16x16b_in_lanes(const __m256i v)
{
  const __m256i lanerev = _mm256_setr_epi16(0x0f0e, 0x0d0c, 0x0b0a, 0x0908,
                                            0x0706, 0x0504, 0x0302, 0x0100,
                                            0x0f0e, 0x0d0c, 0x0b0a, 0x0908,
                                            0x0706, 0x0504, 0x0302, 0x0100);
  return _mm256_shuffle_epi8(v, lanerev);
}

static __m256i reverse_16x16b(const __m256i v)
{
  __m256i tmp = reverse_16x16b_in_lanes(v);
  return        _mm256_permute4x64_epi64(tmp, _MM_SHUFFLE(1, 0, 3, 2));
}

static __m256i m256_from_2xm128(const __m128i lo, const __m128i hi)
{
  __m256i result = _mm256_castsi128_si256 (lo);
  return           _mm256_inserti128_si256(result, hi, 1);
}

// Get a vector consisting of the divisible-by-8 coeffs in DCT's first two
// columns, in order:
// DC00 DC00 DC00 DC00 DC10 DC10 DC10 DC10 | DC08 DC08 DC08 DC08 DC18 DC18 DC18 DC18
static __m256i get_dct_db8_vec(const int16_t *dct_t, uint32_t offset)
{
  const __m256i reorder_mask = _mm256_setr_epi32(0x01000100, 0x01000100, 0x05040504, 0x05040504,
                                                 0x03020302, 0x03020302, 0x07060706, 0x07060706);

  uint16_t coeff00 = (uint16_t)dct_t[offset +  0];
  uint16_t coeff08 = (uint16_t)dct_t[offset +  8];
  uint16_t coeff10 = (uint16_t)dct_t[offset + 32];
  uint16_t coeff18 = (uint16_t)dct_t[offset + 40];

  uint64_t col_db8_packed = (((uint64_t)coeff00) <<  0) |
                            (((uint64_t)coeff08) << 16) |
                            (((uint64_t)coeff10) << 32) |
                            (((uint64_t)coeff18) << 48);

  __m256i col_db8         = _mm256_set1_epi64x (col_db8_packed);
          col_db8         = _mm256_shuffle_epi8(col_db8, reorder_mask);
  __m256i col_db8_shifted = _mm256_slli_epi16  (col_db8, 8);

  return                    _mm256_blend_epi32 (col_db8, col_db8_shifted, 0xaa);
}

// Get first 4 DCT coeffs from DCT rows (rowoff + 4) and (rowoff + 12),
// shifting copies of the coeffs 8 bits left to do half of the 65536-factor
// multiplication required for high words of EEO coeffs
static __m256i get_dct_db4_vec(const int16_t *dct, uint32_t rowoff)
{
  ALIGNED(32) uint64_t buf[4];

  uint64_t r4_0123 = *(uint64_t *)(dct + (rowoff + 0x04) * 32);
  uint64_t rc_0123 = *(uint64_t *)(dct + (rowoff + 0x0c) * 32);

  buf[0] = r4_0123;
  buf[1] = r4_0123 << 8;
  buf[2] = rc_0123;
  buf[3] = rc_0123 << 8;

  return _mm256_load_si256((const __m256i *)buf);
}

// Get first 8 coeffs from DCT rows (rowoff + 2) and (rowoff + 10) and return
// them in a single YMM
static __m256i get_dct_db2_vec(const int16_t *dct, uint32_t rowoff)
{
  __m128i row2 = _mm_load_si128((const __m128i *)(dct + (rowoff +  2) * 32));
  __m128i rowa = _mm_load_si128((const __m128i *)(dct + (rowoff + 10) * 32));

  return m256_from_2xm128(row2, rowa);
}

static void partial_butterfly_32_avx2(const int16_t *src, int16_t *dst, int32_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  const int16_t *dct   = (const int16_t *)(&kvz_g_dct_32  [0][0]);
  const int16_t *dct_t = (const int16_t *)(&kvz_g_dct_32_t[0][0]);

  const __m256i ff           = _mm256_set1_epi32      (-1);
  const __m256i ones         = _mm256_set1_epi16      ( 1);
  const __m128i ff_128       = _mm256_castsi256_si128 (ff);
  const __m256i lolane_smask = _mm256_inserti128_si256(ones, ff_128, 0);
  const __m256i hilane_mask  = _mm256_cmpeq_epi16     (ones, lolane_smask);

  const __m256i eee_lohi_shuf   = _mm256_setr_epi32(0x01000706, 0x09080f0e, 0x03020504, 0x0b0a0d0c,
                                                    0x01000706, 0x09080f0e, 0x03020504, 0x0b0a0d0c);

  const __m256i eee_eee_sgnmask = _mm256_setr_epi32(0x00010001, 0x00010001, 0x00010001, 0x00010001,
                                                    0xffff0001, 0xffff0001, 0xffff0001, 0xffff0001);

  const __m256i dct_c0_db8 = get_dct_db8_vec(dct_t, 0);
  const __m256i dct_c1_db8 = get_dct_db8_vec(dct_t, 16);

  const __m256i dct_r4c_db4[2] = {
    get_dct_db4_vec(dct, 0x00),
    get_dct_db4_vec(dct, 0x10),
  };

  const __m256i dct_r_db2[4] = {
    get_dct_db2_vec(dct, 0),
    get_dct_db2_vec(dct, 4),
    get_dct_db2_vec(dct, 16),
    get_dct_db2_vec(dct, 20),
  };

  __m256i res_tp[2 * 32];
  for (uint32_t i = 0; i < 32; i++) {
    __m256i lo = _mm256_load_si256((const __m256i *)src + 2 * i + 0);
    __m256i hi = _mm256_load_si256((const __m256i *)src + 2 * i + 1);

    __m256i hi_rev = reverse_16x16b(hi);

    __m256i e_lo, e_hi, o_lo, o_hi;
    add_16_16_hilo(lo, hi_rev, &e_lo, &e_hi);
    sub_16_16_hilo(lo, hi_rev, &o_lo, &o_hi);

    __m256i erev_lo   = reverse_16x16b(e_lo);
    __m256i erev_hi   = reverse_16x16b(e_hi);

    // Hack! Negate low lanes to do subtractions there, but retain non-negated
    // low lanes for overflow detection because the original value is
    // essentially what was subtracted there by adding its complement
    __m256i erev_lo_n = _mm256_sign_epi16(erev_lo, lolane_smask);
    __m256i erev_hi_n = _mm256_sign_epi16(erev_hi, lolane_smask);

    // eo0 eo1 eo2 eo3 eo4 eo5 eo6 eo7 | ee7 ee6 ee5 ee4 ee3 ee2 ee1 ee0
    __m256i eo_ee_lo  = _mm256_add_epi16 (e_lo,    erev_lo_n);
    __m256i eo_ee_hi  = _mm256_add_epi16 (e_hi,    erev_hi_n);

    __m256i eo_ee_hi2 = get_overflows    (e_lo, erev_lo, eo_ee_lo, hilane_mask);
            eo_ee_hi  = _mm256_add_epi16 (eo_ee_hi, eo_ee_hi2);

    __m256i eo_eo_lo       = _mm256_permute4x64_epi64(eo_ee_lo, _MM_SHUFFLE(1, 0, 1, 0));
    __m256i eo_eo_hi       = _mm256_permute4x64_epi64(eo_ee_hi, _MM_SHUFFLE(1, 0, 1, 0));

    // Rev: ee7 ee6 ee5 ee4 ee7 ee6 ee5 ee4 | ee3 ee2 ee1 ee0 ee3 ee2 ee1 ee0
    // Fwd: ee0 ee1 ee2 ee3 ee0 ee1 ee2 ee3 | ee4 ee5 ee6 ee7 ee4 ee5 ee6 ee7
    __m256i ee_ee_rev_lo   = _mm256_permute4x64_epi64(eo_ee_lo, _MM_SHUFFLE(3, 3, 2, 2));
    __m256i ee_ee_lo       = _mm256_permute4x64_epi64(eo_ee_lo, _MM_SHUFFLE(2, 2, 3, 3));
            ee_ee_lo       = reverse_16x16b_in_lanes (ee_ee_lo);

    __m256i ee_ee_rev_hi   = _mm256_permute4x64_epi64(eo_ee_hi, _MM_SHUFFLE(3, 3, 2, 2));
    __m256i ee_ee_hi       = _mm256_permute4x64_epi64(eo_ee_hi, _MM_SHUFFLE(2, 2, 3, 3));
            ee_ee_hi       = reverse_16x16b_in_lanes (ee_ee_hi);

    __m256i ee_ee_rev_lo_n = _mm256_sign_epi16(ee_ee_rev_lo, lolane_smask);
    __m256i ee_ee_rev_hi_n = _mm256_sign_epi16(ee_ee_rev_hi, lolane_smask);

    // eeo0 eeo1 eeo2 eeo3 eeo0 eeo1 eeo2 eeo3 | eee3 eee2 eee1 eee0 eee3 eee2 eee1 eee0
    __m256i eeo_eee_lo     = _mm256_add_epi16 (ee_ee_lo, ee_ee_rev_lo_n);
    __m256i eeo_eee_hi     = _mm256_add_epi16 (ee_ee_hi, ee_ee_rev_hi_n);

    __m256i eeo_eee_hi2    = get_overflows(ee_ee_lo, ee_ee_rev_lo, eeo_eee_lo, hilane_mask);
            eeo_eee_hi     = _mm256_add_epi16(eeo_eee_hi, eeo_eee_hi2);

    // Multiply these guys by 256 and also the corresponding DCT coefficients
    // by 256, to multiply their product by 65536. Neither these values nor
    // the coeffs will exceed 255, so this is overflow safe.
    __m256i eeo_eee_hi_shed = _mm256_slli_epi16(eeo_eee_hi, 8);

    // Discard eeo's (low lane), duplicate the high lane, and reorder eee's
    // in-lane. Finally invert the eee2/eee3 components in the high lane.
    __m256i eee_lohi       = _mm256_blend_epi32(eeo_eee_lo, eeo_eee_hi_shed, 0xc0);
            eee_lohi       = _mm256_permute4x64_epi64(eee_lohi, _MM_SHUFFLE(3, 2, 3, 2));

    __m256i eee_lh_ordered = _mm256_shuffle_epi8(eee_lohi,       eee_lohi_shuf);
    __m256i eee_lh_final   = _mm256_sign_epi16  (eee_lh_ordered, eee_eee_sgnmask);

    // D00, D08, D10 and D18 in four parts fit for hadding
    __m256i d00_d08 = _mm256_madd_epi16(eee_lh_final, dct_c0_db8);
    __m256i d10_d18 = _mm256_madd_epi16(eee_lh_final, dct_c1_db8);

    __m256i eeo_lohi = _mm256_unpacklo_epi64   (eeo_eee_lo, eeo_eee_hi_shed);
            eeo_lohi = _mm256_permute4x64_epi64(eeo_lohi,   _MM_SHUFFLE(1, 0, 1, 0));

    __m256i d_db4[2] = {
      _mm256_madd_epi16(eeo_lohi, dct_r4c_db4[0]),
      _mm256_madd_epi16(eeo_lohi, dct_r4c_db4[1]),
    };

    __m256i db2_parts[4];
    for (uint32_t j = 0; j < 4; j++) {
      const __m256i dr2a = dct_r_db2[j];
            __m256i lo   = _mm256_madd_epi16(dr2a, eo_eo_lo);
            __m256i hi   = _mm256_madd_epi16(dr2a, eo_eo_hi);
                    hi   = _mm256_slli_epi32(hi,   16);

      db2_parts[j]       = _mm256_add_epi32 (lo,   hi);
    }

    __m256i odd_parts[16];
    for (uint32_t j = 0; j < 16; j++) {
      __m256i drow_lo = _mm256_load_si256((const __m256i *)(dct + (j * 2 + 1) * 32));
      __m256i odds_lo = _mm256_madd_epi16(o_lo,    drow_lo);
      __m256i odds_hi = _mm256_madd_epi16(o_hi,    drow_lo);
              odds_hi = _mm256_slli_epi32(odds_hi, 16);

      odd_parts[j]    = _mm256_add_epi32 (odds_lo, odds_hi);
    }

    // Rearrange odds so that parts belonging to any single one are all inside
    // one lane - combine 01 | 09 ; 03 | 0b ; 05 | 0d ; 07 | 0f
    //                    11 | 19 ; 13 | 1b ; 15 | 1d ; 17 | 1f
    __m256i odd_parts2[8];
    for (uint32_t j = 0; j < 8; j++) {

      // Turn 0, 1, 2, 3, 4, 5, 6, 7 into:
      //      0, 1, 2, 3, 8, 9, a, b
      uint32_t j_lo  = j & 0x03;
      uint32_t j_hi  = j & 0x04;
      uint32_t id_lo = j_lo  | (j_hi << 1);
      uint32_t id_hi = id_lo | 4;

      __m256i odd_lo = _mm256_permute2x128_si256(odd_parts[id_lo],
                                                 odd_parts[id_hi],
                                                 0x20);

      __m256i odd_hi = _mm256_permute2x128_si256(odd_parts[id_lo],
                                                 odd_parts[id_hi],
                                                 0x31);

      odd_parts2[j]  = _mm256_add_epi32(odd_lo, odd_hi);
    }

    // First stage HADDs...
    __m256i d0001_0809    = _mm256_hadd_epi32(d00_d08,      odd_parts2[0]);
    __m256i d1011_1819    = _mm256_hadd_epi32(d10_d18,      odd_parts2[4]);

    __m256i d0405_0c0d    = _mm256_hadd_epi32(d_db4[0],     odd_parts2[2]);
    __m256i d1415_1c1d    = _mm256_hadd_epi32(d_db4[1],     odd_parts2[6]);

    __m256i d0203_0a0b    = _mm256_hadd_epi32(db2_parts[0], odd_parts2[1]);
    __m256i d0607_0e0f    = _mm256_hadd_epi32(db2_parts[1], odd_parts2[3]);

    __m256i d1213_1a1b    = _mm256_hadd_epi32(db2_parts[2], odd_parts2[5]);
    __m256i d1617_1e1f    = _mm256_hadd_epi32(db2_parts[3], odd_parts2[7]);

    // .. and second stage
    __m256i d0123_89ab_lo = _mm256_hadd_epi32(d0001_0809,   d0203_0a0b);
    __m256i d0123_89ab_hi = _mm256_hadd_epi32(d1011_1819,   d1213_1a1b);

    __m256i d4567_cdef_lo = _mm256_hadd_epi32(d0405_0c0d,   d0607_0e0f);
    __m256i d4567_cdef_hi = _mm256_hadd_epi32(d1415_1c1d,   d1617_1e1f);

            d0123_89ab_lo = truncate(d0123_89ab_lo, debias, shift);
            d0123_89ab_hi = truncate(d0123_89ab_hi, debias, shift);

            d4567_cdef_lo = truncate(d4567_cdef_lo, debias, shift);
            d4567_cdef_hi = truncate(d4567_cdef_hi, debias, shift);

    __m256i final_lo      = _mm256_packs_epi32(d0123_89ab_lo, d4567_cdef_lo);
    __m256i final_hi      = _mm256_packs_epi32(d0123_89ab_hi, d4567_cdef_hi);

    _mm256_store_si256(res_tp + (i * 2) + 0, final_lo);
    _mm256_store_si256(res_tp + (i * 2) + 1, final_hi);
  }
  transpose_32x32((const int16_t *)res_tp, dst);
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

static void matrix_dct_32x32_avx2(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int32_t shift_1st = kvz_g_convert_to_bit[32] + 1 + (bitdepth - 8);
  int32_t shift_2nd = kvz_g_convert_to_bit[32] + 8;
  ALIGNED(64) int16_t tmp[32 * 32];

  partial_butterfly_32_avx2(input, tmp,    shift_1st);
  partial_butterfly_32_avx2(tmp,   output, shift_2nd);
}

// Macro that generates 2D transform functions with clipping values.
// Sets correct shift values and matrices according to transform type and
// block size. Performs matrix multiplication horizontally and vertically.
#define TRANSFORM(type, n) static void matrix_ ## type ## _ ## n ## x ## n ## _avx2(int8_t bitdepth, const int16_t *input, int16_t *output)\
{\
  int32_t shift_1st = kvz_g_convert_to_bit[n] + 1 + (bitdepth - 8); \
  int32_t shift_2nd = kvz_g_convert_to_bit[n] + 8; \
  ALIGNED(64) int16_t tmp[n * n];\
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
  ALIGNED(64) int16_t tmp[n * n];\
  const int16_t *tdct = &kvz_g_ ## type ## _ ## n ## _t[0][0];\
  const int16_t *dct = &kvz_g_ ## type ## _ ## n [0][0];\
\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tdct, input, tmp, shift_1st);\
  mul_clip_matrix_ ## n ## x ## n ## _avx2(tmp, dct, output, shift_2nd);\
}\

// Ha, we've got a tailored implementation for these
// TRANSFORM(dst, 4);
// ITRANSFORM(dst, 4);
// TRANSFORM(dct, 4);
// ITRANSFORM(dct, 4);
// TRANSFORM(dct, 8);
// ITRANSFORM(dct, 8);
// TRANSFORM(dct, 16);
// ITRANSFORM(dct, 16);

// TRANSFORM(dct, 32);

// Generate all the transform functions

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
