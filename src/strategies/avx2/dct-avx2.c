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

// 16x16 matrix multiplication with value clipping.
// Parameters: Two 16x16 matrices containing 16-bit values in consecutive addresses,
//             destination for the result and the shift value for clipping.
static void mul_clip_matrix_16x16_avx2(const int16_t *left, const int16_t *right, int16_t *dst, const int32_t shift)
{
  const int32_t add    = 1 << (shift - 1);
  const __m256i debias = _mm256_set1_epi32(add);

  __m256i sliced_right[16];
  for (int32_t dry = 0; dry < 16; dry += 2) {
      __m256i right_up        = _mm256_load_si256((const __m256i *)right + dry + 0);
      __m256i right_dn        = _mm256_load_si256((const __m256i *)right + dry + 1);

      __m256i right_slices_lo = _mm256_unpacklo_epi16(right_up, right_dn);
      __m256i right_slices_hi = _mm256_unpackhi_epi16(right_up, right_dn);

      sliced_right[dry + 0] = right_slices_lo;
      sliced_right[dry + 1] = right_slices_hi;
  }
  for (int32_t dry = 0; dry < 16; dry += 2) {
    __m256i accum1 = _mm256_setzero_si256();
    __m256i accum2 = _mm256_setzero_si256();
    __m256i accum3 = _mm256_setzero_si256();
    __m256i accum4 = _mm256_setzero_si256();

    for (int32_t lx = 0; lx < 16; lx += 2) {
      const int32_t *curr_left_up = (const int32_t *)(left + (dry + 0) * 16 + lx);
      const int32_t *curr_left_dn = (const int32_t *)(left + (dry + 1) * 16 + lx);

      __m256i left_slice_lo   = _mm256_set1_epi32(*curr_left_up);
      __m256i left_slice_hi   = _mm256_set1_epi32(*curr_left_dn);

      __m256i right_slices_lo = sliced_right[lx + 0];
      __m256i right_slices_hi = sliced_right[lx + 1];

      __m256i prod1 = _mm256_madd_epi16(left_slice_lo, right_slices_lo);
      __m256i prod2 = _mm256_madd_epi16(left_slice_hi, right_slices_lo);
      __m256i prod3 = _mm256_madd_epi16(left_slice_lo, right_slices_hi);
      __m256i prod4 = _mm256_madd_epi16(left_slice_hi, right_slices_hi);

      accum1 = _mm256_add_epi32(accum1, prod1);
      accum2 = _mm256_add_epi32(accum2, prod2);
      accum3 = _mm256_add_epi32(accum3, prod3);
      accum4 = _mm256_add_epi32(accum4, prod4);
    }
    __m256i accum1_tr = truncate(accum1, debias, shift);
    __m256i accum2_tr = truncate(accum2, debias, shift);
    __m256i accum3_tr = truncate(accum3, debias, shift);
    __m256i accum4_tr = truncate(accum4, debias, shift);

    __m256i out_up = _mm256_packs_epi32(accum1_tr, accum3_tr);
    __m256i out_dn = _mm256_packs_epi32(accum2_tr, accum4_tr);

    _mm256_store_si256((__m256i *)dst + dry + 0, out_up);
    _mm256_store_si256((__m256i *)dst + dry + 1, out_dn);
  }
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

// Generate all the transform functions
TRANSFORM(dct, 32);

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
