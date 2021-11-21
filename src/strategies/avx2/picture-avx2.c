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

/*
 * \file
 */

#include "global.h"

#if COMPILE_INTEL_AVX2
#include "kvazaar.h"
#if KVZ_BIT_DEPTH == 8
#include "strategies/avx2/picture-avx2.h"
#include "strategies/avx2/reg_sad_pow2_widths-avx2.h"

#include <immintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <string.h>
#include "strategies/strategies-picture.h"
#include "strategyselector.h"
#include "strategies/generic/picture-generic.h"

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
uint32_t kvz_reg_sad_avx2(const uint8_t * const data1, const uint8_t * const data2,
                          const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  if (width == 0)
    return 0;
  if (width == 4)
    return reg_sad_w4(data1, data2, height, stride1, stride2);
  if (width == 8)
    return reg_sad_w8(data1, data2, height, stride1, stride2);
  if (width == 12)
    return reg_sad_w12(data1, data2, height, stride1, stride2);
  if (width == 16)
    return reg_sad_w16(data1, data2, height, stride1, stride2);
  if (width == 24)
    return reg_sad_w24(data1, data2, height, stride1, stride2);
  if (width == 32)
    return reg_sad_w32(data1, data2, height, stride1, stride2);
  if (width == 64)
    return reg_sad_w64(data1, data2, height, stride1, stride2);
  else
    return reg_sad_arbitrary(data1, data2, width, height, stride1, stride2);
}

/**
* \brief Calculate SAD for 8x8 bytes in continuous memory.
*/
static INLINE __m256i inline_8bit_sad_8x8_avx2(const __m256i *const a, const __m256i *const b)
{
  __m256i sum0, sum1;
  sum0 = _mm256_sad_epu8(_mm256_load_si256(a + 0), _mm256_load_si256(b + 0));
  sum1 = _mm256_sad_epu8(_mm256_load_si256(a + 1), _mm256_load_si256(b + 1));

  return _mm256_add_epi32(sum0, sum1);
}


/**
* \brief Calculate SAD for 16x16 bytes in continuous memory.
*/
static INLINE __m256i inline_8bit_sad_16x16_avx2(const __m256i *const a, const __m256i *const b)
{
  const unsigned size_of_8x8 = 8 * 8 / sizeof(__m256i);

  // Calculate in 4 chunks of 16x4.
  __m256i sum0, sum1, sum2, sum3;
  sum0 = inline_8bit_sad_8x8_avx2(a + 0 * size_of_8x8, b + 0 * size_of_8x8);
  sum1 = inline_8bit_sad_8x8_avx2(a + 1 * size_of_8x8, b + 1 * size_of_8x8);
  sum2 = inline_8bit_sad_8x8_avx2(a + 2 * size_of_8x8, b + 2 * size_of_8x8);
  sum3 = inline_8bit_sad_8x8_avx2(a + 3 * size_of_8x8, b + 3 * size_of_8x8);

  sum0 = _mm256_add_epi32(sum0, sum1);
  sum2 = _mm256_add_epi32(sum2, sum3);

  return _mm256_add_epi32(sum0, sum2);
}


/**
* \brief Get sum of the low 32 bits of four 64 bit numbers from __m256i as uint32_t.
*/
static INLINE uint32_t m256i_horizontal_sum(const __m256i sum)
{
  // Add the high 128 bits to low 128 bits.
  __m128i mm128_result = _mm_add_epi32(_mm256_castsi256_si128(sum), _mm256_extractf128_si256(sum, 1));
  // Add the high 64 bits  to low 64 bits.
  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, mm128_result);
  return result[0] + result[2];
}


static unsigned sad_8bit_8x8_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_8x8_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_16x16_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_16x16_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_32x32_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;

  const unsigned size_of_8x8 = 8 * 8 / sizeof(__m256i);
  const unsigned size_of_32x32 = 32 * 32 / sizeof(__m256i);

  // Looping 512 bytes at a time seems faster than letting VC figure it out
  // through inlining, like inline_8bit_sad_16x16_avx2 does.
  __m256i sum0 = inline_8bit_sad_8x8_avx2(a, b);
  for (unsigned i = size_of_8x8; i < size_of_32x32; i += size_of_8x8) {
    __m256i sum1 = inline_8bit_sad_8x8_avx2(a + i, b + i);
    sum0 = _mm256_add_epi32(sum0, sum1);
  }

  return m256i_horizontal_sum(sum0);
}


static unsigned sad_8bit_64x64_avx2(const uint8_t * buf1, const uint8_t * buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;

  const unsigned size_of_8x8 = 8 * 8 / sizeof(__m256i);
  const unsigned size_of_64x64 = 64 * 64 / sizeof(__m256i);

  // Looping 512 bytes at a time seems faster than letting VC figure it out
  // through inlining, like inline_8bit_sad_16x16_avx2 does.
  __m256i sum0 = inline_8bit_sad_8x8_avx2(a, b);
  for (unsigned i = size_of_8x8; i < size_of_64x64; i += size_of_8x8) {
    __m256i sum1 = inline_8bit_sad_8x8_avx2(a + i, b + i);
    sum0 = _mm256_add_epi32(sum0, sum1);
  }

  return m256i_horizontal_sum(sum0);
}

static unsigned satd_4x4_8bit_avx2(const uint8_t *org, const uint8_t *cur)
{

  __m128i original = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)org));
  __m128i current = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)cur));

  __m128i diff_lo = _mm_sub_epi16(current, original);

  original = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(org + 8)));
  current = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(cur + 8)));

  __m128i diff_hi = _mm_sub_epi16(current, original);


  //Hor
  __m128i row0 = _mm_hadd_epi16(diff_lo, diff_hi);
  __m128i row1 = _mm_hsub_epi16(diff_lo, diff_hi);

  __m128i row2 = _mm_hadd_epi16(row0, row1);
  __m128i row3 = _mm_hsub_epi16(row0, row1);

  //Ver
  row0 = _mm_hadd_epi16(row2, row3);
  row1 = _mm_hsub_epi16(row2, row3);

  row2 = _mm_hadd_epi16(row0, row1);
  row3 = _mm_hsub_epi16(row0, row1);

  //Abs and sum
  row2 = _mm_abs_epi16(row2);
  row3 = _mm_abs_epi16(row3);

  row3 = _mm_add_epi16(row2, row3);

  row3 = _mm_add_epi16(row3, _mm_shuffle_epi32(row3, _MM_SHUFFLE(1, 0, 3, 2) ));
  row3 = _mm_add_epi16(row3, _mm_shuffle_epi32(row3, _MM_SHUFFLE(0, 1, 0, 1) ));
  row3 = _mm_add_epi16(row3, _mm_shufflelo_epi16(row3, _MM_SHUFFLE(0, 1, 0, 1) ));

  unsigned sum = _mm_extract_epi16(row3, 0);
  unsigned satd = (sum + 1) >> 1;

  return satd;
}


static void satd_8bit_4x4_dual_avx2(
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out) 
{

  __m256i original = _mm256_broadcastsi128_si256(_mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)orig)));
  __m256i pred = _mm256_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)preds[0]));
  pred = _mm256_inserti128_si256(pred, _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)preds[1])), 1);

  __m256i diff_lo = _mm256_sub_epi16(pred, original);

  original = _mm256_broadcastsi128_si256(_mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(orig + 8))));
  pred = _mm256_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(preds[0] + 8)));
  pred = _mm256_inserti128_si256(pred, _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(preds[1] + 8))), 1);

  __m256i diff_hi = _mm256_sub_epi16(pred, original);

  //Hor
  __m256i row0 = _mm256_hadd_epi16(diff_lo, diff_hi);
  __m256i row1 = _mm256_hsub_epi16(diff_lo, diff_hi);

  __m256i row2 = _mm256_hadd_epi16(row0, row1);
  __m256i row3 = _mm256_hsub_epi16(row0, row1);

  //Ver
  row0 = _mm256_hadd_epi16(row2, row3);
  row1 = _mm256_hsub_epi16(row2, row3);

  row2 = _mm256_hadd_epi16(row0, row1);
  row3 = _mm256_hsub_epi16(row0, row1);

  //Abs and sum
  row2 = _mm256_abs_epi16(row2);
  row3 = _mm256_abs_epi16(row3);

  row3 = _mm256_add_epi16(row2, row3);

  row3 = _mm256_add_epi16(row3, _mm256_shuffle_epi32(row3, _MM_SHUFFLE(1, 0, 3, 2) ));
  row3 = _mm256_add_epi16(row3, _mm256_shuffle_epi32(row3, _MM_SHUFFLE(0, 1, 0, 1) ));
  row3 = _mm256_add_epi16(row3, _mm256_shufflelo_epi16(row3, _MM_SHUFFLE(0, 1, 0, 1) ));

  unsigned sum1 = _mm_extract_epi16(_mm256_castsi256_si128(row3), 0);
  sum1 = (sum1 + 1) >> 1;

  unsigned sum2 = _mm_extract_epi16(_mm256_extracti128_si256(row3, 1), 0);
  sum2 = (sum2 + 1) >> 1;

  satds_out[0] = sum1;
  satds_out[1] = sum2;
}

static INLINE void hor_transform_row_avx2(__m128i* row){
  
  __m128i mask_pos = _mm_set1_epi16(1);
  __m128i mask_neg = _mm_set1_epi16(-1);
  __m128i sign_mask = _mm_unpacklo_epi64(mask_pos, mask_neg);
  __m128i temp = _mm_shuffle_epi32(*row, _MM_SHUFFLE(1, 0, 3, 2));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);

  sign_mask = _mm_unpacklo_epi32(mask_pos, mask_neg);
  temp = _mm_shuffle_epi32(*row, _MM_SHUFFLE(2, 3, 0, 1));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);

  sign_mask = _mm_unpacklo_epi16(mask_pos, mask_neg);
  temp = _mm_shufflelo_epi16(*row, _MM_SHUFFLE(2,3,0,1));
  temp = _mm_shufflehi_epi16(temp, _MM_SHUFFLE(2,3,0,1));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);
}

static INLINE void hor_transform_row_dual_avx2(__m256i* row){
  
  __m256i mask_pos = _mm256_set1_epi16(1);
  __m256i mask_neg = _mm256_set1_epi16(-1);
  __m256i sign_mask = _mm256_unpacklo_epi64(mask_pos, mask_neg);
  __m256i temp = _mm256_shuffle_epi32(*row, _MM_SHUFFLE(1, 0, 3, 2));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);

  sign_mask = _mm256_unpacklo_epi32(mask_pos, mask_neg);
  temp = _mm256_shuffle_epi32(*row, _MM_SHUFFLE(2, 3, 0, 1));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);

  sign_mask = _mm256_unpacklo_epi16(mask_pos, mask_neg);
  temp = _mm256_shufflelo_epi16(*row, _MM_SHUFFLE(2,3,0,1));
  temp = _mm256_shufflehi_epi16(temp, _MM_SHUFFLE(2,3,0,1));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);
}

static INLINE void add_sub_avx2(__m128i *out, __m128i *in, unsigned out_idx0, unsigned out_idx1, unsigned in_idx0, unsigned in_idx1)
{
  out[out_idx0] = _mm_add_epi16(in[in_idx0], in[in_idx1]);
  out[out_idx1] = _mm_sub_epi16(in[in_idx0], in[in_idx1]);
}

static INLINE void ver_transform_block_avx2(__m128i (*rows)[8]){

  __m128i temp0[8];
  add_sub_avx2(temp0, (*rows), 0, 1, 0, 1);
  add_sub_avx2(temp0, (*rows), 2, 3, 2, 3);
  add_sub_avx2(temp0, (*rows), 4, 5, 4, 5);
  add_sub_avx2(temp0, (*rows), 6, 7, 6, 7);

  __m128i temp1[8];
  add_sub_avx2(temp1, temp0, 0, 1, 0, 2);
  add_sub_avx2(temp1, temp0, 2, 3, 1, 3);
  add_sub_avx2(temp1, temp0, 4, 5, 4, 6);
  add_sub_avx2(temp1, temp0, 6, 7, 5, 7);

  add_sub_avx2((*rows), temp1, 0, 1, 0, 4);
  add_sub_avx2((*rows), temp1, 2, 3, 1, 5);
  add_sub_avx2((*rows), temp1, 4, 5, 2, 6);
  add_sub_avx2((*rows), temp1, 6, 7, 3, 7);
  
}

static INLINE void add_sub_dual_avx2(__m256i *out, __m256i *in, unsigned out_idx0, unsigned out_idx1, unsigned in_idx0, unsigned in_idx1)
{
  out[out_idx0] = _mm256_add_epi16(in[in_idx0], in[in_idx1]);
  out[out_idx1] = _mm256_sub_epi16(in[in_idx0], in[in_idx1]);
}


static INLINE void ver_transform_block_dual_avx2(__m256i (*rows)[8]){

  __m256i temp0[8];
  add_sub_dual_avx2(temp0, (*rows), 0, 1, 0, 1);
  add_sub_dual_avx2(temp0, (*rows), 2, 3, 2, 3);
  add_sub_dual_avx2(temp0, (*rows), 4, 5, 4, 5);
  add_sub_dual_avx2(temp0, (*rows), 6, 7, 6, 7);

  __m256i temp1[8];
  add_sub_dual_avx2(temp1, temp0, 0, 1, 0, 2);
  add_sub_dual_avx2(temp1, temp0, 2, 3, 1, 3);
  add_sub_dual_avx2(temp1, temp0, 4, 5, 4, 6);
  add_sub_dual_avx2(temp1, temp0, 6, 7, 5, 7);

  add_sub_dual_avx2((*rows), temp1, 0, 1, 0, 4);
  add_sub_dual_avx2((*rows), temp1, 2, 3, 1, 5);
  add_sub_dual_avx2((*rows), temp1, 4, 5, 2, 6);
  add_sub_dual_avx2((*rows), temp1, 6, 7, 3, 7);
  
}

INLINE static void haddwd_accumulate_avx2(__m128i *accumulate, __m128i *ver_row)
{
  __m128i abs_value = _mm_abs_epi16(*ver_row);
  *accumulate = _mm_add_epi32(*accumulate, _mm_madd_epi16(abs_value, _mm_set1_epi16(1)));
}

INLINE static void haddwd_accumulate_dual_avx2(__m256i *accumulate, __m256i *ver_row)
{
  __m256i abs_value = _mm256_abs_epi16(*ver_row);
  *accumulate = _mm256_add_epi32(*accumulate, _mm256_madd_epi16(abs_value, _mm256_set1_epi16(1)));
}

INLINE static unsigned sum_block_avx2(__m128i *ver_row)
{
  __m128i sad = _mm_setzero_si128();
  haddwd_accumulate_avx2(&sad, ver_row + 0);
  haddwd_accumulate_avx2(&sad, ver_row + 1);
  haddwd_accumulate_avx2(&sad, ver_row + 2);
  haddwd_accumulate_avx2(&sad, ver_row + 3); 
  haddwd_accumulate_avx2(&sad, ver_row + 4);
  haddwd_accumulate_avx2(&sad, ver_row + 5);
  haddwd_accumulate_avx2(&sad, ver_row + 6);
  haddwd_accumulate_avx2(&sad, ver_row + 7);

  sad = _mm_add_epi32(sad, _mm_shuffle_epi32(sad, _MM_SHUFFLE(1, 0, 3, 2)));
  sad = _mm_add_epi32(sad, _mm_shuffle_epi32(sad, _MM_SHUFFLE(0, 1, 0, 1)));

  return _mm_cvtsi128_si32(sad);
}

INLINE static void sum_block_dual_avx2(__m256i *ver_row, unsigned *sum0, unsigned *sum1)
{
  __m256i sad = _mm256_setzero_si256();
  haddwd_accumulate_dual_avx2(&sad, ver_row + 0);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 1);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 2);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 3); 
  haddwd_accumulate_dual_avx2(&sad, ver_row + 4);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 5);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 6);
  haddwd_accumulate_dual_avx2(&sad, ver_row + 7);

  sad = _mm256_add_epi32(sad, _mm256_shuffle_epi32(sad, _MM_SHUFFLE(1, 0, 3, 2)));
  sad = _mm256_add_epi32(sad, _mm256_shuffle_epi32(sad, _MM_SHUFFLE(0, 1, 0, 1)));

  *sum0 = _mm_cvtsi128_si32(_mm256_extracti128_si256(sad, 0));
  *sum1 = _mm_cvtsi128_si32(_mm256_extracti128_si256(sad, 1));
}

INLINE static __m128i diff_row_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  __m128i buf1_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)buf1));
  __m128i buf2_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)buf2));
  return _mm_sub_epi16(buf1_row, buf2_row);
}

INLINE static __m256i diff_row_dual_avx2(const uint8_t *buf1, const uint8_t *buf2, const uint8_t *orig)
{
  __m128i temp1 = _mm_loadl_epi64((__m128i*)buf1);
  __m128i temp2 = _mm_loadl_epi64((__m128i*)buf2);
  __m128i temp3 = _mm_loadl_epi64((__m128i*)orig);
  __m256i buf1_row = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(temp1, temp2));
  __m256i buf2_row = _mm256_cvtepu8_epi16(_mm_broadcastq_epi64(temp3));

  return _mm256_sub_epi16(buf1_row, buf2_row);
}

INLINE static void diff_blocks_avx2(__m128i (*row_diff)[8],
                                                           const uint8_t * buf1, unsigned stride1,
                                                           const uint8_t * orig, unsigned stride_orig)
{
  (*row_diff)[0] = diff_row_avx2(buf1 + 0 * stride1, orig + 0 * stride_orig);
  (*row_diff)[1] = diff_row_avx2(buf1 + 1 * stride1, orig + 1 * stride_orig);
  (*row_diff)[2] = diff_row_avx2(buf1 + 2 * stride1, orig + 2 * stride_orig);
  (*row_diff)[3] = diff_row_avx2(buf1 + 3 * stride1, orig + 3 * stride_orig);
  (*row_diff)[4] = diff_row_avx2(buf1 + 4 * stride1, orig + 4 * stride_orig);
  (*row_diff)[5] = diff_row_avx2(buf1 + 5 * stride1, orig + 5 * stride_orig);
  (*row_diff)[6] = diff_row_avx2(buf1 + 6 * stride1, orig + 6 * stride_orig);
  (*row_diff)[7] = diff_row_avx2(buf1 + 7 * stride1, orig + 7 * stride_orig);

}

INLINE static void diff_blocks_dual_avx2(__m256i (*row_diff)[8],
                                                           const uint8_t * buf1, unsigned stride1,
                                                           const uint8_t * buf2, unsigned stride2,
                                                           const uint8_t * orig, unsigned stride_orig)
{
  (*row_diff)[0] = diff_row_dual_avx2(buf1 + 0 * stride1, buf2 + 0 * stride2, orig + 0 * stride_orig);
  (*row_diff)[1] = diff_row_dual_avx2(buf1 + 1 * stride1, buf2 + 1 * stride2, orig + 1 * stride_orig);
  (*row_diff)[2] = diff_row_dual_avx2(buf1 + 2 * stride1, buf2 + 2 * stride2, orig + 2 * stride_orig);
  (*row_diff)[3] = diff_row_dual_avx2(buf1 + 3 * stride1, buf2 + 3 * stride2, orig + 3 * stride_orig);
  (*row_diff)[4] = diff_row_dual_avx2(buf1 + 4 * stride1, buf2 + 4 * stride2, orig + 4 * stride_orig);
  (*row_diff)[5] = diff_row_dual_avx2(buf1 + 5 * stride1, buf2 + 5 * stride2, orig + 5 * stride_orig);
  (*row_diff)[6] = diff_row_dual_avx2(buf1 + 6 * stride1, buf2 + 6 * stride2, orig + 6 * stride_orig);
  (*row_diff)[7] = diff_row_dual_avx2(buf1 + 7 * stride1, buf2 + 7 * stride2, orig + 7 * stride_orig);

}

INLINE static void hor_transform_block_avx2(__m128i (*row_diff)[8])
{
  hor_transform_row_avx2((*row_diff) + 0);
  hor_transform_row_avx2((*row_diff) + 1);
  hor_transform_row_avx2((*row_diff) + 2);
  hor_transform_row_avx2((*row_diff) + 3);
  hor_transform_row_avx2((*row_diff) + 4);
  hor_transform_row_avx2((*row_diff) + 5);
  hor_transform_row_avx2((*row_diff) + 6);
  hor_transform_row_avx2((*row_diff) + 7);
}

INLINE static void hor_transform_block_dual_avx2(__m256i (*row_diff)[8])
{
  hor_transform_row_dual_avx2((*row_diff) + 0);
  hor_transform_row_dual_avx2((*row_diff) + 1);
  hor_transform_row_dual_avx2((*row_diff) + 2);
  hor_transform_row_dual_avx2((*row_diff) + 3);
  hor_transform_row_dual_avx2((*row_diff) + 4);
  hor_transform_row_dual_avx2((*row_diff) + 5);
  hor_transform_row_dual_avx2((*row_diff) + 6);
  hor_transform_row_dual_avx2((*row_diff) + 7);
}

static void kvz_satd_8bit_8x8_general_dual_avx2(const uint8_t * buf1, unsigned stride1,
                                                const uint8_t * buf2, unsigned stride2,
                                                const uint8_t * orig, unsigned stride_orig,
                                                unsigned *sum0, unsigned *sum1)
{
  __m256i temp[8];

  diff_blocks_dual_avx2(&temp, buf1, stride1, buf2, stride2, orig, stride_orig);
  hor_transform_block_dual_avx2(&temp);
  ver_transform_block_dual_avx2(&temp);
  
  sum_block_dual_avx2(temp, sum0, sum1);

  *sum0 = (*sum0 + 2) >> 2;
  *sum1 = (*sum1 + 2) >> 2;
}

/**
* \brief  Calculate SATD between two 4x4 blocks inside bigger arrays.
*/
static unsigned kvz_satd_4x4_subblock_8bit_avx2(const uint8_t * buf1,
                                                const int32_t     stride1,
                                                const uint8_t * buf2,
                                                const int32_t     stride2)
{
  // TODO: AVX2 implementation
  return kvz_satd_4x4_subblock_generic(buf1, stride1, buf2, stride2);
}

static void kvz_satd_4x4_subblock_quad_avx2(const uint8_t *preds[4],
                                       const int stride,
                                       const uint8_t *orig,
                                       const int orig_stride,
                                       unsigned costs[4])
{
  // TODO: AVX2 implementation
  kvz_satd_4x4_subblock_quad_generic(preds, stride, orig, orig_stride, costs);
}

static unsigned satd_8x8_subblock_8bit_avx2(const uint8_t * buf1, unsigned stride1, const uint8_t * buf2, unsigned stride2)
{
  __m128i temp[8];

  diff_blocks_avx2(&temp, buf1, stride1, buf2, stride2);
  hor_transform_block_avx2(&temp);
  ver_transform_block_avx2(&temp);
  
  unsigned sad = sum_block_avx2(temp);

  unsigned result = (sad + 2) >> 2;
  return result;
}

static void satd_8x8_subblock_quad_avx2(const uint8_t **preds,
  const int stride,
  const uint8_t *orig,
  const int orig_stride,
  unsigned *costs)
{
  kvz_satd_8bit_8x8_general_dual_avx2(preds[0], stride, preds[1], stride, orig, orig_stride, &costs[0], &costs[1]);
  kvz_satd_8bit_8x8_general_dual_avx2(preds[2], stride, preds[3], stride, orig, orig_stride, &costs[2], &costs[3]);
}

SATD_NxN(8bit_avx2,  8)
SATD_NxN(8bit_avx2, 16)
SATD_NxN(8bit_avx2, 32)
SATD_NxN(8bit_avx2, 64)
SATD_ANY_SIZE(8bit_avx2)

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NXN_DUAL_AVX2(n) \
static void satd_8bit_ ## n ## x ## n ## _dual_avx2( \
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out)  \
{ \
  unsigned x, y; \
  satds_out[0] = 0; \
  satds_out[1] = 0; \
  unsigned sum1 = 0; \
  unsigned sum2 = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  kvz_satd_8bit_8x8_general_dual_avx2(&preds[0][row + x], (n), &preds[1][row + x], (n), &orig[row + x], (n), &sum1, &sum2); \
  satds_out[0] += sum1; \
  satds_out[1] += sum2; \
    } \
    } \
  satds_out[0] >>= (KVZ_BIT_DEPTH-8); \
  satds_out[1] >>= (KVZ_BIT_DEPTH-8); \
}

static void satd_8bit_8x8_dual_avx2(
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out) 
{ 
  unsigned x, y; 
  satds_out[0] = 0;
  satds_out[1] = 0;
  unsigned sum1 = 0;
  unsigned sum2 = 0;
  for (y = 0; y < (8); y += 8) { 
  unsigned row = y * (8); 
  for (x = 0; x < (8); x += 8) { 
  kvz_satd_8bit_8x8_general_dual_avx2(&preds[0][row + x], (8), &preds[1][row + x], (8), &orig[row + x], (8), &sum1, &sum2); 
  satds_out[0] += sum1;
  satds_out[1] += sum2;
      } 
      } 
  satds_out[0] >>= (KVZ_BIT_DEPTH-8);
  satds_out[1] >>= (KVZ_BIT_DEPTH-8);
}

//SATD_NXN_DUAL_AVX2(8) //Use the non-macro version
SATD_NXN_DUAL_AVX2(16)
SATD_NXN_DUAL_AVX2(32)
SATD_NXN_DUAL_AVX2(64)

#define SATD_ANY_SIZE_MULTI_AVX2(suffix, num_parallel_blocks) \
  static cost_pixel_any_size_multi_func satd_any_size_## suffix; \
  static void satd_any_size_ ## suffix ( \
      int width, int height, \
      const uint8_t **preds, \
      const int stride, \
      const uint8_t *orig, \
      const int orig_stride, \
      unsigned num_modes, \
      unsigned *costs_out, \
      int8_t *valid) \
  { \
    unsigned sums[num_parallel_blocks] = { 0 }; \
    const uint8_t *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };\
    const uint8_t *orig_ptr = orig; \
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

SATD_ANY_SIZE_MULTI_AVX2(quad_avx2, 4)


static unsigned pixels_calc_ssd_avx2(const uint8_t *const ref, const uint8_t *const rec,
                 const int ref_stride, const int rec_stride,
                 const int width)
{
  __m256i ssd_part;
  __m256i diff = _mm256_setzero_si256();
  __m128i sum;

  __m256i ref_epi16;
  __m256i rec_epi16;

  __m128i ref_row0, ref_row1, ref_row2, ref_row3;
  __m128i rec_row0, rec_row1, rec_row2, rec_row3;

  int ssd;

  switch (width) {

  case 4:

    ref_row0 = _mm_cvtsi32_si128(*(int32_t*)&(ref[0 * ref_stride]));
    ref_row1 = _mm_cvtsi32_si128(*(int32_t*)&(ref[1 * ref_stride]));
    ref_row2 = _mm_cvtsi32_si128(*(int32_t*)&(ref[2 * ref_stride]));
    ref_row3 = _mm_cvtsi32_si128(*(int32_t*)&(ref[3 * ref_stride]));

    ref_row0 = _mm_unpacklo_epi32(ref_row0, ref_row1);
    ref_row1 = _mm_unpacklo_epi32(ref_row2, ref_row3);
    ref_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(ref_row0, ref_row1) );

    rec_row0 = _mm_cvtsi32_si128(*(int32_t*)&(rec[0 * rec_stride]));
    rec_row1 = _mm_cvtsi32_si128(*(int32_t*)&(rec[1 * rec_stride]));
    rec_row2 = _mm_cvtsi32_si128(*(int32_t*)&(rec[2 * rec_stride]));
    rec_row3 = _mm_cvtsi32_si128(*(int32_t*)&(rec[3 * rec_stride]));

    rec_row0 = _mm_unpacklo_epi32(rec_row0, rec_row1);
    rec_row1 = _mm_unpacklo_epi32(rec_row2, rec_row3);
    rec_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(rec_row0, rec_row1) );

    diff = _mm256_sub_epi16(ref_epi16, rec_epi16);
    ssd_part =  _mm256_madd_epi16(diff, diff);

    sum = _mm_add_epi32(_mm256_castsi256_si128(ssd_part), _mm256_extracti128_si256(ssd_part, 1));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));

    ssd = _mm_cvtsi128_si32(sum);

    return ssd >> (2*(KVZ_BIT_DEPTH-8));
    break;

  default:

    ssd_part = _mm256_setzero_si256();
    for (int y = 0; y < width; y += 8) {
      for (int x = 0; x < width; x += 8) {
        for (int i = 0; i < 8; i += 2) {
          ref_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)&(ref[x + (y + i) * ref_stride])), _mm_loadl_epi64((__m128i*)&(ref[x + (y + i + 1) * ref_stride]))));
          rec_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)&(rec[x + (y + i) * rec_stride])), _mm_loadl_epi64((__m128i*)&(rec[x + (y + i + 1) * rec_stride]))));
          diff = _mm256_sub_epi16(ref_epi16, rec_epi16);
          ssd_part = _mm256_add_epi32(ssd_part, _mm256_madd_epi16(diff, diff));
        }
      }
    }

    sum = _mm_add_epi32(_mm256_castsi256_si128(ssd_part), _mm256_extracti128_si256(ssd_part, 1));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));

    ssd = _mm_cvtsi128_si32(sum);

    return ssd >> (2*(KVZ_BIT_DEPTH-8));
    break;
  }
}

static INLINE void scatter_ymm_4x8_8bit(kvz_pixel * dst, __m256i ymm, unsigned dst_stride)
{
  __m128i ymm_lo = _mm256_castsi256_si128(ymm);
  __m128i ymm_hi = _mm256_extracti128_si256(ymm, 1);
  *(uint32_t *)dst = _mm_cvtsi128_si32(ymm_lo); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 1); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 2); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 3); dst += dst_stride;
  *(uint32_t *)dst = _mm_cvtsi128_si32(ymm_hi); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 1); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 2); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 3);
}

static INLINE void scatter_ymm_8x4_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m256d ymm_as_m256d = _mm256_castsi256_pd(ymm);
  __m128d ymm_lo = _mm256_castpd256_pd128(ymm_as_m256d);
  __m128d ymm_hi = _mm256_extractf128_pd(ymm_as_m256d, 1);
  _mm_storel_pd((double*)dst, ymm_lo); dst += dst_stride;
  _mm_storeh_pd((double*)dst, ymm_lo); dst += dst_stride;
  _mm_storel_pd((double*)dst, ymm_hi); dst += dst_stride;
  _mm_storeh_pd((double*)dst, ymm_hi);
}

static INLINE void scatter_ymm_16x2_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m128i ymm_lo = _mm256_castsi256_si128(ymm);
  __m128i ymm_hi = _mm256_extracti128_si256(ymm, 1);
  _mm_storeu_si128((__m128i *)dst, ymm_lo); dst += dst_stride;
  _mm_storeu_si128((__m128i *)dst, ymm_hi);
}

static INLINE void scatter_ymm_12x2_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m256i mask_a = _mm256_setr_epi32(-1, -1, -1, 0, 0, 0, 0, 0);
  __m256i mask_b = _mm256_setr_epi32(0, 0, 0, -1, -1, -1, 0, 0);
  _mm256_maskstore_epi32((int32_t*)dst, mask_a, ymm); dst += dst_stride - 3 * 4;
  _mm256_maskstore_epi32((int32_t*)dst, mask_b, ymm);
}

static INLINE void bipred_average_px_px_template_avx2(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  bool has_pow2_width = _mm_popcnt_u32(pu_w) == 1;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {

      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_L0 = _mm256_loadu_si256((__m256i*)&px_L0[i]);
      __m256i sample_L1 = _mm256_loadu_si256((__m256i*)&px_L1[i]);
      __m256i avg       = _mm256_avg_epu8(sample_L0, sample_L1);

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i *)&dst[y * dst_stride + x], avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 are not used to simplify the loop
      __m256i mask      = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m256i sample_L0 = _mm256_maskload_epi64((const long long*)&px_L0[i], mask);
      __m256i sample_L1 = _mm256_maskload_epi64((const long long*)&px_L1[i], mask);
      __m256i avg       = _mm256_avg_epu8(sample_L0, sample_L1);

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      __m128i sample_L0, sample_L1, avg;
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          sample_L0 = _mm_loadu_si128((__m128i*)&px_L0[i]);
          sample_L1 = _mm_loadu_si128((__m128i*)&px_L1[i]);
          avg       = _mm_avg_epu8(sample_L0, sample_L1);
          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m128i mask      = _mm_setr_epi32(-1, -1, -1, 0);
          __m128i sample_L0 = _mm_maskload_epi32((const int*)(&px_L0[i]), mask);
          __m128i sample_L1 = _mm_maskload_epi32((const int*)(&px_L1[i]), mask);
          __m128i avg       = _mm_avg_epu8(sample_L0, sample_L1);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  }
}

static INLINE void bipred_average_px_px_avx2(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_px_px_template_avx2(dst, px_L0, px_L1,  4, pu_h, dst_stride); break;
      case  8: bipred_average_px_px_template_avx2(dst, px_L0, px_L1,  8, pu_h, dst_stride); break;
      case 16: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 16, pu_h, dst_stride); break;
      case 32: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 32, pu_h, dst_stride); break;
      case 64: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 64, pu_h, dst_stride); break;

      case  6: bipred_average_px_px_template_avx2(dst, px_L0, px_L1,  6, pu_h, dst_stride); break;
      case 12: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 12, pu_h, dst_stride); break;
      case 24: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 24, pu_h, dst_stride); break;
      case 48: bipred_average_px_px_template_avx2(dst, px_L0, px_L1, 48, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
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
}

static INLINE void bipred_average_im_im_template_avx2(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t scalar_offset = 1 << (shift - 1);
  __m256i offset = _mm256_set1_epi32(scalar_offset);

  bool has_pow2_width = _mm_popcnt_u32(pu_w) == 1;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {
      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_L0_a_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
      __m256i sample_L1_a_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);
      __m256i sample_L0_b_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i + 16]);
      __m256i sample_L1_b_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i + 16]);

      __m256i sample_L0_L1_a_lo = _mm256_unpacklo_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_a_hi = _mm256_unpackhi_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_b_lo = _mm256_unpacklo_epi16(sample_L0_b_16bit, sample_L1_b_16bit);
      __m256i sample_L0_L1_b_hi = _mm256_unpackhi_epi16(sample_L0_b_16bit, sample_L1_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_L0_L1_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_L0_L1_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i*)&dst[y * dst_stride + x], avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 are not used to simplify the loop
      __m256i mask              = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m256i sample_L0_a_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
      __m256i sample_L1_a_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);
      __m256i sample_L0_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im_L0[i + 16]));
      __m256i sample_L1_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im_L1[i + 16]));

      __m256i sample_L0_L1_a_lo = _mm256_unpacklo_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_a_hi = _mm256_unpackhi_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_b_lo = _mm256_unpacklo_epi16(sample_L0_b_16bit, sample_L1_b_16bit);
      __m256i sample_L0_L1_b_hi = _mm256_unpackhi_epi16(sample_L0_b_16bit, sample_L1_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_L0_L1_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_L0_L1_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          __m256i sample_L0_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
          __m256i sample_L1_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);

          __m256i sample_L0_L1_lo = _mm256_unpacklo_epi16(sample_L0_16bit, sample_L1_16bit);
          __m256i sample_L0_L1_hi = _mm256_unpackhi_epi16(sample_L0_16bit, sample_L1_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_L0_L1_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_L0_L1_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m256i mask            = _mm256_setr_epi64x(-1, -1, -1, 0);
          __m256i sample_L0_16bit = _mm256_maskload_epi64((const long long*)(&im_L0[i]), mask);
          __m256i sample_L1_16bit = _mm256_maskload_epi64((const long long*)(&im_L1[i]), mask);

          __m256i sample_L0_L1_lo = _mm256_unpacklo_epi16(sample_L0_16bit, sample_L1_16bit);
          __m256i sample_L0_L1_hi = _mm256_unpackhi_epi16(sample_L0_16bit, sample_L1_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_lo, all_ones);
          __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_hi, all_ones);

          avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
          avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);

          avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
          avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  }
}

static void bipred_average_im_im_avx2(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_im_im_template_avx2(dst, im_L0, im_L1,  4, pu_h, dst_stride); break;
      case  8: bipred_average_im_im_template_avx2(dst, im_L0, im_L1,  8, pu_h, dst_stride); break;
      case 16: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 16, pu_h, dst_stride); break;
      case 32: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 32, pu_h, dst_stride); break;
      case 64: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 64, pu_h, dst_stride); break;

      case  6: bipred_average_im_im_template_avx2(dst, im_L0, im_L1,  6, pu_h, dst_stride); break;
      case 12: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 12, pu_h, dst_stride); break;
      case 24: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 24, pu_h, dst_stride); break;
      case 48: bipred_average_im_im_template_avx2(dst, im_L0, im_L1, 48, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
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
}

static INLINE void bipred_average_px_im_template_avx2(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t scalar_offset = 1 << (shift - 1);
  __m256i offset = _mm256_set1_epi32(scalar_offset);

  bool has_pow2_width = _mm_popcnt_u32(pu_w) == 1;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {

      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_px_a_16bit = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i*)&px[i]));
      __m256i sample_px_b_16bit = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i*)&px[i + 16]));
      sample_px_a_16bit         = _mm256_slli_epi16(sample_px_a_16bit, 14 - KVZ_BIT_DEPTH);
      sample_px_b_16bit         = _mm256_slli_epi16(sample_px_b_16bit, 14 - KVZ_BIT_DEPTH);
      __m256i sample_im_a_16bit = _mm256_loadu_si256((__m256i*)&im[i]);
      __m256i sample_im_b_16bit = _mm256_loadu_si256((__m256i*)&im[i + 16]);

      __m256i sample_px_im_a_lo = _mm256_unpacklo_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_a_hi = _mm256_unpackhi_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_b_lo = _mm256_unpacklo_epi16(sample_px_b_16bit, sample_im_b_16bit);
      __m256i sample_px_im_b_hi = _mm256_unpackhi_epi16(sample_px_b_16bit, sample_im_b_16bit);

      __m256i all_ones  = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_px_im_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_px_im_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_px_im_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_px_im_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i*)&dst[y * dst_stride + x], avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 / 32 bits of the 128 are not used to simplify the loop
      __m256i mask              = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m128i sample_px_a_8bit  = _mm_loadu_si128((__m128i*)&px[i]);
      __m128i sample_px_b_8bit  = _mm_loadl_epi64((__m128i*)&px[i + 16]);
      __m256i sample_px_a_16bit = _mm256_cvtepu8_epi16(sample_px_a_8bit);
      __m256i sample_px_b_16bit = _mm256_cvtepu8_epi16(sample_px_b_8bit);
      sample_px_a_16bit         = _mm256_slli_epi16(sample_px_a_16bit, 14 - KVZ_BIT_DEPTH);
      sample_px_b_16bit         = _mm256_slli_epi16(sample_px_b_16bit, 14 - KVZ_BIT_DEPTH);
      __m256i sample_im_a_16bit = _mm256_loadu_si256((__m256i*)&im[i]);
      __m256i sample_im_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im[i + 16]));

      __m256i sample_px_im_a_lo = _mm256_unpacklo_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_a_hi = _mm256_unpackhi_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_b_lo = _mm256_unpacklo_epi16(sample_px_b_16bit, sample_im_b_16bit);
      __m256i sample_px_im_b_hi = _mm256_unpackhi_epi16(sample_px_b_16bit, sample_im_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_px_im_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_px_im_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_px_im_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_px_im_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          assert(0 && "Unexpected block width.");
          break;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          __m128i sample_px_8bit  = _mm_loadu_si128((__m128i*)&px[i]);
          __m256i sample_px_16bit = _mm256_cvtepu8_epi16(sample_px_8bit);
          sample_px_16bit         = _mm256_slli_epi16(sample_px_16bit, 14 - KVZ_BIT_DEPTH);
          __m256i sample_im_16bit = _mm256_loadu_si256((__m256i*)&im[i]);

          __m256i sample_px_im_lo = _mm256_unpacklo_epi16(sample_px_16bit, sample_im_16bit);
          __m256i sample_px_im_hi = _mm256_unpackhi_epi16(sample_px_16bit, sample_im_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_px_im_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_px_im_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m128i mask128         = _mm_setr_epi32(-1, -1, -1, 0);
          __m128i sample_px_8bit  = _mm_maskload_epi32((const int*)(&px[i]), mask128);

          __m256i mask            = _mm256_setr_epi64x(-1, -1, -1, 0);
          __m256i sample_px_16bit = _mm256_cvtepu8_epi16(sample_px_8bit);
          sample_px_16bit         = _mm256_slli_epi16(sample_px_16bit, 14 - KVZ_BIT_DEPTH);
          __m256i sample_im_16bit = _mm256_maskload_epi64((const long long*)(&im[i]), mask);

          __m256i sample_px_im_lo = _mm256_unpacklo_epi16(sample_px_16bit, sample_im_16bit);
          __m256i sample_px_im_hi = _mm256_unpackhi_epi16(sample_px_16bit, sample_im_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_px_im_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_px_im_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  }
}

static void bipred_average_px_im_avx2(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_px_im_template_avx2(dst, px, im,  4, pu_h, dst_stride); break;
      case  8: bipred_average_px_im_template_avx2(dst, px, im,  8, pu_h, dst_stride); break;
      case 16: bipred_average_px_im_template_avx2(dst, px, im, 16, pu_h, dst_stride); break;
      case 32: bipred_average_px_im_template_avx2(dst, px, im, 32, pu_h, dst_stride); break;
      case 64: bipred_average_px_im_template_avx2(dst, px, im, 64, pu_h, dst_stride); break;

      case  6: bipred_average_px_im_template_avx2(dst, px, im,  6, pu_h, dst_stride); break;
      case 12: bipred_average_px_im_template_avx2(dst, px, im, 12, pu_h, dst_stride); break;
      case 24: bipred_average_px_im_template_avx2(dst, px, im, 24, pu_h, dst_stride); break;
      case 48: bipred_average_px_im_template_avx2(dst, px, im, 48, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
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
}

static void bipred_average_avx2(lcu_t *const lcu,
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
      bipred_average_px_px_avx2(lcu->rec.y + pb_offset, px_L0->y, px_L1->y, pu_w, pu_h, LCU_WIDTH);

    } else if ((im_flags_L0 & 1) && (im_flags_L1 & 1)) {
      bipred_average_im_im_avx2(lcu->rec.y + pb_offset, im_L0->y, im_L1->y, pu_w, pu_h, LCU_WIDTH);

    } else {
      kvz_pixel *src_px    = (im_flags_L0 & 1) ? px_L1->y : px_L0->y;
      kvz_pixel_im *src_im = (im_flags_L0 & 1) ? im_L0->y : im_L1->y;
      bipred_average_px_im_avx2(lcu->rec.y + pb_offset, src_px, src_im, pu_w, pu_h, LCU_WIDTH);
    }
  }
  if (predict_chroma) {
    unsigned pb_offset = SUB_SCU(pu_y) / 2 * LCU_WIDTH_C + SUB_SCU(pu_x) / 2;
    unsigned pb_w = pu_w / 2;
    unsigned pb_h = pu_h / 2;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_px_px_avx2(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_px_px_avx2(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, LCU_WIDTH_C);

    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_im_im_avx2(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_im_im_avx2(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, LCU_WIDTH_C);

    } else {
      kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_px_im_avx2(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, LCU_WIDTH_C);
      bipred_average_px_im_avx2(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, LCU_WIDTH_C);
    }
  }
}

static optimized_sad_func_ptr_t get_optimized_sad_avx2(int32_t width)
{
  if (width == 0)
    return reg_sad_w0;
  if (width == 4)
    return reg_sad_w4;
  if (width == 8)
    return reg_sad_w8;
  if (width == 12)
    return reg_sad_w12;
  if (width == 16)
    return reg_sad_w16;
  if (width == 24)
    return reg_sad_w24;
  if (width == 32)
    return reg_sad_w32;
  if (width == 64)
    return reg_sad_w64;
  else
    return NULL;
}

static uint32_t ver_sad_avx2(const uint8_t *pic_data, const uint8_t *ref_data,
                             int32_t width, int32_t height, uint32_t stride)
{
  if (width == 0)
    return 0;
  if (width == 4)
    return ver_sad_w4(pic_data, ref_data, height, stride);
  if (width == 8)
    return ver_sad_w8(pic_data, ref_data, height, stride);
  if (width == 12)
    return ver_sad_w12(pic_data, ref_data, height, stride);
  if (width == 16)
    return ver_sad_w16(pic_data, ref_data, height, stride);
  else
    return ver_sad_arbitrary(pic_data, ref_data, width, height, stride);
}

static uint32_t hor_sad_avx2(const uint8_t *pic_data, const uint8_t *ref_data,
                             int32_t width, int32_t height, uint32_t pic_stride,
                             uint32_t ref_stride, uint32_t left, uint32_t right)
{
  if (width == 4)
    return hor_sad_sse41_w4(pic_data, ref_data, height,
                            pic_stride, ref_stride, left, right);
  if (width == 8)
    return hor_sad_sse41_w8(pic_data, ref_data, height,
                            pic_stride, ref_stride, left, right);
  if (width == 16)
    return hor_sad_sse41_w16(pic_data, ref_data, height,
                             pic_stride, ref_stride, left, right);
  if (width == 32)
    return hor_sad_avx2_w32 (pic_data, ref_data, height,
                             pic_stride, ref_stride, left, right);
  else
    return hor_sad_sse41_arbitrary(pic_data, ref_data, width, height,
                                   pic_stride, ref_stride, left, right);
}

static double pixel_var_avx2_largebuf(const uint8_t *buf, const uint32_t len)
{
  const float len_f  = (float)len;
  const __m256i zero = _mm256_setzero_si256();

  int64_t sum;
  size_t i;
  __m256i sums = zero;
  for (i = 0; i + 31 < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));
    __m256i curr_sum = _mm256_sad_epu8(curr, zero);
            sums = _mm256_add_epi64(sums, curr_sum);
  }
  __m128i sum_lo = _mm256_castsi256_si128  (sums);
  __m128i sum_hi = _mm256_extracti128_si256(sums,   1);
  __m128i sum_3  = _mm_add_epi64           (sum_lo, sum_hi);
  __m128i sum_4  = _mm_shuffle_epi32       (sum_3,  _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sum_5  = _mm_add_epi64           (sum_3,  sum_4);

  _mm_storel_epi64((__m128i *)&sum, sum_5);

  // Remaining len mod 32 pixels
  for (; i < len; ++i) {
    sum += buf[i];
  }

  float   mean_f = (float)sum / len_f;
  __m256  mean   = _mm256_set1_ps(mean_f);
  __m256  accum  = _mm256_setzero_ps();

  for (i = 0; i + 31 < len; i += 32) {
    __m128i curr0    = _mm_loadl_epi64((const __m128i *)(buf + i +  0));
    __m128i curr1    = _mm_loadl_epi64((const __m128i *)(buf + i +  8));
    __m128i curr2    = _mm_loadl_epi64((const __m128i *)(buf + i + 16));
    __m128i curr3    = _mm_loadl_epi64((const __m128i *)(buf + i + 24));

    __m256i curr0_32 = _mm256_cvtepu8_epi32(curr0);
    __m256i curr1_32 = _mm256_cvtepu8_epi32(curr1);
    __m256i curr2_32 = _mm256_cvtepu8_epi32(curr2);
    __m256i curr3_32 = _mm256_cvtepu8_epi32(curr3);

    __m256  curr0_f  = _mm256_cvtepi32_ps  (curr0_32);
    __m256  curr1_f  = _mm256_cvtepi32_ps  (curr1_32);
    __m256  curr2_f  = _mm256_cvtepi32_ps  (curr2_32);
    __m256  curr3_f  = _mm256_cvtepi32_ps  (curr3_32);

    __m256  curr0_sd = _mm256_sub_ps       (curr0_f,  mean);
    __m256  curr1_sd = _mm256_sub_ps       (curr1_f,  mean);
    __m256  curr2_sd = _mm256_sub_ps       (curr2_f,  mean);
    __m256  curr3_sd = _mm256_sub_ps       (curr3_f,  mean);

    __m256  curr0_v  = _mm256_mul_ps       (curr0_sd, curr0_sd);
    __m256  curr1_v  = _mm256_mul_ps       (curr1_sd, curr1_sd);
    __m256  curr2_v  = _mm256_mul_ps       (curr2_sd, curr2_sd);
    __m256  curr3_v  = _mm256_mul_ps       (curr3_sd, curr3_sd);

    __m256  curr01   = _mm256_add_ps       (curr0_v,  curr1_v);
    __m256  curr23   = _mm256_add_ps       (curr2_v,  curr3_v);
    __m256  curr     = _mm256_add_ps       (curr01,   curr23);
            accum    = _mm256_add_ps       (accum,    curr);
  }
  __m256d accum_d  = _mm256_castps_pd     (accum);
  __m256d accum2_d = _mm256_permute4x64_pd(accum_d, _MM_SHUFFLE(1, 0, 3, 2));
  __m256  accum2   = _mm256_castpd_ps     (accum2_d);

  __m256  accum3   = _mm256_add_ps        (accum,  accum2);
  __m256  accum4   = _mm256_permute_ps    (accum3, _MM_SHUFFLE(1, 0, 3, 2));
  __m256  accum5   = _mm256_add_ps        (accum3, accum4);
  __m256  accum6   = _mm256_permute_ps    (accum5, _MM_SHUFFLE(2, 3, 0, 1));
  __m256  accum7   = _mm256_add_ps        (accum5, accum6);

  __m128  accum8   = _mm256_castps256_ps128(accum7);
  float   var_sum  = _mm_cvtss_f32         (accum8);

  // Remaining len mod 32 pixels
  for (; i < len; ++i) {
    float diff = buf[i] - mean_f;
    var_sum += diff * diff;
  }

  return  var_sum / len_f;
}

#ifdef INACCURATE_VARIANCE_CALCULATION

// Assumes that u is a power of two
static INLINE uint32_t ilog2(uint32_t u)
{
  return _tzcnt_u32(u);
}

// A B C D | E F G H (8x32b)
//        ==>
// A+B C+D | E+F G+H (4x64b)
static __m256i hsum_epi32_to_epi64(const __m256i v)
{
  const __m256i zero    = _mm256_setzero_si256();
        __m256i v_shufd = _mm256_shuffle_epi32(v, _MM_SHUFFLE(3, 3, 1, 1));
        __m256i sums_32 = _mm256_add_epi32    (v, v_shufd);
        __m256i sums_64 = _mm256_blend_epi32  (sums_32, zero, 0xaa);
  return        sums_64;
}

static double pixel_var_avx2(const uint8_t *buf, const uint32_t len)
{
  assert(sizeof(*buf) == 1);
  assert((len & 31) == 0);

  // Uses Q8.7 numbers to measure mean and deviation, so variances are Q16.14
  const uint64_t sum_maxwid     = ilog2(len) + (8 * sizeof(*buf));
  const __m128i normalize_sum   = _mm_cvtsi32_si128(sum_maxwid - 15); // Normalize mean to [0, 32767], so signed 16-bit subtraction never overflows
  const __m128i debias_sum      = _mm_cvtsi32_si128(1 << (sum_maxwid - 16));
  const float varsum_to_f       = 1.0f / (float)(1 << (14 + ilog2(len)));

  const bool power_of_two = (len & (len - 1)) == 0;
  if (sum_maxwid > 32 || sum_maxwid < 15 || !power_of_two) {
    return pixel_var_avx2_largebuf(buf, len);
  }

  const __m256i zero      = _mm256_setzero_si256();
  const __m256i himask_15 = _mm256_set1_epi16(0x7f00);

  uint64_t vars;
  size_t i;
  __m256i sums = zero;
  for (i = 0; i < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));
    __m256i curr_sum = _mm256_sad_epu8(curr, zero);
            sums = _mm256_add_epi64(sums, curr_sum);
  }
  __m128i sum_lo = _mm256_castsi256_si128  (sums);
  __m128i sum_hi = _mm256_extracti128_si256(sums,   1);
  __m128i sum_3  = _mm_add_epi64           (sum_lo, sum_hi);
  __m128i sum_4  = _mm_shuffle_epi32       (sum_3,  _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sum_5  = _mm_add_epi64           (sum_3,  sum_4);
  __m128i sum_5n = _mm_srl_epi32           (sum_5,  normalize_sum);
          sum_5n = _mm_add_epi32           (sum_5n, debias_sum);

  __m256i sum_n  = _mm256_broadcastw_epi16 (sum_5n);

  __m256i accum = zero;
  for (i = 0; i < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));

    __m256i curr0    = _mm256_slli_epi16  (curr,  7);
    __m256i curr1    = _mm256_srli_epi16  (curr,  1);
            curr0    = _mm256_and_si256   (curr0, himask_15);
            curr1    = _mm256_and_si256   (curr1, himask_15);

    __m256i dev0     = _mm256_sub_epi16   (curr0, sum_n);
    __m256i dev1     = _mm256_sub_epi16   (curr1, sum_n);

    __m256i vars0    = _mm256_madd_epi16  (dev0,  dev0);
    __m256i vars1    = _mm256_madd_epi16  (dev1,  dev1);

    __m256i varsum   = _mm256_add_epi32   (vars0, vars1);
            varsum   = hsum_epi32_to_epi64(varsum);
            accum    = _mm256_add_epi64   (accum, varsum);
  }
  __m256i accum2 = _mm256_permute4x64_epi64(accum,  _MM_SHUFFLE(1, 0, 3, 2));
  __m256i accum3 = _mm256_add_epi64        (accum,  accum2);
  __m256i accum4 = _mm256_permute4x64_epi64(accum3, _MM_SHUFFLE(2, 3, 1, 0));
  __m256i v_tot  = _mm256_add_epi64        (accum3, accum4);
  __m128i vt128  = _mm256_castsi256_si128  (v_tot);

  _mm_storel_epi64((__m128i *)&vars, vt128);

  return (float)vars * varsum_to_f;
}

#else // INACCURATE_VARIANCE_CALCULATION

static double pixel_var_avx2(const uint8_t *buf, const uint32_t len)
{
  return pixel_var_avx2_largebuf(buf, len);
}

#endif // !INACCURATE_VARIANCE_CALCULATION

#endif // KVZ_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2

int kvz_strategy_register_picture_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
#if KVZ_BIT_DEPTH == 8
  // We don't actually use SAD for intra right now, other than 4x4 for
  // transform skip, but we might again one day and this is some of the
  // simplest code to look at for anyone interested in doing more
  // optimizations, so it's worth it to keep this maintained.
  if (bitdepth == 8){

    success &= kvz_strategyselector_register(opaque, "reg_sad", "avx2", 40, &kvz_reg_sad_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "avx2", 40, &sad_8bit_8x8_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "avx2", 40, &sad_8bit_16x16_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "avx2", 40, &sad_8bit_32x32_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "avx2", 40, &sad_8bit_64x64_avx2);

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "avx2", 40, &satd_4x4_8bit_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "avx2", 40, &satd_8x8_8bit_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "avx2", 40, &satd_16x16_8bit_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "avx2", 40, &satd_32x32_8bit_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "avx2", 40, &satd_64x64_8bit_avx2);

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "avx2", 40, &satd_8bit_4x4_dual_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "avx2", 40, &satd_8bit_8x8_dual_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "avx2", 40, &satd_8bit_16x16_dual_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "avx2", 40, &satd_8bit_32x32_dual_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "avx2", 40, &satd_8bit_64x64_dual_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "avx2", 40, &satd_any_size_8bit_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "avx2", 40, &satd_any_size_quad_avx2);

    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "avx2", 40, &pixels_calc_ssd_avx2);
    success &= kvz_strategyselector_register(opaque, "bipred_average", "avx2", 40, &bipred_average_avx2);
    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "avx2", 40, &get_optimized_sad_avx2);
    success &= kvz_strategyselector_register(opaque, "ver_sad", "avx2", 40, &ver_sad_avx2);
    success &= kvz_strategyselector_register(opaque, "hor_sad", "avx2", 40, &hor_sad_avx2);

    success &= kvz_strategyselector_register(opaque, "pixel_var", "avx2", 40, &pixel_var_avx2);

  }
#endif // KVZ_BIT_DEPTH == 8
#endif
  return success;
}
