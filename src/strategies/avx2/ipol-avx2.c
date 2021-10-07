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

#include "strategies/avx2/ipol-avx2.h"

#if COMPILE_INTEL_AVX2 && defined X86_64
#include <immintrin.h>
#include <stdio.h>
#include <string.h>

#include "encoder.h"
#include "kvazaar.h"
#include "search_inter.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/strategies-ipol.h"
#include "strategyselector.h"
#include "strategies/generic/ipol-generic.h"


extern int8_t kvz_g_luma_filter[4][8];
extern int8_t kvz_g_chroma_filter[8][4];

static int32_t kvz_eight_tap_filter_hor_avx2(int8_t *filter, kvz_pixel *data)
{
  __m128i fir = _mm_loadl_epi64((__m128i*)filter);
  __m128i row = _mm_loadl_epi64((__m128i*)data);
  __m128i acc;
  acc = _mm_maddubs_epi16(row, fir);
  __m128i temp = _mm_srli_si128(acc, 4);
  acc = _mm_add_epi16(acc, temp);
  temp = _mm_srli_si128(acc, 2);
  acc = _mm_add_epi16(acc, temp);
  int32_t filtered = _mm_cvtsi128_si32(acc);

  return filtered;
}

static int32_t kvz_eight_tap_filter_hor_16bit_avx2(int8_t *filter, int16_t *data)
{
  __m128i fir = _mm_loadl_epi64((__m128i*)filter);
  fir = _mm_cvtepi8_epi16(fir);
  __m128i row = _mm_loadu_si128((__m128i*)data);
  __m128i acc;
  acc = _mm_madd_epi16(fir, row);
  __m128i temp = _mm_srli_si128(acc, 8);
  acc = _mm_add_epi32(acc, temp);
  temp = _mm_srli_si128(acc, 4);
  acc = _mm_add_epi32(acc, temp);
  int32_t filtered = _mm_cvtsi128_si32(acc);

  return filtered;
}

static void kvz_eight_tap_filter_ver_16bit_1x8_avx2(int8_t *filter, int16_t *data, int16_t stride, kvz_pixel *out)
{
  // Interpolation filter shifts
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  // Filter weights
  __m256i all_taps = _mm256_castsi128_si256(_mm_cvtepi8_epi16(_mm_loadl_epi64((__m128i*)filter)));
  __m256i taps_01_23 = _mm256_shuffle_epi32(all_taps, _MM_SHUFFLE(0, 0, 0, 0));
  __m128i taps_23 = _mm_shuffle_epi32(_mm256_castsi256_si128(all_taps), _MM_SHUFFLE(1, 1, 1, 1));
  __m256i taps_45_67 = _mm256_shuffle_epi32(all_taps, _MM_SHUFFLE(2, 2, 2, 2));
  __m128i taps_67 = _mm_shuffle_epi32(_mm256_castsi256_si128(all_taps), _MM_SHUFFLE(3, 3, 3, 3));

  taps_01_23 = _mm256_inserti128_si256(taps_01_23, taps_23, 1);
  taps_45_67 = _mm256_inserti128_si256(taps_45_67, taps_67, 1);

  __m256i rows02 = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&data[0 * stride]));
  __m128i row2 = _mm_loadu_si128((__m128i*)&data[2 * stride]);
  rows02 = _mm256_inserti128_si256(rows02, row2, 1);

  __m256i rows13 = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&data[1 * stride]));
  __m128i row3 = _mm_loadu_si128((__m128i*)&data[3 * stride]);
  rows13 = _mm256_inserti128_si256(rows13, row3, 1);

  __m256i pairs_01_23_lo = _mm256_unpacklo_epi16(rows02, rows13);
  __m256i pairs_01_23_hi = _mm256_unpackhi_epi16(rows02, rows13);
  __m256i temp_01_23_lo = _mm256_madd_epi16(pairs_01_23_lo, taps_01_23);
  __m256i temp_01_23_hi = _mm256_madd_epi16(pairs_01_23_hi, taps_01_23);

  __m256i rows46 = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&data[4 * stride]));
  __m128i row6 = _mm_loadu_si128((__m128i*)&data[6 * stride]);
  rows46 = _mm256_inserti128_si256(rows46, row6, 1);

  __m256i rows57 = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&data[5 * stride]));
  __m128i row7 = _mm_loadu_si128((__m128i*)&data[7 * stride]);
  rows57 = _mm256_inserti128_si256(rows57, row7, 1);

  __m256i pairs_45_67_lo = _mm256_unpacklo_epi16(rows46, rows57);
  __m256i pairs_45_67_hi = _mm256_unpackhi_epi16(rows46, rows57);
  __m256i temp_45_67_lo = _mm256_madd_epi16(pairs_45_67_lo, taps_45_67);
  __m256i temp_45_67_hi = _mm256_madd_epi16(pairs_45_67_hi, taps_45_67);

  __m256i sum_lo_half = _mm256_add_epi32(temp_01_23_lo, temp_45_67_lo);
  __m256i sum_hi_half = _mm256_add_epi32(temp_01_23_hi, temp_45_67_hi);

  __m128i sum_lo = _mm_add_epi32(_mm256_castsi256_si128(sum_lo_half), _mm256_extracti128_si256(sum_lo_half, 1));
  __m128i sum_hi = _mm_add_epi32(_mm256_castsi256_si128(sum_hi_half), _mm256_extracti128_si256(sum_hi_half, 1));

  sum_lo = _mm_srai_epi32(sum_lo, shift2);
  sum_hi = _mm_srai_epi32(sum_hi, shift2);

  __m128i offset = _mm_set1_epi32(wp_offset1);
  sum_lo = _mm_add_epi32(sum_lo, offset);
  sum_lo = _mm_srai_epi32(sum_lo, wp_shift1);
  sum_hi = _mm_add_epi32(sum_hi, offset);
  sum_hi = _mm_srai_epi32(sum_hi, wp_shift1);
  __m128i filtered = _mm_packus_epi32(sum_lo, sum_hi);
  filtered = _mm_packus_epi16(filtered, filtered);


  _mm_storel_epi64((__m128i*)out, filtered);
}

static void kvz_ipol_8tap_hor_px_im_avx2(int8_t *filter,
  int width,
  int height,
  kvz_pixel *src,
  int16_t src_stride,
  int16_t *dst,
  int16_t dst_stride) {
  __m256i shuf01 = _mm256_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8,
                                    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8);
  __m256i shuf23 = _mm256_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
                                    2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10);
  __m256i shuf45 = _mm256_setr_epi8(4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
                                    4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12);
  __m256i shuf67 = _mm256_setr_epi8(6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
                                    6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14);

  __m256i all_w01 = _mm256_set1_epi16(*(uint16_t *)(filter + 0));
  __m256i all_w23 = _mm256_set1_epi16(*(uint16_t *)(filter + 2));
  __m256i all_w45 = _mm256_set1_epi16(*(uint16_t *)(filter + 4));
  __m256i all_w67 = _mm256_set1_epi16(*(uint16_t *)(filter + 6));

  int y_offset = -KVZ_LUMA_FILTER_OFFSET;
  int x_offset = -KVZ_LUMA_FILTER_OFFSET;

  kvz_pixel *top_left = src + src_stride * y_offset + x_offset;

  int y = 0;
  int x = 0;

  for (y = 0; y < height + KVZ_EXT_PADDING_LUMA; y += 2) {

    for (x = 0; x + 7 < width; x += 8) {

      kvz_pixel *chunk_ptr = top_left + src_stride * y + x;
      __m128i r0 = _mm_loadu_si128((__m128i*)(chunk_ptr + 0 * src_stride));
      __m128i r1 = _mm_loadu_si128((__m128i*)(chunk_ptr + 1 * src_stride));
      __m256i r0_r1 = _mm256_castsi128_si256(r0);
      r0_r1 = _mm256_inserti128_si256(r0_r1, r1, 1);

      __m256i r0_r1_01 = _mm256_shuffle_epi8(r0_r1, shuf01);
      __m256i r0_r1_23 = _mm256_shuffle_epi8(r0_r1, shuf23);
      __m256i r0_r1_45 = _mm256_shuffle_epi8(r0_r1, shuf45);
      __m256i r0_r1_67 = _mm256_shuffle_epi8(r0_r1, shuf67);

      __m256i dot01 = _mm256_maddubs_epi16(r0_r1_01, all_w01);
      __m256i dot23 = _mm256_maddubs_epi16(r0_r1_23, all_w23);
      __m256i dot45 = _mm256_maddubs_epi16(r0_r1_45, all_w45);
      __m256i dot67 = _mm256_maddubs_epi16(r0_r1_67, all_w67);

      __m256i sum0123 = _mm256_add_epi16(dot01, dot23);
      __m256i sum4567 = _mm256_add_epi16(dot45, dot67);
      __m256i sum = _mm256_add_epi16(sum0123, sum4567);

      __m128i *dst_r0 = (__m128i*)(dst + (y + 0) * dst_stride + x);
      __m128i *dst_r1 = (__m128i*)(dst + (y + 1) * dst_stride + x);
      __m128i sum_r0 = _mm256_castsi256_si128(sum);
      __m128i sum_r1 = _mm256_extracti128_si256(sum, 1);
      _mm_storeu_si128(dst_r0, sum_r0);
      _mm_storeu_si128(dst_r1, sum_r1);
    }
  }

  if (x < width) {
    for (int y = 0; y < height + KVZ_EXT_PADDING_LUMA; y += 2) {

      kvz_pixel *chunk_ptr = top_left + src_stride * y + x;
      __m128i r0 = _mm_loadu_si128((__m128i *)(chunk_ptr + 0 * src_stride));
      __m128i r1 = _mm_loadu_si128((__m128i *)(chunk_ptr + 1 * src_stride));
      __m256i r0_r1 = _mm256_castsi128_si256(r0);
      r0_r1 = _mm256_inserti128_si256(r0_r1, r1, 1);

      __m256i r0_r1_01 = _mm256_shuffle_epi8(r0_r1, shuf01);
      __m256i r0_r1_23 = _mm256_shuffle_epi8(r0_r1, shuf23);
      __m256i r0_r1_45 = _mm256_shuffle_epi8(r0_r1, shuf45);
      __m256i r0_r1_67 = _mm256_shuffle_epi8(r0_r1, shuf67);

      __m256i dot01 = _mm256_maddubs_epi16(r0_r1_01, all_w01);
      __m256i dot23 = _mm256_maddubs_epi16(r0_r1_23, all_w23);
      __m256i dot45 = _mm256_maddubs_epi16(r0_r1_45, all_w45);
      __m256i dot67 = _mm256_maddubs_epi16(r0_r1_67, all_w67);

      __m256i sum0123 = _mm256_add_epi16(dot01, dot23);
      __m256i sum4567 = _mm256_add_epi16(dot45, dot67);
      __m256i sum = _mm256_add_epi16(sum0123, sum4567);

      __m128i *dst_r0 = (__m128i*)(dst + (y + 0) * dst_stride + x);
      __m128i *dst_r1 = (__m128i*)(dst + (y + 1) * dst_stride + x);
      __m128i sum_r0 = _mm256_castsi256_si128(sum);
      __m128i sum_r1 = _mm256_extracti128_si256(sum, 1);
      _mm_storel_epi64(dst_r0, sum_r0);
      _mm_storel_epi64(dst_r1, sum_r1);
    }
  }
}

static void kvz_ipol_8tap_ver_im_px_avx2(int8_t *filter,
  int width,
  int height,
  int16_t *src,
  int16_t src_stride,
  kvz_pixel *dst,
  int16_t dst_stride)
{
  // Interpolation filter shifts
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  __m128i weights_8b = _mm_set1_epi64x(*(uint64_t *)filter);
  __m256i weights_16b = _mm256_cvtepi8_epi16(weights_8b);
  __m256i all_w01 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i all_w23 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(1, 1, 1, 1));
  __m256i all_w45 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(2, 2, 2, 2));
  __m256i all_w67 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(3, 3, 3, 3));

  for (int x = 0; x + 3 < width; x += 4) {

    int16_t *strip_ptr = src + 0 * src_stride + x;

    // Initial values
    // Broadcasted rows in both lanes
    // __m256i r0; // Unused
    // __m256i r1; // Unused
    __m256i r2 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 0 * src_stride));
    __m256i r3 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 1 * src_stride));
    __m256i r4 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 2 * src_stride));
    __m256i r5 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 3 * src_stride));
    __m256i r6 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 4 * src_stride));
    __m256i r7 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 5 * src_stride));
    __m256i r8 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 6 * src_stride));

    // Consecutive rows in low and high lanes
    // __m256i r0_r1; // Unused
    // __m256i r1_r2; // Unused
    __m256i r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
    __m256i r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);
    __m256i r4_r5 = _mm256_blend_epi32(r4, r5, 0xF0);
    __m256i r5_r6 = _mm256_blend_epi32(r5, r6, 0xF0);
    __m256i r6_r7 = _mm256_blend_epi32(r6, r7, 0xF0);
    __m256i r7_r8 = _mm256_blend_epi32(r7, r8, 0xF0);

    // Paired samples of consecutive rows
    __m256i r01_r12;
    __m256i r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);
    __m256i r45_r56 = _mm256_unpacklo_epi16(r4_r5, r5_r6);
    __m256i r67_r78 = _mm256_unpacklo_epi16(r6_r7, r7_r8);

    for (int y = 0; y < height; y += 2) {

      strip_ptr = src + y * src_stride + x;

      // Slide window
      r01_r12 = r23_r34;
      r23_r34 = r45_r56;
      r45_r56 = r67_r78;
      r6 = r8;
      r7 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 7 * src_stride));
      r8 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 8 * src_stride));
      r6_r7 = _mm256_blend_epi32(r6, r7, 0xF0);
      r7_r8 = _mm256_blend_epi32(r7, r8, 0xF0);

      r67_r78 = _mm256_unpacklo_epi16(r6_r7, r7_r8);

      __m256i dot01 = _mm256_madd_epi16(r01_r12, all_w01);
      __m256i dot23 = _mm256_madd_epi16(r23_r34, all_w23);
      __m256i dot45 = _mm256_madd_epi16(r45_r56, all_w45);
      __m256i dot67 = _mm256_madd_epi16(r67_r78, all_w67);

      __m256i sum0123 = _mm256_add_epi32(dot01, dot23);
      __m256i sum4567 = _mm256_add_epi32(dot45, dot67);
      __m256i sum = _mm256_add_epi32(sum0123, sum4567);
      sum = _mm256_srai_epi32(sum, shift2);
      sum = _mm256_add_epi32(sum, _mm256_set1_epi32(wp_offset1));
      sum = _mm256_srai_epi32(sum, wp_shift1);
      sum = _mm256_packs_epi32(sum, sum);
      sum = _mm256_packus_epi16(sum, sum);

      kvz_pixel *dst_addr0 = &dst[(y + 0) * dst_stride + x];
      kvz_pixel *dst_addr1 = &dst[(y + 1) * dst_stride + x];
      *(uint32_t*)dst_addr0 = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum));
      *(uint32_t*)dst_addr1 = _mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 1));
    }
  }
}

static void kvz_ipol_8tap_ver_im_hi_avx2(int8_t *filter,
int width,
int height,
int16_t *src,
int16_t src_stride,
int16_t *dst,
int16_t dst_stride)
{
  const int shift2 = 6;

  __m128i weights_8b = _mm_set1_epi64x(*(uint64_t *)filter);
  __m256i weights_16b = _mm256_cvtepi8_epi16(weights_8b);
  __m256i all_w01 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i all_w23 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(1, 1, 1, 1));
  __m256i all_w45 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(2, 2, 2, 2));
  __m256i all_w67 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(3, 3, 3, 3));

  for (int x = 0; x + 3 < width; x += 4) {

    int16_t *strip_ptr = src + 0 * src_stride + x;

    // Initial values
    // Broadcasted rows in both lanes
    // __m256i r0; // Unused
    // __m256i r1; // Unused
    __m256i r2 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 0 * src_stride));
    __m256i r3 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 1 * src_stride));
    __m256i r4 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 2 * src_stride));
    __m256i r5 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 3 * src_stride));
    __m256i r6 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 4 * src_stride));
    __m256i r7 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 5 * src_stride));
    __m256i r8 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 6 * src_stride));

    // Consecutive rows in low and high lanes
    // __m256i r0_r1; // Unused
    // __m256i r1_r2; // Unused
    __m256i r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
    __m256i r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);
    __m256i r4_r5 = _mm256_blend_epi32(r4, r5, 0xF0);
    __m256i r5_r6 = _mm256_blend_epi32(r5, r6, 0xF0);
    __m256i r6_r7 = _mm256_blend_epi32(r6, r7, 0xF0);
    __m256i r7_r8 = _mm256_blend_epi32(r7, r8, 0xF0);

    // Paired samples of consecutive rows
    __m256i r01_r12;
    __m256i r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);
    __m256i r45_r56 = _mm256_unpacklo_epi16(r4_r5, r5_r6);
    __m256i r67_r78 = _mm256_unpacklo_epi16(r6_r7, r7_r8);

    for (int y = 0; y < height; y += 2) {

      strip_ptr = src + y * src_stride + x;

      // Slide window
      r01_r12 = r23_r34;
      r23_r34 = r45_r56;
      r45_r56 = r67_r78;
      r6 = r8;
      r7 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 7 * src_stride));
      r8 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 8 * src_stride));
      r6_r7 = _mm256_blend_epi32(r6, r7, 0xF0);
      r7_r8 = _mm256_blend_epi32(r7, r8, 0xF0);

      r67_r78 = _mm256_unpacklo_epi16(r6_r7, r7_r8);

      __m256i dot01 = _mm256_madd_epi16(r01_r12, all_w01);
      __m256i dot23 = _mm256_madd_epi16(r23_r34, all_w23);
      __m256i dot45 = _mm256_madd_epi16(r45_r56, all_w45);
      __m256i dot67 = _mm256_madd_epi16(r67_r78, all_w67);

      __m256i sum0123 = _mm256_add_epi32(dot01, dot23);
      __m256i sum4567 = _mm256_add_epi32(dot45, dot67);
      __m256i sum = _mm256_add_epi32(sum0123, sum4567);
      sum = _mm256_srai_epi32(sum, shift2); // TODO: -8192 offsetting for extreme values
      sum = _mm256_packs_epi32(sum, sum);

      int16_t *dst_addr0 = &dst[(y + 0) * dst_stride + x];
      int16_t *dst_addr1 = &dst[(y + 1) * dst_stride + x];
      _mm_storel_epi64((__m128i*)dst_addr0, _mm256_castsi256_si128(sum));
      _mm_storel_epi64((__m128i*)dst_addr1, _mm256_extracti128_si256(sum, 1));
    }
  }
}

static void kvz_ipol_4tap_hor_px_im_avx2(int8_t *filter,
  int width,
  int height,
  kvz_pixel *src,
  int16_t src_stride,
  int16_t *dst,
  int16_t dst_stride) {

  __m256i shuf01 = _mm256_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4,
    8, 9, 9, 10, 10, 11, 11, 12,
    0, 1, 1, 2, 2, 3, 3, 4,
    8, 9, 9, 10, 10, 11, 11, 12);

  __m256i shuf23 = _mm256_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6,
    10, 11, 11, 12, 12, 13, 13, 14,
    2, 3, 3, 4, 4, 5, 5, 6,
    10, 11, 11, 12, 12, 13, 13, 14);

  __m256i all_w01 = _mm256_set1_epi16(*(uint16_t*)(filter + 0));
  __m256i all_w23 = _mm256_set1_epi16(*(uint16_t*)(filter + 2));

  int y_offset = -KVZ_CHROMA_FILTER_OFFSET;
  int x_offset = -KVZ_CHROMA_FILTER_OFFSET;

  kvz_pixel *top_left = src + src_stride * y_offset + x_offset;

  int y = 0;
  int x = 0;

  for (y = 0; y < height + KVZ_EXT_PADDING_CHROMA; y += 4) {

    for (x = 0; x + 3 < width; x += 4) {

      kvz_pixel *chunk_ptr = top_left + src_stride * y + x;
      __m128i r0r1 = _mm_loadl_epi64((__m128i*)(chunk_ptr + 0 * src_stride));
      __m128i r2r3 = _mm_loadl_epi64((__m128i*)(chunk_ptr + 2 * src_stride));
      r0r1 = _mm_insert_epi64(r0r1, *(uint64_t*)(chunk_ptr + 1 * src_stride), 1);
      r2r3 = _mm_insert_epi64(r2r3, *(uint64_t*)(chunk_ptr + 3 * src_stride), 1);

      __m256i r0r1_r2r3 = _mm256_castsi128_si256(r0r1);
      r0r1_r2r3 = _mm256_inserti128_si256(r0r1_r2r3, r2r3, 1);

      __m256i r0_r1_01 = _mm256_shuffle_epi8(r0r1_r2r3, shuf01);
      __m256i r0_r1_23 = _mm256_shuffle_epi8(r0r1_r2r3, shuf23);

      __m256i dot01 = _mm256_maddubs_epi16(r0_r1_01, all_w01);
      __m256i dot23 = _mm256_maddubs_epi16(r0_r1_23, all_w23);

      __m256i sum = _mm256_add_epi16(dot01, dot23);

      __m128i *dst_r0 = (__m128i*)(dst + (y + 0) * dst_stride + x);
      __m128i *dst_r1 = (__m128i*)(dst + (y + 1) * dst_stride + x);
      __m128i *dst_r2 = (__m128i*)(dst + (y + 2) * dst_stride + x);
      __m128i *dst_r3 = (__m128i*)(dst + (y + 3) * dst_stride + x);
      __m128i sum_r0r1 = _mm256_castsi256_si128(sum);
      __m128i sum_r2r3 = _mm256_extracti128_si256(sum, 1);
      _mm_storel_epi64(dst_r0, sum_r0r1);
      _mm_storeh_pd((double*)dst_r1, _mm_castsi128_pd(sum_r0r1));
      _mm_storel_epi64(dst_r2, sum_r2r3);
      _mm_storeh_pd((double*)dst_r3, _mm_castsi128_pd(sum_r2r3));
    }
  }
}

static void kvz_ipol_4tap_ver_im_px_avx2(int8_t *filter,
  int width,
  int height,
  int16_t *src,
  int16_t src_stride,
  kvz_pixel *dst,
  int16_t dst_stride)
{
  // Interpolation filter shifts
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  __m128i weights_8b = _mm_set1_epi64x(*(uint64_t*)filter);
  __m256i weights_16b = _mm256_cvtepi8_epi16(weights_8b);
  __m256i all_w01 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i all_w23 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(1, 1, 1, 1));

  for (int x = 0; x + 3 < width; x += 4) {

    int16_t *strip_ptr = src + 0 * src_stride + x;

    // Initial values
    // Broadcasted rows in both lanes
    // __m256i r0; // Unused
    // __m256i r1; // Unused
    __m256i r2 = _mm256_set1_epi64x(*(uint64_t*)(strip_ptr + 0 * src_stride));
    __m256i r3 = _mm256_set1_epi64x(*(uint64_t*)(strip_ptr + 1 * src_stride));
    __m256i r4 = _mm256_set1_epi64x(*(uint64_t*)(strip_ptr + 2 * src_stride));

    // Consecutive rows in low and high lanes
    // __m256i r0_r1; // Unused
    // __m256i r1_r2; // Unused
    __m256i r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
    __m256i r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);

    // Paired samples of consecutive rows
    __m256i r01_r12;
    __m256i r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);

    for (int y = 0; y < height; y += 2) {

      strip_ptr = src + y * src_stride + x;

      // Slide window
      r01_r12 = r23_r34;
      r2 = r4;
      r3 = _mm256_set1_epi64x(*(uint64_t*)(strip_ptr + 3 * src_stride));
      r4 = _mm256_set1_epi64x(*(uint64_t*)(strip_ptr + 4 * src_stride));
      r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
      r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);

      r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);

      __m256i dot01 = _mm256_madd_epi16(r01_r12, all_w01);
      __m256i dot23 = _mm256_madd_epi16(r23_r34, all_w23);

      __m256i sum = _mm256_add_epi32(dot01, dot23);
      sum = _mm256_srai_epi32(sum, shift2);
      sum = _mm256_add_epi32(sum, _mm256_set1_epi32(wp_offset1));
      sum = _mm256_srai_epi32(sum, wp_shift1);
      sum = _mm256_packs_epi32(sum, sum);
      sum = _mm256_packus_epi16(sum, sum);

      kvz_pixel *dst_addr0 = &dst[(y + 0) * dst_stride + x];
      kvz_pixel *dst_addr1 = &dst[(y + 1) * dst_stride + x];
      *(uint32_t*)dst_addr0 = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum));
      *(uint32_t*)dst_addr1 = _mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 1));
    }
  }
}

static void kvz_ipol_4tap_ver_im_hi_avx2(int8_t *filter,
  int width,
  int height,
  int16_t *src,
  int16_t src_stride,
  int16_t *dst,
  int16_t dst_stride)
{
  const int shift2 = 6;

  __m128i weights_8b = _mm_set1_epi64x(*(uint64_t *)filter);
  __m256i weights_16b = _mm256_cvtepi8_epi16(weights_8b);
  __m256i all_w01 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(0, 0, 0, 0));
  __m256i all_w23 = _mm256_shuffle_epi32(weights_16b, _MM_SHUFFLE(1, 1, 1, 1));

  for (int x = 0; x + 3 < width; x += 4) {

    int16_t *strip_ptr = src + 0 * src_stride + x;

    // Initial values
    // Broadcasted rows in both lanes
    // __m256i r0; // Unused
    // __m256i r1; // Unused
    __m256i r2 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 0 * src_stride));
    __m256i r3 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 1 * src_stride));
    __m256i r4 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 2 * src_stride));

    // Consecutive rows in low and high lanes
    // __m256i r0_r1; // Unused
    // __m256i r1_r2; // Unused
    __m256i r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
    __m256i r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);

    // Paired samples of consecutive rows
    __m256i r01_r12;
    __m256i r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);

    for (int y = 0; y < height; y += 2) {

      strip_ptr = src + y * src_stride + x;

      // Slide window
      r01_r12 = r23_r34;
      r2 = r4;
      r3 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 3 * src_stride));
      r4 = _mm256_set1_epi64x(*(uint64_t *)(strip_ptr + 4 * src_stride));
      r2_r3 = _mm256_blend_epi32(r2, r3, 0xF0);
      r3_r4 = _mm256_blend_epi32(r3, r4, 0xF0);

      r23_r34 = _mm256_unpacklo_epi16(r2_r3, r3_r4);

      __m256i dot01 = _mm256_madd_epi16(r01_r12, all_w01);
      __m256i dot23 = _mm256_madd_epi16(r23_r34, all_w23);

      __m256i sum = _mm256_add_epi32(dot01, dot23);
      sum = _mm256_srai_epi32(sum, shift2);
      sum = _mm256_packs_epi32(sum, sum);

      int16_t *dst_addr0 = &dst[(y + 0) * dst_stride + x];
      int16_t *dst_addr1 = &dst[(y + 1) * dst_stride + x];
      _mm_storel_epi64((__m128i *)dst_addr0, _mm256_castsi256_si128(sum));
      _mm_storel_epi64((__m128i *)dst_addr1, _mm256_extracti128_si256(sum, 1));
    }
  }
}

static void kvz_filter_hpel_blocks_hor_ver_luma_avx2(const encoder_control_t * encoder,
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
  int x, y, first_y;

  // Interpolation filter shifts
  int16_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  int8_t *fir0 = kvz_g_luma_filter[0];
  int8_t *fir2 = kvz_g_luma_filter[2];

  int16_t dst_stride = LCU_WIDTH;
  int16_t hor_stride = LCU_WIDTH;

  int16_t *hor_pos0 = hor_intermediate[0];
  int16_t *hor_pos2 = hor_intermediate[1];
  int16_t *col_pos0 = hor_first_cols[0];
  int16_t *col_pos2 = hor_first_cols[2];

  // Horizontally filtered samples from the top row are
  // not needed unless samples for diagonal positions are filtered later.
  first_y = fme_level > 1 ? 0 : 1;

  // HORIZONTAL STEP
  // Integer pixels
  for (y = 0; y < height + KVZ_EXT_PADDING_LUMA + 1; ++y) {
    
    for (x = 0; x + 7 < width; x += 8) {
      int ypos = y - KVZ_LUMA_FILTER_OFFSET;
      int xpos = x + 1;
      __m128i* out = (__m128i*)&hor_pos0[y * hor_stride + x];
      __m128i chunk = _mm_loadl_epi64((__m128i*)&src[src_stride*ypos + xpos]);
      chunk = _mm_cvtepu8_epi16(chunk);
      chunk = _mm_slli_epi16(chunk, 6); // Multiply by 64
      _mm_storeu_si128(out, chunk); //TODO: >> shift1
    }
  }

  // Write the first column in contiguous memory
  x = 0;
  for (y = 0; y < height + KVZ_EXT_PADDING_LUMA + 1; ++y) {
    int ypos = y - KVZ_LUMA_FILTER_OFFSET;
    int32_t first_sample = src[src_stride*ypos + x] << 6 >> shift1;
    col_pos0[y] = first_sample;
  }

  // Half pixels
  kvz_ipol_8tap_hor_px_im_avx2(fir2, width, height + 1, src + 1, src_stride, hor_pos2, hor_stride);

  // Write the first column in contiguous memory
  x = 0;
  for (y = first_y; y < height + KVZ_EXT_PADDING_LUMA + 1; ++y) {
    int ypos = y - KVZ_LUMA_FILTER_OFFSET;
    int xpos = x - KVZ_LUMA_FILTER_OFFSET;
    col_pos2[y] = kvz_eight_tap_filter_hor_avx2(fir2, &src[src_stride*ypos + xpos]) >> shift1;
  }

  // VERTICAL STEP
  kvz_pixel *out_l = filtered[0];
  kvz_pixel *out_r = filtered[1];
  kvz_pixel *out_t = filtered[2];
  kvz_pixel *out_b = filtered[3];

  // Right
  int16_t *im = &hor_pos2[hor_stride];
  kvz_pixel *dst = out_r;
  kvz_ipol_8tap_ver_im_px_avx2(fir0, width, height, im, hor_stride, dst, dst_stride);

  // Left
  // Copy from the right filtered block and filter the extra column
  for (y = 0; y < height; ++y) {
    x = 0;
    *(uint64_t*)&out_l[y * dst_stride + x] = *(uint64_t*)&out_r[y * dst_stride + x] << 8;
    for (x = 8; x < width; x += 8) *(int64_t*)&out_l[y * dst_stride + x] = *(int64_t*)&out_r[y * dst_stride + x - 1];
    x = 0;
    int16_t sample = 64 * col_pos2[y + 1 + KVZ_LUMA_FILTER_OFFSET] >> shift2;
    sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
    out_l[y * dst_stride + x] = sample;
  }

  // Top
  im = hor_pos0;
  dst = out_t;
  kvz_ipol_8tap_ver_im_px_avx2(fir2, width, height, im, hor_stride, dst, dst_stride);

  // Bottom
  // Copy what can be copied from the top filtered values.
  // Then filter the last row from horizontal intermediate buffer.
  for (y = 0; y < height - 1; ++y) {
    for (x = 0; x + 7 < width; x += 8) {
      *(int64_t*)&out_b[(y + 0) * dst_stride + x] = *(int64_t*)&out_t[(y + 1) * dst_stride + x];
    }
  }

  for (x = 0; x + 7 < width; x += 8) {
    kvz_eight_tap_filter_ver_16bit_1x8_avx2(fir2, &hor_pos0[(y + 1) * hor_stride + x], hor_stride, &out_b[y * dst_stride + x]);
  }
}

static void kvz_filter_hpel_blocks_diag_luma_avx2(const encoder_control_t * encoder,
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
  int x, y;

  // Interpolation filter shifts
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  int8_t *fir2 = kvz_g_luma_filter[2];

  int16_t dst_stride = LCU_WIDTH;
  int16_t hor_stride = LCU_WIDTH;

  int16_t *hor_pos2 = hor_intermediate[1];
  int16_t *col_pos2 = hor_first_cols[2];

  // VERTICAL STEP
  kvz_pixel *out_tl = filtered[0];
  kvz_pixel *out_tr = filtered[1];
  kvz_pixel *out_bl = filtered[2];
  kvz_pixel *out_br = filtered[3];

  // Top-Right
  int16_t *im = hor_pos2;
  kvz_pixel *dst = out_tr;
  kvz_ipol_8tap_ver_im_px_avx2(fir2, width, height, im, hor_stride, dst, dst_stride);

  // Top-left
  // Copy from the top-right filtered block and filter the extra column
  for (y = 0; y < height; ++y) {
    x = 0;
    int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(fir2, &col_pos2[y]) >> shift2;
    sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
    out_tl[y * dst_stride + x] = sample;

    for (x = 1; x < width; ++x) out_tl[y * dst_stride + x] = out_tr[y * dst_stride + x - 1];
  }

  // Bottom-right
  // Copy what can be copied from top-right filtered values. Filter the last row.
  for (y = 0; y < height - 1; ++y) {
    for (x = 0; x + 7 < width; x += 8) {
      memcpy(&out_br[y * dst_stride + x], &out_tr[(y + 1) * dst_stride + x], 8);
    }
  }

  for (x = 0; x + 7 < width; x += 8) {
    kvz_eight_tap_filter_ver_16bit_1x8_avx2(fir2, &hor_pos2[(y + 1) * hor_stride + x], hor_stride, &out_br[y * dst_stride + x]);
  }

  // Bottom-left
  // Copy what can be copied from the top-left filtered values.
  // Copy what can be copied from the bottom-right filtered values.
  // Finally filter the last pixel from the column array.
  for (y = 0; y < height - 1; ++y) {
    for (x = 0; x + 7 < width; x += 8) {
      memcpy(&out_bl[y * dst_stride + x], &out_tl[(y + 1) * dst_stride + x], 8);
    }
  }

  for (x = 1; x < width; ++x) out_bl[y * dst_stride + x] = out_br[y * dst_stride + x - 1];
  x = 0;
  int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(fir2, &col_pos2[(y + 1)]) >> shift2;
  sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
  out_bl[y * dst_stride + x] = sample;
}

static void kvz_filter_qpel_blocks_hor_ver_luma_avx2(const encoder_control_t * encoder,
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
  int x, y;

  // Interpolation filter shifts
  int16_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  int8_t *fir0 = kvz_g_luma_filter[0];
  int8_t *fir2 = kvz_g_luma_filter[2];
  int8_t *fir1 = kvz_g_luma_filter[1];
  int8_t *fir3 = kvz_g_luma_filter[3];

  // Horiziontal positions. Positions 0 and 2 have already been calculated in filtered.
  int16_t *hor_pos0 = hor_intermediate[0];
  int16_t *hor_pos2 = hor_intermediate[1];
  int16_t *hor_pos_l = hor_intermediate[3];
  int16_t *hor_pos_r = hor_intermediate[4];
  int8_t *hor_fir_l = hpel_off_x != 0 ? fir1 : fir3;
  int8_t *hor_fir_r = hpel_off_x != 0 ? fir3 : fir1;
  int16_t *col_pos_l = hor_first_cols[1];
  int16_t *col_pos_r = hor_first_cols[3];

  int16_t dst_stride = LCU_WIDTH;
  int16_t hor_stride = LCU_WIDTH;

  int16_t *hor_hpel_pos = hpel_off_x != 0 ? hor_pos2 : hor_pos0;
  int16_t *col_pos_hor = hpel_off_x != 0 ? hor_first_cols[2] : hor_first_cols[0];

  // Specify if integer pixels are filtered from left or/and top integer samples
  int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;

  // HORIZONTAL STEP
  // Left QPEL
  int sample_off_y = hpel_off_y < 0 ? 0 : 1;
  kvz_ipol_8tap_hor_px_im_avx2(hor_fir_l, width, height + 1, src + 1, src_stride, hor_pos_l, hor_stride);

  // Write the first column in contiguous memory
  x = 0;
  for (y = 0; y < height + KVZ_EXT_PADDING_LUMA + 1; ++y) {
    int ypos = y - KVZ_LUMA_FILTER_OFFSET;
    int xpos = x - KVZ_LUMA_FILTER_OFFSET;
    col_pos_l[y] = kvz_eight_tap_filter_hor_avx2(hor_fir_l, &src[src_stride*ypos + xpos]) >> shift1;
  }

  // Right QPEL
  kvz_ipol_8tap_hor_px_im_avx2(hor_fir_r, width, height + 1, src + 1, src_stride, hor_pos_r, hor_stride);

  // Write the first column in contiguous memory
  x = 0;
  for (y = 0; y < height + KVZ_EXT_PADDING_LUMA + 1; ++y) {
    int ypos = y - KVZ_LUMA_FILTER_OFFSET;
    int xpos = x - KVZ_LUMA_FILTER_OFFSET;
    col_pos_r[y] = kvz_eight_tap_filter_hor_avx2(hor_fir_r, &src[src_stride*ypos + xpos]) >> shift1;
  }

  // VERTICAL STEP
  kvz_pixel *out_l = filtered[0];
  kvz_pixel *out_r = filtered[1];
  kvz_pixel *out_t = filtered[2];
  kvz_pixel *out_b = filtered[3];

  int8_t *ver_fir_l = hpel_off_y != 0 ? fir2 : fir0;
  int8_t *ver_fir_r = hpel_off_y != 0 ? fir2 : fir0;
  int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  // Left QPEL (1/4 or 3/4 x positions) 
  // Filter block and then filter column and align if neccessary
  int16_t *im = &hor_pos_l[sample_off_y * hor_stride];
  kvz_pixel *dst = out_l;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_l, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_l) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_l[y * dst_stride + x - 1];
        *(uint64_t*)&out_l[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_l, &col_pos_l[y + sample_off_y]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_l[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_l[y * dst_stride + x] = chunk;
    }
  }

  // Right QPEL (3/4 or 1/4 x positions)
  // Filter block and then filter column and align if neccessary
  im = &hor_pos_r[sample_off_y * hor_stride];
  dst = out_r;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_r, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_r) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_r[y * dst_stride + x - 1];
        *(uint64_t*)&out_r[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_r, &col_pos_r[y + sample_off_y]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_r[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_r[y * dst_stride + x] = chunk;
    }
  }

  
  // Top QPEL (1/4 or 3/4 y positions)
  // Filter block and then filter column and align if neccessary
  int sample_off_x = (hpel_off_x > -1 ? 1 : 0);

  im = &hor_hpel_pos[off_y_fir_t * hor_stride];
  dst = out_t;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride);

  if (!sample_off_x) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_t[y * dst_stride + x - 1];
        *(uint64_t*)&out_t[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_t, &col_pos_hor[y + off_y_fir_t]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_t[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_t[y * dst_stride + x] = chunk;
    }
  }

  // Bottom QPEL (3/4 or 1/4 y positions)
  // Filter block and then filter column and align if neccessary

  im = &hor_hpel_pos[off_y_fir_b * hor_stride];
  dst = out_b;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride);

  if (!sample_off_x) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_b[y * dst_stride + x - 1];
        *(uint64_t*)&out_b[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_b, &col_pos_hor[y + off_y_fir_b]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_b[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_b[y * dst_stride + x] = chunk;
    }
  }
}

static void kvz_filter_qpel_blocks_diag_luma_avx2(const encoder_control_t * encoder,
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
  int x, y;

  // Interpolation filter shifts
  int32_t shift2 = 6;

  // Weighted prediction offset and shift
  int32_t wp_shift1 = 14 - KVZ_BIT_DEPTH;
  int32_t wp_offset1 = 1 << (wp_shift1 - 1);

  int8_t *fir1 = kvz_g_luma_filter[1];
  int8_t *fir3 = kvz_g_luma_filter[3];

  int16_t *hor_pos_l = hor_intermediate[3];
  int16_t *hor_pos_r = hor_intermediate[4];

  int16_t *col_pos_l = hor_first_cols[1];
  int16_t *col_pos_r = hor_first_cols[3];

  int16_t dst_stride = LCU_WIDTH;
  int16_t hor_stride = LCU_WIDTH;

  // VERTICAL STEP
  kvz_pixel *out_tl = filtered[0];
  kvz_pixel *out_tr = filtered[1];
  kvz_pixel *out_bl = filtered[2];
  kvz_pixel *out_br = filtered[3];

  int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  // Specify if integer pixels are filtered from left or/and top integer samples
  int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;

  // Top-left QPEL
  // Filter block and then filter column and align if neccessary
  int16_t *im = &hor_pos_l[off_y_fir_t * hor_stride];
  kvz_pixel *dst = out_tl;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_l) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_tl[y * dst_stride + x - 1];
        *(uint64_t*)&out_tl[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_t, &col_pos_l[y + off_y_fir_t]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_tl[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_tl[y * dst_stride + x] = chunk;
    }
  }

  // Top-right QPEL
  // Filter block and then filter column and align if neccessary

  im = &hor_pos_r[off_y_fir_t * hor_stride];
  dst = out_tr;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_r) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_tr[y * dst_stride + x - 1];
        *(uint64_t*)&out_tr[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_t, &col_pos_r[y + off_y_fir_t]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_tr[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_tr[y * dst_stride + x] = chunk;
    }
  }

  // Bottom-left QPEL
  // Filter block and then filter column and align if neccessary
  im = &hor_pos_l[off_y_fir_b * hor_stride];
  dst = out_bl;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_l) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_bl[y * dst_stride + x - 1];
        *(uint64_t*)&out_bl[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_b, &col_pos_l[y + off_y_fir_b]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_bl[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_bl[y * dst_stride + x] = chunk;
    }
  }

  // Bottom-right QPEL
  // Filter block and then filter column and align if neccessary
  im = &hor_pos_r[off_y_fir_b * hor_stride];
  dst = out_br;
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride);

  if (!off_x_fir_r) {
    for (y = 0; y < height; ++y) {
      for (x = width - 8; x >= 8; x -= 8) {
        uint64_t chunk = *(uint64_t*)&out_br[y * dst_stride + x - 1];
        *(uint64_t*)&out_br[y * dst_stride + x] = chunk;
      }

      x = 0;
      int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(ver_fir_b, &col_pos_r[y + off_y_fir_b]) >> shift2;
      sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
      uint64_t first = sample;
      uint64_t rest = *(uint64_t*)&out_br[y * dst_stride + x];
      uint64_t chunk = (rest << 8) | first;
      *(uint64_t*)&out_br[y * dst_stride + x] = chunk;
    }
  }
}

static void kvz_sample_quarterpel_luma_avx2(const encoder_control_t * const encoder,
  kvz_pixel *src, 
  int16_t src_stride, 
  int width, 
  int height, 
  kvz_pixel *dst, 
  int16_t dst_stride, 
  int8_t hor_flag, 
  int8_t ver_flag, 
  const int16_t mv[2])
{
  // TODO: horizontal and vertical only filtering
  int8_t *hor_fir = kvz_g_luma_filter[mv[0] & 3];
  int8_t *ver_fir = kvz_g_luma_filter[mv[1] & 3];

  // Buffer for intermediate values with one extra row 
  // because the loop writes two rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
  int16_t hor_stride = LCU_WIDTH;

  kvz_ipol_8tap_hor_px_im_avx2(hor_fir, width, height, src, src_stride, hor_intermediate, hor_stride);
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride);
}


static void kvz_sample_quarterpel_luma_hi_avx2(const encoder_control_t * const encoder,
  kvz_pixel *src, 
  int16_t src_stride, 
  int width, 
  int height, 
  int16_t *dst, 
  int16_t dst_stride, 
  int8_t hor_flag, 
  int8_t ver_flag, 
  const int16_t mv[2])
{
  // TODO: horizontal and vertical only filtering
  int8_t *hor_fir = kvz_g_luma_filter[mv[0] & 3];
  int8_t *ver_fir = kvz_g_luma_filter[mv[1] & 3];
  
  // Buffer for intermediate values with one extra row 
  // because the loop writes two rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
  int16_t hor_stride = LCU_WIDTH;

  kvz_ipol_8tap_hor_px_im_avx2(hor_fir, width, height, src, src_stride, hor_intermediate, hor_stride);
  kvz_ipol_8tap_ver_im_hi_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride);
}


static void kvz_sample_octpel_chroma_avx2(const encoder_control_t *const encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  kvz_pixel *dst,
  int16_t dst_stride,
  int8_t hor_flag,
  int8_t ver_flag,
  const int16_t mv[2])
{
  // TODO: Optimizations for rest of the blocks (for example 2x8).
  if (width % 4 != 0) {
    kvz_sample_octpel_chroma_generic(encoder, src, src_stride, width, height, dst, dst_stride, hor_flag, ver_flag, mv);
    return;
  }
  int8_t *hor_fir = kvz_g_chroma_filter[mv[0] & 7];
  int8_t *ver_fir = kvz_g_chroma_filter[mv[1] & 7];

  // Buffer for intermediate values with 3 extra rows 
  // because the loop writes four rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_CHROMA_SIMD];
  int16_t hor_stride = LCU_WIDTH_C;

  kvz_ipol_4tap_hor_px_im_avx2(hor_fir, width, height, src, src_stride, hor_intermediate, hor_stride);
  kvz_ipol_4tap_ver_im_px_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride);
}

static void kvz_sample_octpel_chroma_hi_avx2(const encoder_control_t *const encoder,
  kvz_pixel *src,
  int16_t src_stride,
  int width,
  int height,
  int16_t *dst,
  int16_t dst_stride,
  int8_t hor_flag,
  int8_t ver_flag,
  const int16_t mv[2])
{
  // TODO: Optimizations for rest of the blocks (for example 2x8).
  if (width % 4 != 0) {
    kvz_sample_octpel_chroma_hi_generic(encoder, src, src_stride, width, height, dst, dst_stride, hor_flag, ver_flag, mv);
    return;
  }
  int8_t *hor_fir = kvz_g_chroma_filter[mv[0] & 7];
  int8_t *ver_fir = kvz_g_chroma_filter[mv[1] & 7];

  // Buffer for intermediate values with 3 extra rows 
  // because the loop writes four rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_CHROMA_SIMD];
  int16_t hor_stride = LCU_WIDTH_C;

  kvz_ipol_4tap_hor_px_im_avx2(hor_fir, width, height, src, src_stride, hor_intermediate, hor_stride);
  kvz_ipol_4tap_ver_im_hi_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride);
}

#endif //COMPILE_INTEL_AVX2 && defined X86_64

int kvz_strategy_register_ipol_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2 && defined X86_64
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_hor_ver_luma", "avx2", 40, &kvz_filter_hpel_blocks_hor_ver_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_diag_luma", "avx2", 40, &kvz_filter_hpel_blocks_diag_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_hor_ver_luma", "avx2", 40, &kvz_filter_qpel_blocks_hor_ver_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_diag_luma", "avx2", 40, &kvz_filter_qpel_blocks_diag_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma", "avx2", 40, &kvz_sample_quarterpel_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma", "avx2", 40, &kvz_sample_octpel_chroma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma_hi", "avx2", 40, &kvz_sample_quarterpel_luma_hi_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma_hi", "avx2", 40, &kvz_sample_octpel_chroma_hi_avx2);
  }
#endif //COMPILE_INTEL_AVX2 && defined X86_64
  return success;
}
