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
  // Only used by the block filters: store sum >> 6 as int16 like generic.
  sum_lo = _mm_srai_epi32(_mm_slli_epi32(sum_lo, 16), 16);
  sum_hi = _mm_srai_epi32(_mm_slli_epi32(sum_hi, 16), 16);

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

  if (filter == kvz_g_luma_filter[0]) {
    // Filter row 0 is a single tap of 64 at position 3: the sample shifted
    // left by 6. Rows and store widths follow the general loops below.
    for (int yy = 0; yy < height + KVZ_EXT_PADDING_LUMA; yy += 2) {
      for (int r = 0; r < 2; ++r) {
        const kvz_pixel *in = top_left + src_stride * (yy + r) + 3;
        int16_t *out = dst + (yy + r) * dst_stride;
        int xx = 0;
        for (; xx + 7 < width; xx += 8) {
          _mm_storeu_si128((__m128i *)(out + xx),
                           _mm_slli_epi16(_mm_cvtepu8_epi16(_mm_loadl_epi64((const __m128i *)(in + xx))), 6));
        }
        if (xx < width) {
          uint32_t px;
          memcpy(&px, in + xx, 4);
          _mm_storel_epi64((__m128i *)(out + xx), _mm_slli_epi16(_mm_cvtepu8_epi16(_mm_cvtsi32_si128((int)px)), 6));
        }
      }
    }
    return;
  }

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

// wrap16: store sum >> 6 as int16 before rounding, like the generic block
// filters (int16_t sample = ... >> shift2). The generic sample functions keep
// 32 bits, so they pass 0. Only the half-pel diagonal can exceed int16.
static void kvz_ipol_8tap_ver_im_px_avx2(int8_t *filter,
  int width,
  int height,
  int16_t *src,
  int16_t src_stride,
  kvz_pixel *dst,
  int16_t dst_stride,
  int wrap16)
{
  if (filter == kvz_g_luma_filter[0]) {
    // Filter row 0 is a single tap of 64 at position 3: (64 * h) >> 6 == h,
    // with or without the int16 store, so the result is clip((h + 32) >> 6).
    // Rows are written in pairs like the general loop below.
    const __m128i offset = _mm_set1_epi16(32);
    for (int y = 0; y < height; y += 2) {
      for (int r = 0; r < 2; ++r) {
        const int16_t *in = src + (y + r + 3) * src_stride;
        kvz_pixel *out = dst + (y + r) * dst_stride;
        int x = 0;
        for (; x + 7 < width; x += 8) {
          __m128i v = _mm_srai_epi16(_mm_add_epi16(_mm_loadu_si128((const __m128i *)(in + x)), offset), 6);
          _mm_storel_epi64((__m128i *)(out + x), _mm_packus_epi16(v, v));
        }
        for (; x + 3 < width; x += 4) {
          __m128i v = _mm_srai_epi16(_mm_add_epi16(_mm_loadl_epi64((const __m128i *)(in + x)), offset), 6);
          *(uint32_t *)(out + x) = (uint32_t)_mm_cvtsi128_si32(_mm_packus_epi16(v, v));
        }
      }
    }
    return;
  }

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
      if (wrap16) sum = _mm256_srai_epi32(_mm256_slli_epi32(sum, 16), 16);
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
      sum = _mm256_srai_epi32(sum, shift2);
      // Generic assigns the 32-bit value to an int16_t, which wraps on extreme
      // input. Wrap the same way so the saturating pack below is exact.
      sum = _mm256_srai_epi32(_mm256_slli_epi32(sum, 16), 16);
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

  // A 4-tap filter row is 4 bytes: load only those (the upper half of the
  // widened vector is not used).
  int32_t filter_4b;
  memcpy(&filter_4b, filter, sizeof(filter_4b));
  __m128i weights_8b = _mm_set1_epi32(filter_4b);
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

  // A 4-tap filter row is 4 bytes: load only those (the upper half of the
  // widened vector is not used).
  int32_t filter_4b;
  memcpy(&filter_4b, filter, sizeof(filter_4b));
  __m128i weights_8b = _mm_set1_epi32(filter_4b);
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
  kvz_ipol_8tap_ver_im_px_avx2(fir0, width, height, im, hor_stride, dst, dst_stride, 1);

  // Left
  // Copy from the right filtered block and filter the extra column
  for (y = 0; y < height; ++y) {
    x = 0;
    *(uint64_t*)&out_l[y * dst_stride + x] = *(uint64_t*)&out_r[y * dst_stride + x] << 8;
    for (x = 8; x < width; x += 8) *(int64_t*)&out_l[y * dst_stride + x] = *(int64_t*)&out_r[y * dst_stride + x - 1];
    x = 0;
    int16_t sample = 64 * col_pos2[y + 1 + KVZ_LUMA_FILTER_OFFSET] >> shift2;
    sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
    out_l[y * dst_stride + x] = (kvz_pixel)sample;
  }

  // Top
  im = hor_pos0;
  dst = out_t;
  kvz_ipol_8tap_ver_im_px_avx2(fir2, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(fir2, width, height, im, hor_stride, dst, dst_stride, 1);

  // Top-left
  // Copy from the top-right filtered block and filter the extra column
  for (y = 0; y < height; ++y) {
    x = 0;
    int16_t sample = kvz_eight_tap_filter_hor_16bit_avx2(fir2, &col_pos2[y]) >> shift2;
    sample = kvz_fast_clip_16bit_to_pixel((sample + wp_offset1) >> wp_shift1);
    out_tl[y * dst_stride + x] = (kvz_pixel)sample;

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
  out_bl[y * dst_stride + x] = (kvz_pixel)sample;
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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_l, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_r, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_t, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir_b, width, height, im, hor_stride, dst, dst_stride, 1);

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
  kvz_ipol_8tap_ver_im_px_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride, 0);
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
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  int8_t *hor_fir = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  int8_t *ver_fir = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  // Buffer for intermediate values with 3 extra rows 
  // because the loop writes four rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
  int16_t hor_stride = LCU_WIDTH >> shift_w;

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
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  int8_t* hor_fir = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  int8_t* ver_fir = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  // Buffer for intermediate values with 3 extra rows 
  // because the loop writes four rows each iteration.
  ALIGNED(64) int16_t hor_intermediate[KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
  int16_t hor_stride = LCU_WIDTH >> shift_w;

  kvz_ipol_4tap_hor_px_im_avx2(hor_fir, width, height, src, src_stride, hor_intermediate, hor_stride);
  kvz_ipol_4tap_ver_im_hi_avx2(ver_fir, width, height, hor_intermediate, hor_stride, dst, dst_stride);
}


#if KVZ_BIT_DEPTH > 8

/*
 * 10 and 12-bit interpolation.
 *
 * kvz_pixel is uint16_t and the samples are at most PIXEL_MAX < 2^15, so they
 * can be used directly as int16 operands of _mm256_madd_epi16. The products
 * and sums of both passes need 32-bit lanes: the horizontal sum is up to
 * 88 * PIXEL_MAX before the shift by KVZ_BIT_DEPTH - 8 (the shifted value fits
 * int16), and the vertical sum over int16 intermediates is up to about 2^21.
 *
 * The intermediate buffers of the block filters (hor_intermediate and
 * hor_first_cols) have the same layout and contents as in the generic
 * functions, so these functions can be mixed with the generic ones. Like the
 * 8-bit versions, the block filters require the width and height to be
 * multiples of 8, which search_frac guarantees by rounding the block size up.
 *
 * Only the samples the generic functions read are loaded.
 */

#define HBD_SHIFT1 (KVZ_BIT_DEPTH - 8)
#define HBD_WP_SHIFT (14 - KVZ_BIT_DEPTH)
#define HBD_WP_OFFSET (1 << (HBD_WP_SHIFT - 1))

// Rounding of the vertical pass.
enum {
  // clip(((sum >> 6) + offset) >> shift) in 32 bits, as in the sample functions.
  HBD_VER_PX = 0,
  // The block filters store sum >> 6 in an int16_t before the rounding, which
  // wraps for the extreme half-pel diagonal sums (88 * 88 + 24 * 24 times
  // PIXEL_MAX >> (KVZ_BIT_DEPTH - 2) > 32767). All other filter combinations
  // stay below 30000, where HBD_VER_PX gives the same result.
  HBD_VER_PX_WRAP = 1,
  // (int16_t)(sum >> 6), the _hi sample functions.
  HBD_VER_HI = 2,
};

// Taps k and k + 1 as a pair of int16 in every 32-bit lane.
static INLINE __m256i hbd_tap_pair(const int8_t *filter, int k)
{
  const uint32_t lo = (uint16_t)(int16_t)filter[k];
  const uint32_t hi = (uint16_t)(int16_t)filter[k + 1];
  return _mm256_set1_epi32((int32_t)(lo | (hi << 16)));
}

static INLINE kvz_pixel hbd_clip_px(int32_t value)
{
  return (kvz_pixel)CLIP(0, PIXEL_MAX, value);
}

/*
 * Eight horizontal sums per 128-bit lane: out[j] = sum_k f[k] * p[j + k] with
 * the lane's samples p[0..7] in v0 and p[8..14] in v1. The even outputs come
 * from v0 shifted by 0, 2, 4 and 6 samples, the odd ones from v0 shifted by
 * 1, 3, 5 and 7 samples. Returns outputs 0..3 in *lo and 4..7 in *hi.
 */
static FORCE_INLINE void hbd_hor_sums(__m256i v0, __m256i v1, const __m256i *w, int taps,
                                __m256i *lo, __m256i *hi)
{
  __m256i even = _mm256_madd_epi16(v0, w[0]);
  __m256i odd = _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 2), w[0]);
  even = _mm256_add_epi32(even, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 4), w[1]));
  odd = _mm256_add_epi32(odd, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 6), w[1]));
  if (taps == 8) {
    even = _mm256_add_epi32(even, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 8), w[2]));
    odd = _mm256_add_epi32(odd, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 10), w[2]));
    even = _mm256_add_epi32(even, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 12), w[3]));
    odd = _mm256_add_epi32(odd, _mm256_madd_epi16(_mm256_alignr_epi8(v1, v0, 14), w[3]));
  }
  *lo = _mm256_unpacklo_epi32(even, odd);
  *hi = _mm256_unpackhi_epi32(even, odd);
}

// (sum >> shift1) as int16. The shifted sums always fit int16.
static INLINE __m256i hbd_hor_pack(__m256i lo, __m256i hi)
{
  return _mm256_packs_epi32(_mm256_srai_epi32(lo, HBD_SHIFT1), _mm256_srai_epi32(hi, HBD_SHIFT1));
}

/*
 * v1 for eight outputs: samples p[8 .. taps + 6]. The load ends at the last
 * sample used and is shifted down so that the lane starts at p[8].
 */
static FORCE_INLINE __m128i hbd_load_v1(const kvz_pixel *p, int taps)
{
  const __m128i t = _mm_loadu_si128((const __m128i *)(p + taps - 1));
  return taps == 8 ? _mm_srli_si128(t, 2) : _mm_srli_si128(t, 10);
}

// v0 and v1 for four outputs (samples p[0 .. taps + 2]).
static FORCE_INLINE __m128i hbd_load_4out(const kvz_pixel *p, int taps, __m128i *v1)
{
  if (taps == 8) {
    *v1 = _mm_srli_si128(_mm_loadl_epi64((const __m128i *)(p + 7)), 2);
    return _mm_loadu_si128((const __m128i *)p);
  }
  *v1 = _mm_setzero_si128();
  // p[0..3] and p[3..6]; sample 3 is in both.
  return _mm_or_si128(_mm_loadl_epi64((const __m128i *)p),
                      _mm_slli_si128(_mm_loadl_epi64((const __m128i *)(p + 3)), 6));
}

static INLINE __m256i hbd_combine(__m128i lane0, __m128i lane1)
{
  return _mm256_inserti128_si256(_mm256_castsi128_si256(lane0), lane1, 1);
}

static int16_t hbd_hor_scalar(const int8_t *filter, int taps, const kvz_pixel *p)
{
  int32_t sum = 0;
  for (int k = 0; k < taps; ++k) sum += filter[k] * p[k];
  return (int16_t)(sum >> HBD_SHIFT1);
}

/*
 * Horizontal pass: dst[y][x] = (sum_k f[k] * src[y][x + k]) >> shift1 for
 * `rows` rows. src points to the first tap of output (0, 0).
 */
static FORCE_INLINE void hbd_hor_pass_taps(const int8_t *filter, const int taps, int width, int rows,
                                           const kvz_pixel *src, int src_stride,
                                           int16_t *dst, int dst_stride)
{
  const int center = taps / 2 - 1;

  if (filter[center] == 64) {
    // Integer position: (64 * p) >> shift1 == p << (6 - shift1).
    for (int y = 0; y < rows; ++y) {
      const kvz_pixel *s = src + y * src_stride + center;
      int16_t *d = dst + y * dst_stride;
      int x = 0;
      for (; x + 16 <= width; x += 16) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(s + x));
        _mm256_storeu_si256((__m256i *)(d + x), _mm256_slli_epi16(v, 6 - HBD_SHIFT1));
      }
      if (x + 8 <= width) {
        __m128i v = _mm_loadu_si128((const __m128i *)(s + x));
        _mm_storeu_si128((__m128i *)(d + x), _mm_slli_epi16(v, 6 - HBD_SHIFT1));
        x += 8;
      }
      if (x + 4 <= width) {
        __m128i v = _mm_loadl_epi64((const __m128i *)(s + x));
        _mm_storel_epi64((__m128i *)(d + x), _mm_slli_epi16(v, 6 - HBD_SHIFT1));
        x += 4;
      }
      for (; x < width; ++x) d[x] = (int16_t)(s[x] << (6 - HBD_SHIFT1));
    }
    return;
  }

  __m256i w[4];
  for (int k = 0; k < taps / 2; ++k) w[k] = hbd_tap_pair(filter, 2 * k);

  int x = 0;

  // 16 outputs of one row: lane 0 has outputs 0..7 and lane 1 outputs 8..15.
  for (; x + 16 <= width; x += 16) {
    for (int y = 0; y < rows; ++y) {
      const kvz_pixel *p = src + y * src_stride + x;
      __m256i v0 = _mm256_loadu_si256((const __m256i *)p);
      __m256i v1 = _mm256_loadu_si256((const __m256i *)(p + taps - 1));
      v1 = taps == 8 ? _mm256_srli_si256(v1, 2) : _mm256_srli_si256(v1, 10);
      __m256i lo, hi;
      hbd_hor_sums(v0, v1, w, taps, &lo, &hi);
      _mm256_storeu_si256((__m256i *)(dst + y * dst_stride + x), hbd_hor_pack(lo, hi));
    }
  }

  // Narrower columns: two rows at a time, one in each lane. With an odd row
  // count the last row is computed twice.
  if (x + 8 <= width) {
    for (int y = 0; y < rows; y += 2) {
      const int y1 = y + 1 < rows ? y + 1 : y;
      const kvz_pixel *p0 = src + y * src_stride + x;
      const kvz_pixel *p1 = src + y1 * src_stride + x;
      __m256i v0 = hbd_combine(_mm_loadu_si128((const __m128i *)p0), _mm_loadu_si128((const __m128i *)p1));
      __m256i v1 = hbd_combine(hbd_load_v1(p0, taps), hbd_load_v1(p1, taps));
      __m256i lo, hi;
      hbd_hor_sums(v0, v1, w, taps, &lo, &hi);
      __m256i res = hbd_hor_pack(lo, hi);
      _mm_storeu_si128((__m128i *)(dst + y * dst_stride + x), _mm256_castsi256_si128(res));
      _mm_storeu_si128((__m128i *)(dst + y1 * dst_stride + x), _mm256_extracti128_si256(res, 1));
    }
    x += 8;
  }

  if (x + 4 <= width) {
    for (int y = 0; y < rows; y += 2) {
      const int y1 = y + 1 < rows ? y + 1 : y;
      __m128i v1_0, v1_1;
      __m128i v0_0 = hbd_load_4out(src + y * src_stride + x, taps, &v1_0);
      __m128i v0_1 = hbd_load_4out(src + y1 * src_stride + x, taps, &v1_1);
      __m256i lo, hi;
      hbd_hor_sums(hbd_combine(v0_0, v0_1), hbd_combine(v1_0, v1_1), w, taps, &lo, &hi);
      __m256i res = hbd_hor_pack(lo, hi);
      _mm_storel_epi64((__m128i *)(dst + y * dst_stride + x), _mm256_castsi256_si128(res));
      _mm_storel_epi64((__m128i *)(dst + y1 * dst_stride + x), _mm256_extracti128_si256(res, 1));
    }
    x += 4;
  }

  // Chroma blocks of width 2 and 6.
  for (; x < width; ++x) {
    for (int y = 0; y < rows; ++y) {
      dst[y * dst_stride + x] = hbd_hor_scalar(filter, taps, src + y * src_stride + x);
    }
  }
}

static void hbd_hor_pass_avx2(const int8_t *filter, int taps, int width, int rows,
                              const kvz_pixel *src, int src_stride,
                              int16_t *dst, int dst_stride)
{
  if (taps == 8) hbd_hor_pass_taps(filter, 8, width, rows, src, src_stride, dst, dst_stride);
  else hbd_hor_pass_taps(filter, 4, width, rows, src, src_stride, dst, dst_stride);
}

// Rounds the 32-bit vertical sums of two rows of eight (or four) columns and
// packs them to 16 bits. Lane 0 holds the first row and lane 1 the second.
static FORCE_INLINE __m256i hbd_ver_round(__m256i lo, __m256i hi, const int mode)
{
  const __m256i max = _mm256_set1_epi16(PIXEL_MAX);
  if (mode == HBD_VER_PX) {
    // ((sum >> 6) + offset) >> shift == (sum + (offset << 6)) >> (shift + 6)
    const __m256i rnd = _mm256_set1_epi32(HBD_WP_OFFSET << 6);
    lo = _mm256_srai_epi32(_mm256_add_epi32(lo, rnd), 6 + HBD_WP_SHIFT);
    hi = _mm256_srai_epi32(_mm256_add_epi32(hi, rnd), 6 + HBD_WP_SHIFT);
    return _mm256_min_epu16(_mm256_packus_epi32(lo, hi), max);
  }
  // int16_t sample = sum >> 6
  lo = _mm256_srai_epi32(_mm256_slli_epi32(_mm256_srai_epi32(lo, 6), 16), 16);
  hi = _mm256_srai_epi32(_mm256_slli_epi32(_mm256_srai_epi32(hi, 6), 16), 16);
  if (mode == HBD_VER_HI) return _mm256_packs_epi32(lo, hi);
  const __m256i off = _mm256_set1_epi32(HBD_WP_OFFSET);
  lo = _mm256_srai_epi32(_mm256_add_epi32(lo, off), HBD_WP_SHIFT);
  hi = _mm256_srai_epi32(_mm256_add_epi32(hi, off), HBD_WP_SHIFT);
  return _mm256_min_epu16(_mm256_packus_epi32(lo, hi), max);
}

static void hbd_ver_scalar(const int8_t *filter, int taps, const int16_t *src, int src_stride,
                           kvz_pixel *dst, int16_t *dst_hi, int mode)
{
  int32_t sum = 0;
  for (int k = 0; k < taps; ++k) sum += filter[k] * src[k * src_stride];
  if (mode == HBD_VER_PX) {
    *dst = hbd_clip_px(((sum >> 6) + HBD_WP_OFFSET) >> HBD_WP_SHIFT);
  } else {
    const int16_t sample = (int16_t)(sum >> 6);
    if (mode == HBD_VER_HI) *dst_hi = sample;
    else *dst = hbd_clip_px((sample + HBD_WP_OFFSET) >> HBD_WP_SHIFT);
  }
}

/*
 * One column strip of eight (n = 8) or four (n = 4) columns, two output rows
 * per iteration. Every input row is broadcast to both lanes and
 * R(k) = [row k | row k + 1]. Unpacking R(k) with R(k + 1) gives the sample
 * pairs of rows (k, k + 1) for the first output row in lane 0 and of rows
 * (k + 1, k + 2) for the second output row in lane 1, and each pair is used
 * for taps / 2 consecutive iterations. With an odd height the last row is
 * computed twice; no row after src row height + taps - 2 is read.
 */
static FORCE_INLINE void hbd_ver_strip(const int16_t *s, int src_stride, int height,
                                       const __m256i *w, const int taps, const int n, const int mode,
                                       kvz_pixel *dst, int16_t *dst_hi, int dst_stride)
{
#define HBD_BROW(r) (n == 8 ? \
    _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(s + (r) * src_stride))) : \
    _mm256_broadcastq_epi64(_mm_loadl_epi64((const __m128i *)(s + (r) * src_stride))))

  const int last_row = height + taps - 2;
  __m256i l0, l1 = _mm256_setzero_si256(), l2 = _mm256_setzero_si256();
  __m256i h0, h1 = _mm256_setzero_si256(), h2 = _mm256_setzero_si256();
  __m256i b0 = HBD_BROW(0);
  __m256i b1 = HBD_BROW(1);
  __m256i b2 = HBD_BROW(2);
  __m256i r0 = _mm256_blend_epi32(b0, b1, 0xF0);
  __m256i r1 = _mm256_blend_epi32(b1, b2, 0xF0);
  l0 = _mm256_unpacklo_epi16(r0, r1);
  h0 = _mm256_unpackhi_epi16(r0, r1);
  if (taps == 8) {
    b1 = HBD_BROW(3);
    b0 = HBD_BROW(4);
    r0 = _mm256_blend_epi32(b2, b1, 0xF0);
    r1 = _mm256_blend_epi32(b1, b0, 0xF0);
    l1 = _mm256_unpacklo_epi16(r0, r1);
    h1 = _mm256_unpackhi_epi16(r0, r1);
    b1 = HBD_BROW(5);
    b2 = HBD_BROW(6);
    r0 = _mm256_blend_epi32(b0, b1, 0xF0);
    r1 = _mm256_blend_epi32(b1, b2, 0xF0);
    l2 = _mm256_unpacklo_epi16(r0, r1);
    h2 = _mm256_unpackhi_epi16(r0, r1);
  }
  // b2 is the last row loaded.

  for (int y = 0; y < height; y += 2) {
    b0 = b2;
    b1 = HBD_BROW(y + taps - 1);
    b2 = HBD_BROW(MIN(y + taps, last_row));
    r0 = _mm256_blend_epi32(b0, b1, 0xF0);
    r1 = _mm256_blend_epi32(b1, b2, 0xF0);
    const __m256i ln = _mm256_unpacklo_epi16(r0, r1);
    const __m256i hn = n == 8 ? _mm256_unpackhi_epi16(r0, r1) : ln;

    __m256i lo, hi;
    if (taps == 8) {
      lo = _mm256_add_epi32(_mm256_add_epi32(_mm256_madd_epi16(l0, w[0]), _mm256_madd_epi16(l1, w[1])),
                            _mm256_add_epi32(_mm256_madd_epi16(l2, w[2]), _mm256_madd_epi16(ln, w[3])));
      l0 = l1;
      l1 = l2;
      l2 = ln;
      if (n == 8) {
        hi = _mm256_add_epi32(_mm256_add_epi32(_mm256_madd_epi16(h0, w[0]), _mm256_madd_epi16(h1, w[1])),
                              _mm256_add_epi32(_mm256_madd_epi16(h2, w[2]), _mm256_madd_epi16(hn, w[3])));
        h0 = h1;
        h1 = h2;
        h2 = hn;
      } else {
        hi = lo;
      }
    } else {
      lo = _mm256_add_epi32(_mm256_madd_epi16(l0, w[0]), _mm256_madd_epi16(ln, w[1]));
      l0 = ln;
      if (n == 8) {
        hi = _mm256_add_epi32(_mm256_madd_epi16(h0, w[0]), _mm256_madd_epi16(hn, w[1]));
        h0 = hn;
      } else {
        hi = lo;
      }
    }

    const __m256i res = hbd_ver_round(lo, hi, mode);
    const __m128i res0 = _mm256_castsi256_si128(res);
    const __m128i res1 = _mm256_extracti128_si256(res, 1);
    const int second = y + 1 < height;
    void *d = mode == HBD_VER_HI ? (void *)(dst_hi + y * dst_stride) : (void *)(dst + y * dst_stride);
    void *d1 = mode == HBD_VER_HI ? (void *)(dst_hi + (y + 1) * dst_stride) : (void *)(dst + (y + 1) * dst_stride);
    if (n == 8) {
      _mm_storeu_si128((__m128i *)d, res0);
      if (second) _mm_storeu_si128((__m128i *)d1, res1);
    } else {
      _mm_storel_epi64((__m128i *)d, res0);
      if (second) _mm_storel_epi64((__m128i *)d1, res1);
    }
  }
#undef HBD_BROW
}

static FORCE_INLINE void hbd_ver_pass_taps(const int8_t *filter, const int taps, int width, int height,
                                           const int16_t *src, int src_stride,
                                           kvz_pixel *dst, int16_t *dst_hi, int dst_stride, const int mode)
{
  __m256i w[4];
  for (int k = 0; k < taps / 2; ++k) w[k] = hbd_tap_pair(filter, 2 * k);

  int x = 0;
  for (; x + 8 <= width; x += 8) {
    hbd_ver_strip(src + x, src_stride, height, w, taps, 8, mode,
                  mode == HBD_VER_HI ? NULL : dst + x, mode == HBD_VER_HI ? dst_hi + x : NULL, dst_stride);
  }
  if (x + 4 <= width) {
    hbd_ver_strip(src + x, src_stride, height, w, taps, 4, mode,
                  mode == HBD_VER_HI ? NULL : dst + x, mode == HBD_VER_HI ? dst_hi + x : NULL, dst_stride);
    x += 4;
  }
  // Chroma blocks of width 2 and 6.
  for (; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      hbd_ver_scalar(filter, taps, src + y * src_stride + x, src_stride,
                     mode == HBD_VER_HI ? NULL : dst + y * dst_stride + x,
                     mode == HBD_VER_HI ? dst_hi + y * dst_stride + x : NULL, mode);
    }
  }
}

/*
 * Vertical pass over int16 intermediates. src points to the first tap row of
 * output row 0 and has height + taps - 1 rows. Writes pixels to dst, or int16
 * values to dst_hi for HBD_VER_HI.
 */
static void hbd_ver_pass_avx2(const int8_t *filter, int taps, int width, int height,
                              const int16_t *src, int src_stride,
                              kvz_pixel *dst, int16_t *dst_hi, int dst_stride, int mode)
{
  const int center = taps / 2 - 1;

  if (filter[center] == 64) {
    // Integer position: (64 * h) >> 6 == h, so the int16 store does not wrap,
    // and h + offset fits int16.
    for (int y = 0; y < height; ++y) {
      const int16_t *s = src + (y + center) * src_stride;
      if (mode == HBD_VER_HI) {
        memcpy(dst_hi + y * dst_stride, s, width * sizeof(int16_t));
        continue;
      }
      kvz_pixel *d = dst + y * dst_stride;
      int x = 0;
      for (; x + 16 <= width; x += 16) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(s + x));
        v = _mm256_srai_epi16(_mm256_add_epi16(v, _mm256_set1_epi16(HBD_WP_OFFSET)), HBD_WP_SHIFT);
        v = _mm256_min_epi16(_mm256_max_epi16(v, _mm256_setzero_si256()), _mm256_set1_epi16(PIXEL_MAX));
        _mm256_storeu_si256((__m256i *)(d + x), v);
      }
      for (; x + 8 <= width; x += 8) {
        __m128i v = _mm_loadu_si128((const __m128i *)(s + x));
        v = _mm_srai_epi16(_mm_add_epi16(v, _mm_set1_epi16(HBD_WP_OFFSET)), HBD_WP_SHIFT);
        v = _mm_min_epi16(_mm_max_epi16(v, _mm_setzero_si128()), _mm_set1_epi16(PIXEL_MAX));
        _mm_storeu_si128((__m128i *)(d + x), v);
      }
      for (; x < width; ++x) d[x] = hbd_clip_px((s[x] + HBD_WP_OFFSET) >> HBD_WP_SHIFT);
    }
    return;
  }

  if (taps == 8) {
    switch (mode) {
      case HBD_VER_PX: hbd_ver_pass_taps(filter, 8, width, height, src, src_stride, dst, dst_hi, dst_stride, HBD_VER_PX); break;
      case HBD_VER_PX_WRAP: hbd_ver_pass_taps(filter, 8, width, height, src, src_stride, dst, dst_hi, dst_stride, HBD_VER_PX_WRAP); break;
      default: hbd_ver_pass_taps(filter, 8, width, height, src, src_stride, dst, dst_hi, dst_stride, HBD_VER_HI); break;
    }
  } else {
    if (mode == HBD_VER_PX) hbd_ver_pass_taps(filter, 4, width, height, src, src_stride, dst, dst_hi, dst_stride, HBD_VER_PX);
    else hbd_ver_pass_taps(filter, 4, width, height, src, src_stride, dst, dst_hi, dst_stride, HBD_VER_HI);
  }
}

/*
 * The first column of the block filters: col[y] = (sum_k f[k] * src[y][k]) >> shift1
 * for y < rows, src pointing to the first tap of row 0. Four rows at a time.
 */
static void hbd_hor_first_col_avx2(const int8_t *filter, const kvz_pixel *src, int src_stride,
                                   int rows, int16_t *col)
{
  const __m128i f = _mm_cvtepi8_epi16(_mm_loadl_epi64((const __m128i *)filter));
  int y = 0;
  for (; y + 4 <= rows; y += 4) {
    const kvz_pixel *p = src + y * src_stride;
    __m128i m0 = _mm_madd_epi16(_mm_loadu_si128((const __m128i *)(p + 0 * src_stride)), f);
    __m128i m1 = _mm_madd_epi16(_mm_loadu_si128((const __m128i *)(p + 1 * src_stride)), f);
    __m128i m2 = _mm_madd_epi16(_mm_loadu_si128((const __m128i *)(p + 2 * src_stride)), f);
    __m128i m3 = _mm_madd_epi16(_mm_loadu_si128((const __m128i *)(p + 3 * src_stride)), f);
    __m128i sum = _mm_hadd_epi32(_mm_hadd_epi32(m0, m1), _mm_hadd_epi32(m2, m3));
    sum = _mm_srai_epi32(sum, HBD_SHIFT1);
    _mm_storel_epi64((__m128i *)(col + y), _mm_packs_epi32(sum, sum));
  }
  for (; y < rows; ++y) col[y] = hbd_hor_scalar(filter, 8, src + y * src_stride);
}

// One sample of a first column: clip(((int16_t)(sum >> 6) + offset) >> shift).
static kvz_pixel hbd_col_sample(const int8_t *filter, const int16_t *col)
{
  int32_t sum = 0;
  for (int k = 0; k < 8; ++k) sum += filter[k] * col[k];
  const int16_t sample = (int16_t)(sum >> 6);
  return hbd_clip_px((sample + HBD_WP_OFFSET) >> HBD_WP_SHIFT);
}

/*
 * Vertical 8-tap filter of a first column, written to column 0 of the n rows
 * of dst: dst[y][0] = hbd_col_sample(filter, &col[y]). n is a multiple of 8.
 */
static void hbd_col_filter_avx2(const int8_t *filter, const int16_t *col, int n,
                                kvz_pixel *dst, int dst_stride)
{
  __m256i w[4];
  for (int k = 0; k < 4; ++k) w[k] = hbd_tap_pair(filter, 2 * k);

  for (int i = 0; i < n; i += 16) {
    const int16_t *p = col + i;
    __m256i v0, v1;
    if (i + 16 <= n) {
      v0 = _mm256_loadu_si256((const __m256i *)p);
      v1 = _mm256_srli_si256(_mm256_loadu_si256((const __m256i *)(p + 7)), 2);
    } else {
      // The last eight outputs; lane 1 is not used.
      __m128i a0 = _mm_loadu_si128((const __m128i *)p);
      __m128i a1 = hbd_load_v1((const kvz_pixel *)p, 8);
      v0 = hbd_combine(a0, a0);
      v1 = hbd_combine(a1, a1);
    }
    __m256i lo, hi;
    hbd_hor_sums(v0, v1, w, 8, &lo, &hi);
    // hbd_ver_round packs lane by lane: the outputs are in order.
    ALIGNED(32) kvz_pixel out[16];
    _mm256_store_si256((__m256i *)out, hbd_ver_round(lo, hi, HBD_VER_PX_WRAP));
    const int count = MIN(16, n - i);
    for (int j = 0; j < count; ++j) dst[(i + j) * dst_stride] = out[j];
  }
}

// dst[y][x] = src[y][x - 1] for 1 <= x < width; dst[y][0] is set by the caller.
// width is a multiple of 8. src[0][-1] must be readable.
static void hbd_copy_shifted(const kvz_pixel *src, kvz_pixel *dst, int width, int height, int stride)
{
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; x += 8) {
      _mm_storeu_si128((__m128i *)(dst + y * stride + x), _mm_loadu_si128((const __m128i *)(src + y * stride + x - 1)));
    }
  }
}

// Copies height rows of width (a multiple of 8) samples.
static void hbd_copy_rows(const kvz_pixel *src, kvz_pixel *dst, int width, int height, int stride)
{
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; x += 8) {
      _mm_storeu_si128((__m128i *)(dst + y * stride + x), _mm_loadu_si128((const __m128i *)(src + y * stride + x)));
    }
  }
}

/*
 * A quarter-pel block of the block filters: the vertical filter over the
 * horizontal intermediate rows starting at im. When shift is set, output
 * column x comes from intermediate column x - 1 and column 0 from the first
 * column buffer col. The vertical pass then reads column -1 (the end of the
 * previous row or buffer) for column 0, which is overwritten.
 */
static void hbd_qpel_block_avx2(const int8_t *filter, const int16_t *im, const int16_t *col,
                                int shift, int width, int height, kvz_pixel *dst)
{
  hbd_ver_pass_avx2(filter, 8, width, height, im - shift, LCU_WIDTH, dst, NULL, LCU_WIDTH, HBD_VER_PX);
  if (shift) hbd_col_filter_avx2(filter, col, height, dst, LCU_WIDTH);
}

static void hbd_filter_hpel_blocks_hor_ver_luma_avx2(const encoder_control_t *encoder,
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
  const int8_t *fir0 = kvz_g_luma_filter[0];
  const int8_t *fir2 = kvz_g_luma_filter[2];
  const int stride = LCU_WIDTH;

  // Horizontally filtered samples from the top row are
  // not needed unless samples for diagonal positions are filtered later.
  const int first_y = fme_level > 1 ? 0 : 1;
  const int rows = height + KVZ_EXT_PADDING_LUMA + 1;

  // First taps of hor_intermediate[][0] and of the first columns.
  const kvz_pixel *col_base = src - KVZ_LUMA_FILTER_OFFSET * src_stride - KVZ_LUMA_FILTER_OFFSET;
  const kvz_pixel *im_base = col_base + 1;

  int16_t *col_pos0 = hor_first_cols[0];
  int16_t *col_pos2 = hor_first_cols[2];

  // HORIZONTAL STEP
  // Integer pixels
  hbd_hor_pass_avx2(fir0, 8, width, rows, im_base, src_stride, hor_intermediate[0], stride);
  for (int y = 0; y < rows; ++y) {
    col_pos0[y] = (int16_t)(col_base[y * src_stride + KVZ_LUMA_FILTER_OFFSET] << (6 - HBD_SHIFT1));
  }

  // Half pixels
  hbd_hor_pass_avx2(fir2, 8, width, rows - first_y, im_base + first_y * src_stride, src_stride,
                    hor_intermediate[1] + first_y * stride, stride);
  hbd_hor_first_col_avx2(fir2, col_base + first_y * src_stride, src_stride, rows - first_y, col_pos2 + first_y);

  // VERTICAL STEP
  kvz_pixel *out_l = filtered[0];
  kvz_pixel *out_r = filtered[1];
  kvz_pixel *out_t = filtered[2];
  kvz_pixel *out_b = filtered[3];

  // Right: only the horizontal filter (rows 4 .. height + 3 of the intermediate)
  hbd_ver_pass_avx2(fir0, 8, width, height, hor_intermediate[1] + stride, stride, out_r, NULL, stride, HBD_VER_PX);

  // Left: the right block shifted by one column and the extra column
  hbd_copy_shifted(out_r, out_l, width, height, stride);
  for (int y = 0; y < height; ++y) {
    out_l[y * stride] = hbd_clip_px((col_pos2[y + KVZ_LUMA_FILTER_OFFSET + 1] + HBD_WP_OFFSET) >> HBD_WP_SHIFT);
  }

  // Top: only the vertical filter. The integer intermediate is the sample
  // times 1 << (6 - shift1), so filtering it and rounding by 6 + wp_shift
  // equals the generic (sum >> shift1 + offset) >> wp_shift.
  hbd_ver_pass_avx2(fir2, 8, width, height, hor_intermediate[0], stride, out_t, NULL, stride, HBD_VER_PX);

  // Bottom: the top block shifted by one row and the last row
  hbd_copy_rows(out_t + stride, out_b, width, height - 1, stride);
  hbd_ver_pass_avx2(fir2, 8, width, 1, hor_intermediate[0] + height * stride, stride,
                    out_b + (height - 1) * stride, NULL, stride, HBD_VER_PX);
}

static void hbd_filter_hpel_blocks_diag_luma_avx2(const encoder_control_t *encoder,
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
  const int8_t *fir2 = kvz_g_luma_filter[2];
  const int stride = LCU_WIDTH;

  const int16_t *hor_pos2 = hor_intermediate[1];
  const int16_t *col_pos2 = hor_first_cols[2];

  kvz_pixel *out_tl = filtered[0];
  kvz_pixel *out_tr = filtered[1];
  kvz_pixel *out_bl = filtered[2];
  kvz_pixel *out_br = filtered[3];

  // Top-right
  hbd_ver_pass_avx2(fir2, 8, width, height, hor_pos2, stride, out_tr, NULL, stride, HBD_VER_PX_WRAP);

  // Top-left: the top-right block shifted by one column and the filtered extra column
  hbd_copy_shifted(out_tr, out_tl, width, height, stride);
  hbd_col_filter_avx2(fir2, col_pos2, height, out_tl, stride);

  // Bottom-right: the top-right block shifted by one row and the last row
  hbd_copy_rows(out_tr + stride, out_br, width, height - 1, stride);
  hbd_ver_pass_avx2(fir2, 8, width, 1, hor_pos2 + height * stride, stride,
                    out_br + (height - 1) * stride, NULL, stride, HBD_VER_PX_WRAP);

  // Bottom-left: the top-left block shifted by one row, the last row from
  // the bottom-right block and the last sample from the column buffer
  hbd_copy_rows(out_tl + stride, out_bl, width, height - 1, stride);
  hbd_copy_shifted(out_br + (height - 1) * stride, out_bl + (height - 1) * stride, width, 1, stride);
  out_bl[(height - 1) * stride] = hbd_col_sample(fir2, col_pos2 + height);
}

static void hbd_filter_qpel_blocks_hor_ver_luma_avx2(const encoder_control_t *encoder,
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
  const int8_t *fir0 = kvz_g_luma_filter[0];
  const int8_t *fir1 = kvz_g_luma_filter[1];
  const int8_t *fir2 = kvz_g_luma_filter[2];
  const int8_t *fir3 = kvz_g_luma_filter[3];
  const int stride = LCU_WIDTH;

  // Horizontal positions. Positions 0 and 2 have already been calculated.
  int16_t *hor_pos0 = hor_intermediate[0];
  int16_t *hor_pos2 = hor_intermediate[1];
  int16_t *hor_pos_l = hor_intermediate[3];
  int16_t *hor_pos_r = hor_intermediate[4];
  const int8_t *hor_fir_l = hpel_off_x != 0 ? fir1 : fir3;
  const int8_t *hor_fir_r = hpel_off_x != 0 ? fir3 : fir1;
  int16_t *col_pos_l = hor_first_cols[1];
  int16_t *col_pos_r = hor_first_cols[3];

  const int16_t *hor_hpel_pos = hpel_off_x != 0 ? hor_pos2 : hor_pos0;
  const int16_t *col_pos_hor = hpel_off_x != 0 ? hor_first_cols[2] : hor_first_cols[0];

  // Specify if integer pixels are filtered from left or/and top integer samples
  const int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  const int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  const int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  const int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;
  const int sample_off_y = hpel_off_y < 0 ? 0 : 1;
  const int sample_off_x = hpel_off_x > -1 ? 1 : 0;

  const int rows = height + KVZ_EXT_PADDING_LUMA + 1;
  const kvz_pixel *col_base = src - KVZ_LUMA_FILTER_OFFSET * src_stride - KVZ_LUMA_FILTER_OFFSET;
  const kvz_pixel *im_base = col_base + 1;

  // HORIZONTAL STEP
  hbd_hor_pass_avx2(hor_fir_l, 8, width, rows, im_base, src_stride, hor_pos_l, stride);
  hbd_hor_first_col_avx2(hor_fir_l, col_base, src_stride, rows, col_pos_l);
  hbd_hor_pass_avx2(hor_fir_r, 8, width, rows, im_base, src_stride, hor_pos_r, stride);
  hbd_hor_first_col_avx2(hor_fir_r, col_base, src_stride, rows, col_pos_r);

  // VERTICAL STEP
  const int8_t *ver_fir_l = hpel_off_y != 0 ? fir2 : fir0;
  const int8_t *ver_fir_r = hpel_off_y != 0 ? fir2 : fir0;
  const int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  const int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  // Left QPEL (1/4 or 3/4 x positions)
  hbd_qpel_block_avx2(ver_fir_l, hor_pos_l + sample_off_y * stride, col_pos_l + sample_off_y,
                      !off_x_fir_l, width, height, filtered[0]);
  // Right QPEL (3/4 or 1/4 x positions)
  hbd_qpel_block_avx2(ver_fir_r, hor_pos_r + sample_off_y * stride, col_pos_r + sample_off_y,
                      !off_x_fir_r, width, height, filtered[1]);
  // Top QPEL (1/4 or 3/4 y positions)
  hbd_qpel_block_avx2(ver_fir_t, hor_hpel_pos + off_y_fir_t * stride, col_pos_hor + off_y_fir_t,
                      !sample_off_x, width, height, filtered[2]);
  // Bottom QPEL (3/4 or 1/4 y positions)
  hbd_qpel_block_avx2(ver_fir_b, hor_hpel_pos + off_y_fir_b * stride, col_pos_hor + off_y_fir_b,
                      !sample_off_x, width, height, filtered[3]);
}

static void hbd_filter_qpel_blocks_diag_luma_avx2(const encoder_control_t *encoder,
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
  const int8_t *fir1 = kvz_g_luma_filter[1];
  const int8_t *fir3 = kvz_g_luma_filter[3];
  const int stride = LCU_WIDTH;

  const int16_t *hor_pos_l = hor_intermediate[3];
  const int16_t *hor_pos_r = hor_intermediate[4];
  const int16_t *col_pos_l = hor_first_cols[1];
  const int16_t *col_pos_r = hor_first_cols[3];

  const int8_t *ver_fir_t = hpel_off_y != 0 ? fir1 : fir3;
  const int8_t *ver_fir_b = hpel_off_y != 0 ? fir3 : fir1;

  // Specify if integer pixels are filtered from left or/and top integer samples
  const int off_x_fir_l = hpel_off_x < 1 ? 0 : 1;
  const int off_x_fir_r = hpel_off_x < 0 ? 0 : 1;
  const int off_y_fir_t = hpel_off_y < 1 ? 0 : 1;
  const int off_y_fir_b = hpel_off_y < 0 ? 0 : 1;

  // Top-left QPEL
  hbd_qpel_block_avx2(ver_fir_t, hor_pos_l + off_y_fir_t * stride, col_pos_l + off_y_fir_t,
                      !off_x_fir_l, width, height, filtered[0]);
  // Top-right QPEL
  hbd_qpel_block_avx2(ver_fir_t, hor_pos_r + off_y_fir_t * stride, col_pos_r + off_y_fir_t,
                      !off_x_fir_r, width, height, filtered[1]);
  // Bottom-left QPEL
  hbd_qpel_block_avx2(ver_fir_b, hor_pos_l + off_y_fir_b * stride, col_pos_l + off_y_fir_b,
                      !off_x_fir_l, width, height, filtered[2]);
  // Bottom-right QPEL
  hbd_qpel_block_avx2(ver_fir_b, hor_pos_r + off_y_fir_b * stride, col_pos_r + off_y_fir_b,
                      !off_x_fir_r, width, height, filtered[3]);
}

/*
 * Motion compensation: horizontal pass into an int16 buffer and vertical
 * pass. When the vertical filter is the integer one, only the rows it uses
 * are filtered horizontally.
 */
static void hbd_sample_avx2(const int8_t *hor_fir, const int8_t *ver_fir, int taps,
                            const kvz_pixel *src, int src_stride, int width, int height,
                            kvz_pixel *dst, int16_t *dst_hi, int dst_stride, int mode)
{
  ALIGNED(32) int16_t im[KVZ_EXT_BLOCK_W_LUMA * LCU_WIDTH];
  const int offset = taps / 2 - 1;
  const kvz_pixel *base = src - offset * src_stride - offset;

  if (ver_fir[offset] == 64) {
    hbd_hor_pass_avx2(hor_fir, taps, width, height, base + offset * src_stride, src_stride,
                      im + offset * LCU_WIDTH, LCU_WIDTH);
  } else {
    hbd_hor_pass_avx2(hor_fir, taps, width, height + taps - 1, base, src_stride, im, LCU_WIDTH);
  }
  hbd_ver_pass_avx2(ver_fir, taps, width, height, im, LCU_WIDTH, dst, dst_hi, dst_stride, mode);
}

static void hbd_sample_quarterpel_luma_avx2(const encoder_control_t *const encoder,
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
  hbd_sample_avx2(kvz_g_luma_filter[mv[0] & 3], kvz_g_luma_filter[mv[1] & 3], 8,
                  src, src_stride, width, height, dst, NULL, dst_stride, HBD_VER_PX);
}

static void hbd_sample_quarterpel_luma_hi_avx2(const encoder_control_t *const encoder,
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
  hbd_sample_avx2(kvz_g_luma_filter[mv[0] & 3], kvz_g_luma_filter[mv[1] & 3], 8,
                  src, src_stride, width, height, NULL, dst, dst_stride, HBD_VER_HI);
}

static void hbd_sample_octpel_chroma_avx2(const encoder_control_t *const encoder,
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
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  const int8_t *hor_fir = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  const int8_t *ver_fir = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  hbd_sample_avx2(hor_fir, ver_fir, 4, src, src_stride, width, height, dst, NULL, dst_stride, HBD_VER_PX);
}

static void hbd_sample_octpel_chroma_hi_avx2(const encoder_control_t *const encoder,
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
  const uint8_t shift_w = encoder->cfg.chroma_shift_w;
  const uint8_t shift_h = encoder->cfg.chroma_shift_h;
  const int mv_mask_x = shift_w ? 7 : 3;
  const int mv_mask_y = shift_h ? 7 : 3;
  const int8_t *hor_fir = kvz_g_chroma_filter[(mv[0] & mv_mask_x) << (1 - shift_w)];
  const int8_t *ver_fir = kvz_g_chroma_filter[(mv[1] & mv_mask_y) << (1 - shift_h)];

  hbd_sample_avx2(hor_fir, ver_fir, 4, src, src_stride, width, height, NULL, dst, dst_stride, HBD_VER_HI);
}

#endif // KVZ_BIT_DEPTH > 8

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
#if KVZ_BIT_DEPTH > 8
  if (bitdepth == KVZ_BIT_DEPTH) {
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_hor_ver_luma", "avx2", 40, &hbd_filter_hpel_blocks_hor_ver_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_hpel_blocks_diag_luma", "avx2", 40, &hbd_filter_hpel_blocks_diag_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_hor_ver_luma", "avx2", 40, &hbd_filter_qpel_blocks_hor_ver_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "filter_qpel_blocks_diag_luma", "avx2", 40, &hbd_filter_qpel_blocks_diag_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma", "avx2", 40, &hbd_sample_quarterpel_luma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma", "avx2", 40, &hbd_sample_octpel_chroma_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_quarterpel_luma_hi", "avx2", 40, &hbd_sample_quarterpel_luma_hi_avx2);
    success &= kvz_strategyselector_register(opaque, "sample_octpel_chroma_hi", "avx2", 40, &hbd_sample_octpel_chroma_hi_avx2);
  }
#endif // KVZ_BIT_DEPTH > 8
#endif //COMPILE_INTEL_AVX2 && defined X86_64
  return success;
}
