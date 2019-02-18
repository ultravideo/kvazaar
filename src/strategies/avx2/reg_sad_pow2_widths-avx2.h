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

#ifndef REG_SAD_POW2_WIDTHS_AVX2_H_
#define REG_SAD_POW2_WIDTHS_AVX2_H_

#include "strategies/sse41/reg_sad_pow2_widths-sse41.h"
#include "kvazaar.h"

static INLINE uint32_t reg_sad_w32(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1));
    __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2));
    __m256i c = _mm256_loadu_si256((const __m256i *)(data1 + (y + 1) * stride1));
    __m256i d = _mm256_loadu_si256((const __m256i *)(data2 + (y + 1) * stride2));
    __m256i e = _mm256_loadu_si256((const __m256i *)(data1 + (y + 2) * stride1));
    __m256i f = _mm256_loadu_si256((const __m256i *)(data2 + (y + 2) * stride2));
    __m256i g = _mm256_loadu_si256((const __m256i *)(data1 + (y + 3) * stride1));
    __m256i h = _mm256_loadu_si256((const __m256i *)(data2 + (y + 3) * stride2));

    __m256i curr_sads_ab = _mm256_sad_epu8(a, b);
    __m256i curr_sads_cd = _mm256_sad_epu8(c, d);
    __m256i curr_sads_ef = _mm256_sad_epu8(e, f);
    __m256i curr_sads_gh = _mm256_sad_epu8(g, h);

    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ab);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_cd);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ef);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1));
      __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2));

      __m256i curr_sads = _mm256_sad_epu8(a, b);
      avx_inc = _mm256_add_epi64(avx_inc, curr_sads);
    }
  }

  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

static INLINE uint32_t reg_sad_w64(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;

  const int32_t height_twoline_groups = height & ~1;
  const int32_t height_residual_lines = height &  1;

  for (y = 0; y < height_twoline_groups; y += 2) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1));
    __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2));
    __m256i c = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1 + 32));
    __m256i d = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2 + 32));

    __m256i e = _mm256_loadu_si256((const __m256i *)(data1 + (y + 1) * stride1));
    __m256i f = _mm256_loadu_si256((const __m256i *)(data2 + (y + 1) * stride2));
    __m256i g = _mm256_loadu_si256((const __m256i *)(data1 + (y + 1) * stride1 + 32));
    __m256i h = _mm256_loadu_si256((const __m256i *)(data2 + (y + 1) * stride2 + 32));

    __m256i curr_sads_ab = _mm256_sad_epu8(a, b);
    __m256i curr_sads_cd = _mm256_sad_epu8(c, d);
    __m256i curr_sads_ef = _mm256_sad_epu8(e, f);
    __m256i curr_sads_gh = _mm256_sad_epu8(g, h);

    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ab);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_cd);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ef);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1));
      __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2));
      __m256i c = _mm256_loadu_si256((const __m256i *)(data1 + (y + 0) * stride1 + 32));
      __m256i d = _mm256_loadu_si256((const __m256i *)(data2 + (y + 0) * stride2 + 32));

      __m256i curr_sads_ab = _mm256_sad_epu8(a, b);
      __m256i curr_sads_cd = _mm256_sad_epu8(c, d);
      avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ab);
      avx_inc = _mm256_add_epi64(avx_inc, curr_sads_cd);
    }
  }

  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t hor_sad_avx2_w32(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                 int32_t height, uint32_t pic_stride, uint32_t ref_stride,
                                 uint32_t left, uint32_t right)
{
  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  __m256i ns = _mm256_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                                8,  9,  10, 11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20, 21, 22, 23,
                                24, 25, 26, 27, 28, 29, 30, 31);
  __m256i epol_mask;
  int32_t border_pix_off;

  if (left) {
    border_pix_off          = left;
    __m256i first_valid_idx = _mm256_set1_epi8(left);

    epol_mask = _mm256_cmpgt_epi8(first_valid_idx, ns);
  } else {
    border_pix_off          = 31 - right;
    __m256i last_valid_idx  = _mm256_set1_epi8(border_pix_off);

    epol_mask = _mm256_cmpgt_epi8(ns, last_valid_idx);
  }

  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;
  for (y = 0; y < height_fourline_groups; y += 4) {
    __m256i a = _mm256_loadu_si256((__m256i *)(pic_data + (y + 0) * pic_stride));
    __m256i b = _mm256_loadu_si256((__m256i *)(ref_data + (y + 0) * ref_stride));
    __m256i c = _mm256_loadu_si256((__m256i *)(pic_data + (y + 1) * pic_stride));
    __m256i d = _mm256_loadu_si256((__m256i *)(ref_data + (y + 1) * ref_stride));
    __m256i e = _mm256_loadu_si256((__m256i *)(pic_data + (y + 2) * pic_stride));
    __m256i f = _mm256_loadu_si256((__m256i *)(ref_data + (y + 2) * ref_stride));
    __m256i g = _mm256_loadu_si256((__m256i *)(pic_data + (y + 3) * pic_stride));
    __m256i h = _mm256_loadu_si256((__m256i *)(ref_data + (y + 3) * ref_stride));

    __m256i border_px_b  = _mm256_set1_epi8   (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_pix_off));
    __m256i border_px_d  = _mm256_set1_epi8   (*(uint8_t *)(ref_data + (y + 1) * ref_stride + border_pix_off));
    __m256i border_px_f  = _mm256_set1_epi8   (*(uint8_t *)(ref_data + (y + 2) * ref_stride + border_pix_off));
    __m256i border_px_h  = _mm256_set1_epi8   (*(uint8_t *)(ref_data + (y + 3) * ref_stride + border_pix_off));

    __m256i b_epol       = _mm256_blendv_epi8(b, border_px_b, epol_mask);
    __m256i d_epol       = _mm256_blendv_epi8(d, border_px_d, epol_mask);
    __m256i f_epol       = _mm256_blendv_epi8(f, border_px_f, epol_mask);
    __m256i h_epol       = _mm256_blendv_epi8(h, border_px_h, epol_mask);

    __m256i curr_sads_ab = _mm256_sad_epu8(a, b_epol);
    __m256i curr_sads_cd = _mm256_sad_epu8(c, d_epol);
    __m256i curr_sads_ef = _mm256_sad_epu8(e, f_epol);
    __m256i curr_sads_gh = _mm256_sad_epu8(g, h_epol);

    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ab);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_cd);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ef);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m256i a = _mm256_loadu_si256((__m256i *)(pic_data + y * pic_stride));
      __m256i b = _mm256_loadu_si256((__m256i *)(ref_data + y * ref_stride));

      __m256i border_px_b = _mm256_set1_epi8  (*(uint8_t *)(ref_data + y * ref_stride + border_pix_off));
      __m256i b_epol      = _mm256_blendv_epi8(b, border_px_b, epol_mask);

      __m256i curr_sads_ab = _mm256_sad_epu8(a, b_epol);

      avx_inc = _mm256_add_epi64(avx_inc, curr_sads_ab);
    }
  }

  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

#endif
