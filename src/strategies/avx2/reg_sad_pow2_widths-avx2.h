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

#endif
