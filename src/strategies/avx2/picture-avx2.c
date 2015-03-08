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
#include "picture-avx2.h"
#include "strategyselector.h"

#if COMPILE_INTEL_AVX2
#  include "image.h"
#  include <immintrin.h>


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
* \brief Calculate SAD for 32x32 bytes in continuous memory.
*/
static INLINE __m256i inline_8bit_sad_32x32_avx2(const __m256i *const a, const __m256i *const b)
{
  const unsigned size_of_16x16 = 16 * 16 / sizeof(__m256i);

  // Calculate in 4 chunks of 32x8.
  __m256i sum0, sum1, sum2, sum3;
  sum0 = inline_8bit_sad_16x16_avx2(a + 0 * size_of_16x16, b + 0 * size_of_16x16);
  sum1 = inline_8bit_sad_16x16_avx2(a + 1 * size_of_16x16, b + 1 * size_of_16x16);
  sum2 = inline_8bit_sad_16x16_avx2(a + 2 * size_of_16x16, b + 2 * size_of_16x16);
  sum3 = inline_8bit_sad_16x16_avx2(a + 3 * size_of_16x16, b + 3 * size_of_16x16);

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


static unsigned sad_8bit_8x8_avx2(const pixel_t *buf1, const pixel_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_8x8_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_16x16_avx2(const pixel_t *buf1, const pixel_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_16x16_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_32x32_avx2(const pixel_t *buf1, const pixel_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;

  __m256i sum = inline_8bit_sad_32x32_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


#endif //COMPILE_INTEL_AVX2


int strategy_register_picture_avx2(void* opaque)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  success &= strategyselector_register(opaque, "sad_8bit_8x8", "avx2", 40, &sad_8bit_8x8_avx2);
  success &= strategyselector_register(opaque, "sad_8bit_16x16", "avx2", 40, &sad_8bit_16x16_avx2);
  success &= strategyselector_register(opaque, "sad_8bit_32x32", "avx2", 40, &sad_8bit_32x32_avx2);
#endif
  return success;
}
