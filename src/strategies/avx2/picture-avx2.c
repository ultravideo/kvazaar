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
#include "picture-avx2.h"
#include "strategyselector.h"

#if COMPILE_INTEL_AVX2
#  include "image.h"
#  include <immintrin.h>


static unsigned sad_8bit_8x8_avx2(const pixel *buf1, const pixel *buf2)
{
  __m256i sum;
  {
    // Get SADs for 8x8 pixels and add the results hierarchically into sum0.
    const __m256i *const a = (const __m256i *)buf1;
    const __m256i *const b = (const __m256i *)buf2;

    __m256i sum0, sum1;
    sum0 = _mm256_sad_epu8(_mm256_load_si256(a + 0), _mm256_load_si256(b + 0));
    sum1 = _mm256_sad_epu8(_mm256_load_si256(a + 1), _mm256_load_si256(b + 1));
    sum = _mm256_add_epi32(sum0, sum1);
  }

  // Add the high 128 bits to low 128 bits.
  __m128i mm128_result = _mm_add_epi32(_mm256_castsi256_si128(sum), _mm256_extractf128_si256(sum, 1));
  // Add the high 64 bits  to low 64 bits.
  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, mm128_result);
  return result[0] + result[2];
}


static unsigned sad_8bit_16x16_avx2(const pixel *buf1, const pixel *buf2)
{
  __m256i sum;
  {
    // Get SADs for 16x16 pixels and add the results hierarchically into sum.
    const __m256i *const a = (const __m256i *)buf1;
    const __m256i *const b = (const __m256i *)buf2;

    __m256i sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
    sum0 = _mm256_sad_epu8(_mm256_load_si256(a + 0), _mm256_load_si256(b + 0));
    sum1 = _mm256_sad_epu8(_mm256_load_si256(a + 1), _mm256_load_si256(b + 1));
    sum2 = _mm256_sad_epu8(_mm256_load_si256(a + 2), _mm256_load_si256(b + 2));
    sum3 = _mm256_sad_epu8(_mm256_load_si256(a + 3), _mm256_load_si256(b + 3));
    sum4 = _mm256_sad_epu8(_mm256_load_si256(a + 4), _mm256_load_si256(b + 4));
    sum5 = _mm256_sad_epu8(_mm256_load_si256(a + 5), _mm256_load_si256(b + 5));
    sum6 = _mm256_sad_epu8(_mm256_load_si256(a + 6), _mm256_load_si256(b + 6));
    sum7 = _mm256_sad_epu8(_mm256_load_si256(a + 7), _mm256_load_si256(b + 7));

    sum0 = _mm256_add_epi32(sum0, sum1);
    sum2 = _mm256_add_epi32(sum2, sum3);
    sum4 = _mm256_add_epi32(sum4, sum5);
    sum6 = _mm256_add_epi32(sum6, sum7);

    sum0 = _mm256_add_epi32(sum0, sum2);
    sum4 = _mm256_add_epi32(sum4, sum6);

    sum = _mm256_add_epi32(sum0, sum4);
  }

  // Add the high 128 bits to low 128 bits.
  __m128i mm128_result = _mm_add_epi32(_mm256_castsi256_si128(sum), _mm256_extractf128_si256(sum, 1));
  // Add the high 64 bits  to low 64 bits.
  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, mm128_result);
  return result[0] + result[2];
}


static unsigned sad_8bit_32x32_avx2(const pixel *buf1, const pixel *buf2)
{
  // Do 32x32 in 4 blocks.
  __m256i sum = _mm256_setzero_si256();
  for (int i = 0; i < 32; i += 8) {
    // Get SADs for 32x8 pixels and add the results hierarchically into sum.
    const __m256i *const a = (const __m256i *)buf1 + i;
    const __m256i *const b = (const __m256i *)buf2 + i;

    __m256i sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
    sum0 = _mm256_sad_epu8(_mm256_load_si256(a + 0), _mm256_load_si256(b + 0));
    sum1 = _mm256_sad_epu8(_mm256_load_si256(a + 1), _mm256_load_si256(b + 1));
    sum2 = _mm256_sad_epu8(_mm256_load_si256(a + 2), _mm256_load_si256(b + 2));
    sum3 = _mm256_sad_epu8(_mm256_load_si256(a + 3), _mm256_load_si256(b + 3));
    sum4 = _mm256_sad_epu8(_mm256_load_si256(a + 4), _mm256_load_si256(b + 4));
    sum5 = _mm256_sad_epu8(_mm256_load_si256(a + 5), _mm256_load_si256(b + 5));
    sum6 = _mm256_sad_epu8(_mm256_load_si256(a + 6), _mm256_load_si256(b + 6));
    sum7 = _mm256_sad_epu8(_mm256_load_si256(a + 7), _mm256_load_si256(b + 7));

    sum0 = _mm256_add_epi32(sum0, sum1);
    sum2 = _mm256_add_epi32(sum2, sum3);
    sum4 = _mm256_add_epi32(sum4, sum5);
    sum6 = _mm256_add_epi32(sum6, sum7);

    sum0 = _mm256_add_epi32(sum0, sum2);
    sum4 = _mm256_add_epi32(sum4, sum6);

    sum = _mm256_add_epi32(sum, sum0);
    sum = _mm256_add_epi32(sum, sum4);
  }

  // Add the high 128 bits to low 128 bits.
  __m128i mm128_result = _mm_add_epi32(_mm256_castsi256_si128(sum), _mm256_extractf128_si256(sum, 1));
  // Add the high 64 bits  to low 64 bits.
  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, mm128_result);
  return result[0] + result[2];
}

#endif //COMPILE_INTEL_AVX2


int strategy_register_picture_avx2(void* opaque) {
  bool success = true;
#if COMPILE_INTEL_AVX2
  success &= strategyselector_register(opaque, "sad_8bit_8x8", "avx2", 40, &sad_8bit_8x8_avx2);
  success &= strategyselector_register(opaque, "sad_8bit_16x16", "avx2", 40, &sad_8bit_16x16_avx2);
  success &= strategyselector_register(opaque, "sad_8bit_32x32", "avx2", 40, &sad_8bit_32x32_avx2);
#endif
  return success;
}
