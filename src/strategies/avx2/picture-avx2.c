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


static unsigned sad_8bit_8x8_avx2(const kvz_pixel *buf1, const kvz_pixel *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_8x8_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_16x16_avx2(const kvz_pixel *buf1, const kvz_pixel *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_16x16_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


static unsigned sad_8bit_32x32_avx2(const kvz_pixel *buf1, const kvz_pixel *buf2)
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


static unsigned sad_8bit_64x64_avx2(const kvz_pixel * buf1, const kvz_pixel * buf2)
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

static void hor_add_sub_avx2(__m128i *row0, __m128i *row1){

  __m128i a = _mm_hadd_epi16(*row0, *row1);
  __m128i b = _mm_hsub_epi16(*row0, *row1);

  __m128i c = _mm_hadd_epi16(a, b);
  __m128i d = _mm_hsub_epi16(a, b);

  *row0 = _mm_hadd_epi16(c, d);
  *row1 = _mm_hsub_epi16(c, d);
}

static INLINE void ver_add_sub_avx2(__m128i temp_hor[8], __m128i temp_ver[8]){

  // First stage
  for (int i = 0; i < 8; i += 2){
    temp_ver[i+0] = _mm_hadd_epi16(temp_hor[i + 0], temp_hor[i + 1]);
    temp_ver[i+1] = _mm_hsub_epi16(temp_hor[i + 0], temp_hor[i + 1]);
  }

  // Second stage
  for (int i = 0; i < 8; i += 4){
    temp_hor[i + 0] = _mm_add_epi16(temp_ver[i + 0], temp_ver[i + 2]);
    temp_hor[i + 1] = _mm_add_epi16(temp_ver[i + 1], temp_ver[i + 3]);
    temp_hor[i + 2] = _mm_sub_epi16(temp_ver[i + 0], temp_ver[i + 2]);
    temp_hor[i + 3] = _mm_sub_epi16(temp_ver[i + 1], temp_ver[i + 3]);
  }

  // Third stage
  for (int i = 0; i < 4; ++i){
    temp_ver[i + 0] = _mm_add_epi16(temp_hor[0 + i], temp_hor[4 + i]);
    temp_ver[i + 4] = _mm_sub_epi16(temp_hor[0 + i], temp_hor[4 + i]);
  }
}

static unsigned kvz_satd_8bit_8x8_general_avx2(const kvz_pixel * buf1, unsigned stride1, const kvz_pixel * buf2, unsigned stride2)
{
  __m128i temp_hor[8];
  __m128i temp_ver[8];

  for (int row = 0; row < 8; row += 2){
    for (int i = 0; i < 2; ++i){
      __m128i buf1_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(&buf1[(row + i) * stride1])));
      __m128i buf2_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)(&buf2[(row + i) * stride2])));
      temp_hor[row + i] = _mm_sub_epi16(buf1_row, buf2_row);
    }
    hor_add_sub_avx2(&temp_hor[row], &temp_hor[row + 1]);
  }

  ver_add_sub_avx2(temp_hor, temp_ver);
  
  __m128i sad = _mm_setzero_si128();
  for (int row = 0; row < 8; ++row){
    __m128i abs_value = _mm_abs_epi16(temp_ver[row]);
    sad = _mm_add_epi32(sad, _mm_madd_epi16(abs_value, _mm_set1_epi16(1)));
  }

  sad = _mm_hadd_epi32(sad, sad);
  sad = _mm_hadd_epi32(sad, sad);

  unsigned result = (_mm_cvtsi128_si32(sad) + 2) >> 2;
  return result;
}

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NXN_AVX2(n) \
static unsigned satd_8bit_ ## n ## x ## n ## _avx2( \
  const kvz_pixel * const block1, const kvz_pixel * const block2) \
{ \
  unsigned x, y; \
  unsigned sum = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  sum += kvz_satd_8bit_8x8_general_avx2(&block1[row + x], (n), &block2[row + x], (n)); \
    } \
    } \
  return sum>>(KVZ_BIT_DEPTH-8); \
}

static unsigned satd_8bit_8x8_avx2(
  const kvz_pixel * const block1, const kvz_pixel * const block2) 
{ 
  unsigned x, y; 
  unsigned sum = 0; 
  for (y = 0; y < (8); y += 8) { 
  unsigned row = y * (8); 
  for (x = 0; x < (8); x += 8) { 
  sum += kvz_satd_8bit_8x8_general_avx2(&block1[row + x], (8), &block2[row + x], (8)); 
      } 
      } 
  return sum>>(KVZ_BIT_DEPTH-8); \
}

//SATD_NXN_AVX2(8)
SATD_NXN_AVX2(16)
SATD_NXN_AVX2(32)
SATD_NXN_AVX2(64)

#endif //COMPILE_INTEL_AVX2


int kvz_strategy_register_picture_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  // We don't actually use SAD for intra right now, other than 4x4 for
  // transform skip, but we might again one day and this is some of the
  // simplest code to look at for anyone interested in doing more
  // optimizations, so it's worth it to keep this maintained.
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "avx2", 40, &sad_8bit_8x8_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "avx2", 40, &sad_8bit_16x16_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "avx2", 40, &sad_8bit_32x32_avx2);
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "avx2", 40, &sad_8bit_64x64_avx2);

    success &= kvz_strategyselector_register(opaque, "satd_8x8", "avx2", 40, &satd_8bit_8x8_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "avx2", 40, &satd_8bit_16x16_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "avx2", 40, &satd_8bit_32x32_avx2);
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "avx2", 40, &satd_8bit_64x64_avx2);
  }
#endif
  return success;
}
