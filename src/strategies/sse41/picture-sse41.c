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

#include "strategies/sse41/picture-sse41.h"

#if COMPILE_INTEL_SSE41
#include <immintrin.h>
#include <stdlib.h>

#include "kvazaar.h"
#include "strategyselector.h"


uint32_t kvz_reg_sad_sse41(const kvz_pixel * const data1, const kvz_pixel * const data2,
                           const int32_t width, const int32_t height, const uint32_t stride1,
                           const uint32_t stride2)
{
  int32_t y, x;
  __m128i sse_inc = _mm_setzero_si128();
  
  // Bytes in block in 128-bit blocks per each scanline, and remainder
  const int32_t largeblock_bytes = width & ~15;
  const int32_t residual_bytes   = width &  15;

  const __m128i rds    = _mm_set1_epi8 (residual_bytes);
  const __m128i ns     = _mm_setr_epi8 (0,  1,  2,  3,  4,  5,  6,  7,
                                        8,  9,  10, 11, 12, 13, 14, 15);
  const __m128i rdmask = _mm_cmpgt_epi8(rds, ns);

  for (y = 0; y < height; ++y) {
    for (x = 0; x < largeblock_bytes; x += 16) {
      __m128i a = _mm_loadu_si128((__m128i const*) &data1[y * stride1 + x]);
      __m128i b = _mm_loadu_si128((__m128i const*) &data2[y * stride2 + x]);
      __m128i curr_sads = _mm_sad_epu8(a, b);
      sse_inc = _mm_add_epi32(sse_inc, curr_sads);
    }
    
    if (residual_bytes) {
      __m128i a = _mm_loadu_si128((__m128i const*) &data1[y * stride1 + x]);
      __m128i b = _mm_loadu_si128((__m128i const*) &data2[y * stride2 + x]);

      __m128i b_masked  = _mm_blendv_epi8(a, b, rdmask);
      __m128i curr_sads = _mm_sad_epu8(a, b_masked);
      sse_inc = _mm_add_epi32(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

#endif //COMPILE_INTEL_SSE41


int kvz_strategy_register_picture_sse41(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_INTEL_SSE41
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "reg_sad", "sse41", 20, &kvz_reg_sad_sse41);
  }
#endif
  return success;
}
