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

#include "strategies/sse2/picture-sse2.h"

#if COMPILE_INTEL_SSE2
#include "kvazaar.h"
#if KVZ_BIT_DEPTH == 8
#include <immintrin.h>
#include <stdlib.h>

#include "strategyselector.h"


static unsigned reg_sad_sse2(const uint8_t * const data1, const uint8_t * const data2,
                        const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  int y, x;
  unsigned sad = 0;
  __m128i sse_inc = _mm_setzero_si128 ();
  long long int sse_inc_array[2];
  
  for (y = 0; y < height; ++y) {
    for (x = 0; x <= width-16; x+=16) {
      const __m128i a = _mm_loadu_si128((__m128i const*) &data1[y * stride1 + x]);
      const __m128i b = _mm_loadu_si128((__m128i const*) &data2[y * stride2 + x]);
      sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a,b));
    }

    for (; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }
  _mm_storeu_si128((__m128i*) sse_inc_array, sse_inc);
  sad += sse_inc_array[0] + sse_inc_array[1];

  return sad;
}

static unsigned sad_8bit_4x4_sse2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m128i *const mbuf1 = (const __m128i *)buf1;
  const __m128i *const mbuf2 = (const __m128i *)buf2;

  __m128i sum = _mm_sad_epu8(_mm_load_si128(mbuf1), _mm_load_si128(mbuf2));

  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, sum);
  return result[0] + result[2];
}

#endif // KVZ_BIT_DEPTH == 8
#endif //COMPILE_INTEL_SSE2

int kvz_strategy_register_picture_sse2(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_INTEL_SSE2
#if KVZ_BIT_DEPTH == 8
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "reg_sad", "sse2", 10, &reg_sad_sse2);
    success &= kvz_strategyselector_register(opaque, "sad_4x4", "sse2", 10, &sad_8bit_4x4_sse2);
  }
#endif // KVZ_BIT_DEPTH == 8
#endif
  return success;
}


