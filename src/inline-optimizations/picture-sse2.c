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

#include "../picture.h"
#include <immintrin.h>

#ifdef __SSE2__
static unsigned reg_sad(const pixel * const data1, const pixel * const data2,
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
    
#ifdef __SSE4_1__
    {
      const __m128i a = _mm_loadu_si128((__m128i const*) &data1[y * stride1 + x]);
      const __m128i b = _mm_loadu_si128((__m128i const*) &data2[y * stride2 + x]);
      switch (((width - (width%2)) - x)/2) {
        case 0:
          break;
        case 1:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x01)));
          break;
        case 2:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x03)));
          break;
        case 3:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x07)));
          break;
        case 4:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x0f)));
          break;
        case 5:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x1f)));
          break;
        case 6:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x3f)));
          break;
        case 7:
          sse_inc = _mm_add_epi32(sse_inc, _mm_sad_epu8(a, _mm_blend_epi16(a, b, 0x7f)));
          break;
        default:
          //Should not happen
          assert(0);
      }
      x = (width - (width%2));
    }
#endif //__SSE4_1__

    for (; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }
  _mm_storeu_si128((__m128i*) sse_inc_array, sse_inc);
  sad += sse_inc_array[0] + sse_inc_array[1];

  return sad;
}
#else
#error picture-sse2.c requires __SSE2__
#endif //__SSE2__
