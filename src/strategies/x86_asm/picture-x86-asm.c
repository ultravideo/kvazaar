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
#include <stdlib.h>
#include "strategyselector.h"

#if !defined(KVAZAAR_DISABLE_YASM)

#include "picture-x86-asm-sad.h"
#include "picture-x86-asm-satd.h"

#ifdef __GNUC__
__attribute__((__target__("avx")))
#endif

static unsigned kvz_sad_32x32_avx(const pixel *data1, const pixel *data2)
{
  unsigned sad = 0;
  sad += kvz_sad_16x16_avx(data1, data2);
  sad += kvz_sad_16x16_avx(data1 + 8 * 32, data2 + 8 * 32);
  sad += kvz_sad_16x16_avx(data1 + 16 * 32, data2 + 16 * 32);
  sad += kvz_sad_16x16_avx(data1 + 24 * 32, data2 + 24 * 32);
  return sad;
}

static unsigned kvz_sad_32x32_stride_avx(const pixel *data1, const pixel *data2, unsigned stride)
{
  unsigned sad = 0;
  sad += kvz_sad_16x16_stride_avx(data1, data2, stride);
  sad += kvz_sad_16x16_stride_avx(data1 + 16, data2 + 16, stride);
  sad += kvz_sad_16x16_stride_avx(data1 + 16 * stride, data2 + 16 * stride, stride);
  sad += kvz_sad_16x16_stride_avx(data1 + 16 * stride + 16, data2 + 16 * stride + 16, stride);
  return sad;
}

static unsigned kvz_sad_64x64_avx(const pixel *data1, const pixel *data2)
{
  unsigned sad = 0;
  sad += kvz_sad_32x32_avx(data1, data2);
  sad += kvz_sad_32x32_avx(data1 + 16 * 64, data2 + 16 * 64);
  sad += kvz_sad_32x32_avx(data1 + 32 * 64, data2 + 32 * 64);
  sad += kvz_sad_32x32_avx(data1 + 48 * 64, data2 + 48 * 64);
  return sad;
}

static unsigned kvz_sad_64x64_stride_avx(const pixel *data1, const pixel *data2, unsigned stride)
{
  unsigned sad = 0;
  sad += kvz_sad_32x32_stride_avx(data1, data2, stride);
  sad += kvz_sad_32x32_stride_avx(data1 + 32, data2 + 32, stride);
  sad += kvz_sad_32x32_stride_avx(data1 + 32 * stride, data2 + 32 * stride, stride);
  sad += kvz_sad_32x32_stride_avx(data1 + 32 * stride + 32, data2 + 32 * stride + 32, stride);
  return sad;
}

static unsigned kvz_sad_other_avx(const pixel * const data1, const pixel * const data2,
  const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  int y, x;
  unsigned sad = 0;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }

  return sad;
}

static unsigned reg_sad_x86_asm(const pixel * const data1, const pixel * const data2,
const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  if (width == 4 && height == 4) {
    return kvz_sad_4x4_stride_avx(data1, data2, stride1);
  } else if (width == 8 && height == 8) {
    return kvz_sad_8x8_stride_avx(data1, data2, stride1);
  } else if (width == 16 && height == 16) {
    return kvz_sad_16x16_stride_avx(data1, data2, stride1);
  } else if (width == 32 && height == 32) {
    return kvz_sad_32x32_stride_avx(data1, data2, stride1);
  } else if (width == 64 && height == 64) {
    return kvz_sad_64x64_stride_avx(data1, data2, stride1);
  } else {
    return kvz_sad_other_avx(data1, data2, width, height, stride1, stride2);
  }
}

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define KVZ_SATD_NXN(n, pixel_type, suffix) \
  static unsigned kvz_satd_ ## suffix ## _ ## n ## x ## n ## _stride( \
  const pixel_type * const block1, const pixel_type * const block2) \
{ \
  unsigned x, y; \
  unsigned sum = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  sum += kvz_satd_8x8_stride_avx(&block1[row + x], (n), &block2[row + x], (n)); \
  } \
  } \
  return sum; \
}

// Declare these functions to make sure the signature of the macro matches.
cost_pixel_nxn_func kvz_satd_8bit_4x4_avx;
cost_pixel_nxn_func kvz_satd_8bit_8x8_avx;
cost_pixel_nxn_func kvz_satd_8bit_16x16_avx;
cost_pixel_nxn_func kvz_satd_8bit_32x32_avx;
cost_pixel_nxn_func kvz_satd_8bit_64x64_avx;

#endif //COMPILE_INTEL_AVX && !defined(KVAZAAR_DISABLE_YASM)

int strategy_register_picture_x86_asm_avx(void* opaque) {
  bool success = true;
#if !defined(KVAZAAR_DISABLE_YASM)
  success &= strategyselector_register(opaque, "reg_sad", "x86_asm_avx", 30, &reg_sad_x86_asm);

  success &= strategyselector_register(opaque, "sad_8bit_4x4", "x86_asm_avx", 30, &kvz_sad_4x4_avx);
  success &= strategyselector_register(opaque, "sad_8bit_8x8", "x86_asm_avx", 30, &kvz_sad_8x8_avx);
  success &= strategyselector_register(opaque, "sad_8bit_16x16", "x86_asm_avx", 30, &kvz_sad_16x16_avx);
  success &= strategyselector_register(opaque, "sad_8bit_32x32", "x86_asm_avx", 30, &kvz_sad_32x32_avx);
  success &= strategyselector_register(opaque, "sad_8bit_64x64", "x86_asm_avx", 30, &kvz_sad_64x64_avx);

  success &= strategyselector_register(opaque, "satd_8bit_4x4", "x86_asm_avx", 30, &kvz_satd_4x4_avx);
  success &= strategyselector_register(opaque, "satd_8bit_8x8", "x86_asm_avx", 30, &kvz_satd_8x8_avx);
  success &= strategyselector_register(opaque, "satd_8bit_16x16", "x86_asm_avx", 30, &kvz_satd_16x16_avx);
  success &= strategyselector_register(opaque, "satd_8bit_32x32", "x86_asm_avx", 30, &kvz_satd_32x32_avx);
  success &= strategyselector_register(opaque, "satd_8bit_64x64", "x86_asm_avx", 30, &kvz_satd_64x64_avx);
#endif //COMPILE_INTEL_AVX && !defined(KVAZAAR_DISABLE_YASM)
  return success;
}
