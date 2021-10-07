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

#include "strategies/x86_asm/picture-x86-asm.h"

#if defined(KVZ_COMPILE_ASM)
#include "kvazaar.h"
#if KVZ_BIT_DEPTH == 8
#include <stdlib.h>

#include "strategies/x86_asm/picture-x86-asm-sad.h"
#include "strategies/x86_asm/picture-x86-asm-satd.h"
#include "strategies/sse41/picture-sse41.h"
#include "strategyselector.h"


static unsigned kvz_sad_32x32_avx(const uint8_t *data1, const uint8_t *data2)
{
  unsigned sad = 0;
  sad += kvz_sad_16x16_avx(data1, data2);
  sad += kvz_sad_16x16_avx(data1 + 8 * 32, data2 + 8 * 32);
  sad += kvz_sad_16x16_avx(data1 + 16 * 32, data2 + 16 * 32);
  sad += kvz_sad_16x16_avx(data1 + 24 * 32, data2 + 24 * 32);
  return sad;
}

static unsigned kvz_sad_64x64_avx(const uint8_t *data1, const uint8_t *data2)
{
  unsigned sad = 0;
  sad += kvz_sad_32x32_avx(data1, data2);
  sad += kvz_sad_32x32_avx(data1 + 16 * 64, data2 + 16 * 64);
  sad += kvz_sad_32x32_avx(data1 + 32 * 64, data2 + 32 * 64);
  sad += kvz_sad_32x32_avx(data1 + 48 * 64, data2 + 48 * 64);
  return sad;
}

static unsigned kvz_sad_other_avx(const uint8_t *data1, const uint8_t *data2,
                                  int width, int height,
                                  unsigned stride)
{
  unsigned sad = 0;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      sad += abs(data1[y * stride + x] - data2[y * stride + x]);
    }
  }

  return sad;
}

static unsigned reg_sad_x86_asm(const uint8_t *data1, const uint8_t * data2,
                                const int width, const int height,
                                const unsigned stride1, const unsigned stride2)
{
  if (width == height) {
    if (width == 8) {
      return kvz_sad_8x8_stride_avx(data1, data2, stride1);
    } else if (width == 16) {
      return kvz_sad_16x16_stride_avx(data1, data2, stride1);
    } else if (width == 32) {
      return kvz_sad_32x32_stride_avx(data1, data2, stride1);
    } else if (width == 64) {
      return kvz_sad_64x64_stride_avx(data1, data2, stride1);
    }
  }

  if (width * height >= 16) {
    // Call the vectorized general SAD SSE41 function when the block
    // is big enough to make it worth it.
    return kvz_reg_sad_sse41(data1, data2, width, height, stride1, stride2);
  } else {
    return kvz_sad_other_avx(data1, data2, width, height, stride1);
  }
}

#endif // KVZ_BIT_DEPTH == 8
#endif //defined(KVZ_COMPILE_ASM)

int kvz_strategy_register_picture_x86_asm_avx(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if defined(KVZ_COMPILE_ASM)
#if KVZ_BIT_DEPTH == 8
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "reg_sad", "x86_asm_avx", 30, &reg_sad_x86_asm);

    success &= kvz_strategyselector_register(opaque, "sad_4x4", "x86_asm_avx", 30, &kvz_sad_4x4_avx);
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "x86_asm_avx", 30, &kvz_sad_8x8_avx);
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "x86_asm_avx", 30, &kvz_sad_16x16_avx);
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "x86_asm_avx", 30, &kvz_sad_32x32_avx);
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "x86_asm_avx", 30, &kvz_sad_64x64_avx);

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "x86_asm_avx", 30, &kvz_satd_4x4_avx);
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "x86_asm_avx", 30, &kvz_satd_8x8_avx);
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "x86_asm_avx", 30, &kvz_satd_16x16_avx);
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "x86_asm_avx", 30, &kvz_satd_32x32_avx);
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "x86_asm_avx", 30, &kvz_satd_64x64_avx);
  }
#endif // KVZ_BIT_DEPTH == 8
#endif //!defined(KVZ_COMPILE_ASM)
  return success;
}
