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
#include "strategies/sse41/reg_sad_pow2_widths-sse41.h"

#if COMPILE_INTEL_SSE41
#include <immintrin.h>
#include <stdlib.h>

#include "kvazaar.h"
#include "strategyselector.h"

uint32_t kvz_reg_sad_sse41(const kvz_pixel * const data1, const kvz_pixel * const data2,
                           const int32_t width, const int32_t height, const uint32_t stride1,
                           const uint32_t stride2)
{
  if (width == 0)
    return 0;
  if (width == 4)
    return reg_sad_w4(data1, data2, height, stride1, stride2);
  if (width == 8)
    return reg_sad_w8(data1, data2, height, stride1, stride2);
  if (width == 12)
    return reg_sad_w12(data1, data2, height, stride1, stride2);
  if (width == 16)
    return reg_sad_w16(data1, data2, height, stride1, stride2);
  if (width == 24)
    return reg_sad_w24(data1, data2, height, stride1, stride2);
  else
    return reg_sad_arbitrary(data1, data2, width, height, stride1, stride2);
}

static optimized_sad_func_ptr_t get_optimized_sad_sse41(int32_t width)
{
  if (width == 0)
    return reg_sad_w0;
  if (width == 4)
    return reg_sad_w4;
  if (width == 8)
    return reg_sad_w8;
  if (width == 12)
    return reg_sad_w12;
  if (width == 16)
    return reg_sad_w16;
  if (width == 24)
    return reg_sad_w24;
  else
    return NULL;
}

static uint32_t ver_sad_sse41(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                              int32_t width, int32_t height, uint32_t stride)
{
  if (width == 0)
    return 0;
  if (width == 4)
    return ver_sad_w4(pic_data, ref_data, height, stride);
  if (width == 8)
    return ver_sad_w8(pic_data, ref_data, height, stride);
  if (width == 12)
    return ver_sad_w12(pic_data, ref_data, height, stride);
  if (width == 16)
    return ver_sad_w16(pic_data, ref_data, height, stride);
  else
    return ver_sad_arbitrary(pic_data, ref_data, width, height, stride);
}

static uint32_t hor_sad_sse41(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                              int32_t width, int32_t height, uint32_t pic_stride,
                              uint32_t ref_stride, uint32_t left, uint32_t right)
{
  if (width == 4)
    return hor_sad_sse41_w4(pic_data, ref_data, width, height,
                            pic_stride, ref_stride, left, right);
  if (width == 8)
    return hor_sad_sse41_w8(pic_data, ref_data, width, height,
                            pic_stride, ref_stride, left, right);
  if (width == 16)
    return hor_sad_sse41_w16(pic_data, ref_data, width, height,
                             pic_stride, ref_stride, left, right);
  if (width == 32)
    return hor_sad_sse41_w32(pic_data, ref_data, width, height,
                             pic_stride, ref_stride, left, right);
  else
    return hor_sad_sse41_arbitrary(pic_data, ref_data, width, height,
                                   pic_stride, ref_stride, left, right);
}

#endif //COMPILE_INTEL_SSE41


int kvz_strategy_register_picture_sse41(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_INTEL_SSE41
  if (bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "reg_sad", "sse41", 20, &kvz_reg_sad_sse41);
    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "sse41", 20, &get_optimized_sad_sse41);
    success &= kvz_strategyselector_register(opaque, "ver_sad", "sse41", 20, &ver_sad_sse41);
    success &= kvz_strategyselector_register(opaque, "hor_sad", "sse41", 20, &hor_sad_sse41);
  }
#endif
  return success;
}
