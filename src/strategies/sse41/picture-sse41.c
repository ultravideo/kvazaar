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

static uint32_t hor_sad_sse41_w32(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                  int32_t height, uint32_t pic_stride, uint32_t ref_stride,
                                  uint32_t left, uint32_t right)
{
  const int32_t height_twoline_groups = height & ~1;
  const int32_t height_residual_lines = height &  1;

  __m128i nslo = _mm_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                               8,  9,  10, 11, 12, 13, 14, 15);
  __m128i nshi = _mm_setr_epi8(16, 17, 18, 19, 20, 21, 22, 23,
                               24, 25, 26, 27, 28, 29, 30, 31);
  __m128i epol_masklo, epol_maskhi;
  int32_t border_pix_off;

  if (left) {
    border_pix_off          = left;
    __m128i first_valid_idx = _mm_set1_epi8(left);

    epol_masklo = _mm_cmpgt_epi8(first_valid_idx, nslo);
    epol_maskhi = _mm_cmpgt_epi8(first_valid_idx, nshi);
  } else {
    border_pix_off          = 31 - right;
    __m128i last_valid_idx  = _mm_set1_epi8(border_pix_off);

    epol_masklo = _mm_cmpgt_epi8(nslo, last_valid_idx);
    epol_maskhi = _mm_cmpgt_epi8(nshi, last_valid_idx);
  }

  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;
  for (y = 0; y < height_twoline_groups; y += 2) {
    __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + 0));
    __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + 0));
    __m128i c = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + 16));
    __m128i d = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + 16));
    __m128i e = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * pic_stride + 0));
    __m128i f = _mm_loadu_si128((__m128i *)(ref_data + (y + 1) * ref_stride + 0));
    __m128i g = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * pic_stride + 16));
    __m128i h = _mm_loadu_si128((__m128i *)(ref_data + (y + 1) * ref_stride + 16));

    __m128i border_px_lo = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_pix_off));
    __m128i border_px_hi = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 1) * ref_stride + border_pix_off));
    __m128i b_epol       = _mm_blendv_epi8(b, border_px_lo, epol_masklo);
    __m128i d_epol       = _mm_blendv_epi8(d, border_px_lo, epol_maskhi);
    __m128i f_epol       = _mm_blendv_epi8(f, border_px_hi, epol_masklo);
    __m128i h_epol       = _mm_blendv_epi8(h, border_px_hi, epol_maskhi);

    __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol);
    __m128i curr_sads_ef = _mm_sad_epu8(e, f_epol);
    __m128i curr_sads_gh = _mm_sad_epu8(g, h_epol);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + 0));
    __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + 0));
    __m128i c = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + 16));
    __m128i d = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + 16));

    __m128i border_px = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_pix_off));
    __m128i b_epol    = _mm_blendv_epi8(b, border_px, epol_masklo);
    __m128i d_epol    = _mm_blendv_epi8(d, border_px, epol_maskhi);

    __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
  }

  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static uint32_t hor_sad_sse41(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                              int32_t width, int32_t height, uint32_t pic_stride,
                              uint32_t ref_stride, uint32_t left, uint32_t right)
{
  if (width == 4)
    return hor_sad_sse41_w4(pic_data, ref_data, height,
                            pic_stride, ref_stride, left, right);
  if (width == 8)
    return hor_sad_sse41_w8(pic_data, ref_data, height,
                            pic_stride, ref_stride, left, right);
  if (width == 16)
    return hor_sad_sse41_w16(pic_data, ref_data, height,
                             pic_stride, ref_stride, left, right);
  if (width == 32)
    return hor_sad_sse41_w32(pic_data, ref_data, height,
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
