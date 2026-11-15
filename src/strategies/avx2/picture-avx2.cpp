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

/*
 * \file
 */

#include "strategies/strategies-cpp.h"

extern "C" {
#include "global.h"
#include "kvazaar.h"
#include "strategies/avx2/picture-avx2.h"
#if COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH == 8
#include "strategies/avx2/reg_sad_pow2_widths-avx2.h"
#endif
#include "strategies/strategies-picture.h"
#include "strategyselector.h"
#include "strategies/generic/picture-generic.h"
}

#if COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8

// In a namespace of their own to keep the names apart from the other
// high bit depth functions of this file.
namespace {
namespace hbd_sad {

/*
 * AVX2 versions of the SAD, SSD, variance and bi-prediction functions for
 * 10 and 12-bit pixels (kvz_pixel = uint16_t). Every function except
 * pixel_var computes exactly the same result as the generic version.
 *
 * |a - b| of two pixels fits in an int16 lane, so SAD sums are accumulated in
 * 16-bit lanes and moved to 32-bit lanes with madd every kHbdSadAdds
 * additions, before a lane could exceed INT16_MAX.
 */
constexpr int kHbdSadAdds = 32767 / PIXEL_MAX;
static_assert(kHbdSadAdds >= 8, "SAD accumulation needs at least 8 additions per lane");

KVZ_FORCE_INLINE __m128i hbd_load4(const kvz_pixel *p)
{
  return _mm_loadl_epi64((const __m128i *)p);
}

KVZ_FORCE_INLINE __m128i hbd_load8(const kvz_pixel *p)
{
  return _mm_loadu_si128((const __m128i *)p);
}

KVZ_FORCE_INLINE __m256i hbd_load16(const kvz_pixel *p)
{
  return _mm256_loadu_si256((const __m256i *)p);
}

// Two 4-pixel rows in one 128-bit vector.
KVZ_FORCE_INLINE __m128i hbd_load_2x4(const kvz_pixel *p, unsigned stride)
{
  return _mm_unpacklo_epi64(hbd_load4(p), hbd_load4(p + stride));
}

// Two 8-pixel rows in one 256-bit vector.
KVZ_FORCE_INLINE __m256i hbd_load_2x8(const kvz_pixel *p, unsigned stride)
{
  return _mm256_inserti128_si256(_mm256_castsi128_si256(hbd_load8(p)), hbd_load8(p + stride), 1);
}

// Four 4-pixel rows in one 256-bit vector.
KVZ_FORCE_INLINE __m256i hbd_load_4x4(const kvz_pixel *p, unsigned stride)
{
  return _mm256_inserti128_si256(_mm256_castsi128_si256(hbd_load_2x4(p, stride)),
                                 hbd_load_2x4(p + 2 * stride, stride), 1);
}

KVZ_FORCE_INLINE __m256i hbd_absdiff(__m256i a, __m256i b)
{
  return _mm256_abs_epi16(_mm256_sub_epi16(a, b));
}

KVZ_FORCE_INLINE __m128i hbd_absdiff(__m128i a, __m128i b)
{
  return _mm_abs_epi16(_mm_sub_epi16(a, b));
}

KVZ_FORCE_INLINE __m256i hbd_widen_add(__m256i acc32, __m256i acc16)
{
  return _mm256_add_epi32(acc32, _mm256_madd_epi16(acc16, _mm256_set1_epi16(1)));
}

KVZ_FORCE_INLINE __m128i hbd_widen_add(__m128i acc32, __m128i acc16)
{
  return _mm_add_epi32(acc32, _mm_madd_epi16(acc16, _mm_set1_epi16(1)));
}

KVZ_FORCE_INLINE uint32_t hbd_hsum_epi32(__m128i v)
{
  v = _mm_add_epi32(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2)));
  v = _mm_add_epi32(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(2, 3, 0, 1)));
  return (uint32_t)_mm_cvtsi128_si32(v);
}

KVZ_FORCE_INLINE uint32_t hbd_hsum_epi32(__m256i v)
{
  return hbd_hsum_epi32(_mm_add_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1)));
}

/*
 * SAD of a W-wide block. Rows are processed in pairs: 16-pixel parts are
 * loaded one row at a time, an 8-pixel part as a row pair in one 256-bit
 * vector and a 4-pixel part as a row pair in one 128-bit vector. Only the
 * pixels of the block are read.
 */
template <int W>
uint32_t reg_sad_hbd_w(const kvz_pixel *const data1, const kvz_pixel *const data2,
                       const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  static_assert(W > 0 && W % 4 == 0 && W <= 64, "unsupported width");
  constexpr int k16 = W / 16;
  constexpr bool k8 = (W & 8) != 0;
  constexpr bool k4 = (W & 4) != 0;

  __m256i total = _mm256_setzero_si256();
  __m128i total4 = _mm_setzero_si128();
  int y = 0;

  if constexpr (W == 4) {
    // Four rows per vector.
    const int h4 = height & ~3;
    while (y < h4) {
      const int y_end = h4 - y < 4 * kHbdSadAdds ? h4 : y + 4 * kHbdSadAdds;
      __m256i acc = _mm256_setzero_si256();
      for (; y < y_end; y += 4) {
        acc = _mm256_add_epi16(acc, hbd_absdiff(hbd_load_4x4(data1 + y * stride1, stride1),
                                                hbd_load_4x4(data2 + y * stride2, stride2)));
      }
      total = hbd_widen_add(total, acc);
    }
    __m128i acc = _mm_setzero_si128();
    if (height - y >= 2) {
      acc = hbd_absdiff(hbd_load_2x4(data1 + y * stride1, stride1),
                        hbd_load_2x4(data2 + y * stride2, stride2));
      y += 2;
    }
    if (y < height) {
      acc = _mm_add_epi16(acc, hbd_absdiff(hbd_load4(data1 + y * stride1),
                                           hbd_load4(data2 + y * stride2)));
    }
    total4 = hbd_widen_add(total4, acc);
  } else {
    // Additions to a lane of acc per row pair.
    constexpr int kAddsPerPair = 2 * k16 + (k8 ? 1 : 0);
    constexpr int kPairs = kHbdSadAdds / kAddsPerPair;
    static_assert(kPairs >= 1, "too wide for 16-bit accumulation");

    const int h2 = height & ~1;
    while (y < h2) {
      const int y_end = h2 - y < 2 * kPairs ? h2 : y + 2 * kPairs;
      __m256i acc = _mm256_setzero_si256();
      __m128i acc4 = _mm_setzero_si128();
      for (; y < y_end; y += 2) {
        const kvz_pixel *a = data1 + y * stride1;
        const kvz_pixel *b = data2 + y * stride2;
        KVZ_UNROLL
        for (int i = 0; i < k16; ++i) {
          acc = _mm256_add_epi16(acc, hbd_absdiff(hbd_load16(a + 16 * i), hbd_load16(b + 16 * i)));
          acc = _mm256_add_epi16(acc, hbd_absdiff(hbd_load16(a + stride1 + 16 * i),
                                                  hbd_load16(b + stride2 + 16 * i)));
        }
        if constexpr (k8) {
          acc = _mm256_add_epi16(acc, hbd_absdiff(hbd_load_2x8(a + 16 * k16, stride1),
                                                  hbd_load_2x8(b + 16 * k16, stride2)));
        }
        if constexpr (k4) {
          acc4 = _mm_add_epi16(acc4, hbd_absdiff(hbd_load_2x4(a + W - 4, stride1),
                                                 hbd_load_2x4(b + W - 4, stride2)));
        }
      }
      total = hbd_widen_add(total, acc);
      if constexpr (k4) {
        total4 = hbd_widen_add(total4, acc4);
      }
    }
    if (y < height) {
      // Last row of an odd height.
      const kvz_pixel *a = data1 + y * stride1;
      const kvz_pixel *b = data2 + y * stride2;
      if constexpr (k16 > 0) {
        __m256i acc = _mm256_setzero_si256();
        KVZ_UNROLL
        for (int i = 0; i < k16; ++i) {
          acc = _mm256_add_epi16(acc, hbd_absdiff(hbd_load16(a + 16 * i), hbd_load16(b + 16 * i)));
        }
        total = hbd_widen_add(total, acc);
      }
      if constexpr (k8 || k4) {
        __m128i acc = _mm_setzero_si128();
        if constexpr (k8) {
          acc = hbd_absdiff(hbd_load8(a + 16 * k16), hbd_load8(b + 16 * k16));
        }
        if constexpr (k4) {
          acc = _mm_add_epi16(acc, hbd_absdiff(hbd_load4(a + W - 4), hbd_load4(b + W - 4)));
        }
        total4 = hbd_widen_add(total4, acc);
      }
    }
  }
  return hbd_hsum_epi32(total) + hbd_hsum_epi32(total4);
}

uint32_t reg_sad_hbd_w0(const kvz_pixel *const, const kvz_pixel *const,
                        const int32_t, const uint32_t, const uint32_t)
{
  return 0;
}

/*
 * SAD of a block of any width. With kHor, every row of data1 is compared
 * with the single pixel data2[y * stride2] (hor_sad). Each vector of absolute
 * differences is widened to 32 bits right away.
 */
template <bool kHor>
uint32_t sad_hbd_any(const kvz_pixel *data1, const kvz_pixel *data2,
                     int width, int height, unsigned stride1, unsigned stride2)
{
  __m256i acc = _mm256_setzero_si256();
  __m128i acc4 = _mm_setzero_si128();
  uint32_t tail = 0;
  for (int y = 0; y < height; ++y) {
    const kvz_pixel *a = data1 + y * stride1;
    const kvz_pixel *b = data2 + y * stride2;
    const __m256i ref16 = _mm256_set1_epi16(kHor ? (short)b[0] : 0);
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      acc = hbd_widen_add(acc, hbd_absdiff(hbd_load16(a + x), kHor ? ref16 : hbd_load16(b + x)));
    }
    if (x + 8 <= width) {
      acc4 = hbd_widen_add(acc4, hbd_absdiff(hbd_load8(a + x), kHor ? _mm256_castsi256_si128(ref16)
                                                                    : hbd_load8(b + x)));
      x += 8;
    }
    if (x + 4 <= width) {
      // With kHor the upper half of the reference vector is not zero, so
      // zero the result instead.
      __m128i d = hbd_absdiff(hbd_load4(a + x), kHor ? _mm256_castsi256_si128(ref16) : hbd_load4(b + x));
      acc4 = hbd_widen_add(acc4, _mm_move_epi64(d));
      x += 4;
    }
    for (; x < width; ++x) {
      tail += abs((int)a[x] - (int)(kHor ? b[0] : b[x]));
    }
  }
  return hbd_hsum_epi32(acc) + hbd_hsum_epi32(acc4) + tail;
}

uint32_t reg_sad_hbd_avx2(const kvz_pixel *const data1, const kvz_pixel *const data2,
                          const int width, const int height,
                          const unsigned stride1, const unsigned stride2)
{
  switch (width) {
    case 0:  return 0;
    case 4:  return reg_sad_hbd_w<4>(data1, data2, height, stride1, stride2);
    case 8:  return reg_sad_hbd_w<8>(data1, data2, height, stride1, stride2);
    case 12: return reg_sad_hbd_w<12>(data1, data2, height, stride1, stride2);
    case 16: return reg_sad_hbd_w<16>(data1, data2, height, stride1, stride2);
    case 24: return reg_sad_hbd_w<24>(data1, data2, height, stride1, stride2);
    case 32: return reg_sad_hbd_w<32>(data1, data2, height, stride1, stride2);
    case 48: return reg_sad_hbd_w<48>(data1, data2, height, stride1, stride2);
    case 64: return reg_sad_hbd_w<64>(data1, data2, height, stride1, stride2);
    default: return sad_hbd_any<false>(data1, data2, width, height, stride1, stride2);
  }
}

optimized_sad_func_ptr_t get_optimized_sad_hbd_avx2(int32_t width)
{
  switch (width) {
    case 0:  return reg_sad_hbd_w0;
    case 4:  return reg_sad_hbd_w<4>;
    case 8:  return reg_sad_hbd_w<8>;
    case 12: return reg_sad_hbd_w<12>;
    case 16: return reg_sad_hbd_w<16>;
    case 24: return reg_sad_hbd_w<24>;
    case 32: return reg_sad_hbd_w<32>;
    case 48: return reg_sad_hbd_w<48>;
    case 64: return reg_sad_hbd_w<64>;
    default: return NULL;
  }
}

template <int N>
unsigned sad_nxn_hbd_avx2(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return reg_sad_hbd_w<N>(block1, block2, N, N, N) >> (KVZ_BIT_DEPTH - 8);
}

template <int N>
void sad_nxn_dual_hbd_avx2(const pred_buffer preds, const kvz_pixel *const orig,
                           unsigned /* num_modes */, unsigned *costs_out)
{
  costs_out[0] = reg_sad_hbd_w<N>(preds[0], orig, N, N, N) >> (KVZ_BIT_DEPTH - 8);
  costs_out[1] = reg_sad_hbd_w<N>(preds[1], orig, N, N, N) >> (KVZ_BIT_DEPTH - 8);
}

// The reference row is the same for every row of the block.
uint32_t ver_sad_hbd_avx2(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                          int32_t width, int32_t height, uint32_t stride)
{
  return reg_sad_hbd_avx2(pic_data, ref_data, width, height, stride, 0);
}

uint32_t hor_sad_hbd_avx2(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                          int32_t width, int32_t height, uint32_t pic_stride,
                          uint32_t ref_stride, uint32_t left, uint32_t right)
{
  uint32_t result = 0;
  if (left) {
    result += sad_hbd_any<true>(pic_data, ref_data + left, left, height, pic_stride, ref_stride);
    result += reg_sad_hbd_avx2(pic_data + left, ref_data + left, width - left,
                               height, pic_stride, ref_stride);
  } else if (right) {
    result += reg_sad_hbd_avx2(pic_data, ref_data, width - right,
                               height, pic_stride, ref_stride);
    result += sad_hbd_any<true>(pic_data + width - right, ref_data + width - right - 1,
                                right, height, pic_stride, ref_stride);
  } else {
    result += reg_sad_hbd_avx2(pic_data, ref_data, width, height, pic_stride, ref_stride);
  }
  return result;
}

/*
 * SSD of a W x W block. madd of a difference with itself adds two squares,
 * at most 2 * PIXEL_MAX^2, to a 32-bit lane; the lanes are moved to 64-bit
 * lanes before they could overflow as unsigned.
 */
constexpr uint32_t kHbdSsdAdds = 0xFFFFFFFFu / (2u * PIXEL_MAX * PIXEL_MAX);

KVZ_FORCE_INLINE __m256i hbd_sqdiff(__m256i a, __m256i b)
{
  const __m256i d = _mm256_sub_epi16(a, b);
  return _mm256_madd_epi16(d, d);
}

KVZ_FORCE_INLINE __m256i hbd_widen_add64(__m256i acc64, __m256i acc32)
{
  const __m256i zero = _mm256_setzero_si256();
  acc64 = _mm256_add_epi64(acc64, _mm256_unpacklo_epi32(acc32, zero));
  return _mm256_add_epi64(acc64, _mm256_unpackhi_epi32(acc32, zero));
}

template <int W>
unsigned pixels_calc_ssd_hbd_w(const kvz_pixel *const ref, const kvz_pixel *const rec,
                               const int ref_stride, const int rec_stride)
{
  // Rows per vector and vectors per row.
  constexpr int kRowsPerVec = W >= 16 ? 1 : 16 / W;
  constexpr int kVecsPerRow = W >= 16 ? W / 16 : 1;
  constexpr int kRowsPerFlush = (int)(kHbdSsdAdds / kVecsPerRow) * kRowsPerVec;
  static_assert(kRowsPerFlush >= kRowsPerVec, "SSD accumulation overflows");

  __m256i acc64 = _mm256_setzero_si256();
  for (int y0 = 0; y0 < W; y0 += kRowsPerFlush) {
    const int y1 = W - y0 < kRowsPerFlush ? W : y0 + kRowsPerFlush;
    __m256i acc = _mm256_setzero_si256();
    for (int y = y0; y < y1; y += kRowsPerVec) {
      const kvz_pixel *a = ref + y * ref_stride;
      const kvz_pixel *b = rec + y * rec_stride;
      if constexpr (W == 4) {
        acc = _mm256_add_epi32(acc, hbd_sqdiff(hbd_load_4x4(a, ref_stride), hbd_load_4x4(b, rec_stride)));
      } else if constexpr (W == 8) {
        acc = _mm256_add_epi32(acc, hbd_sqdiff(hbd_load_2x8(a, ref_stride), hbd_load_2x8(b, rec_stride)));
      } else {
        KVZ_UNROLL
        for (int i = 0; i < kVecsPerRow; ++i) {
          acc = _mm256_add_epi32(acc, hbd_sqdiff(hbd_load16(a + 16 * i), hbd_load16(b + 16 * i)));
        }
      }
    }
    acc64 = hbd_widen_add64(acc64, acc);
  }
  __m128i s = _mm_add_epi64(_mm256_castsi256_si128(acc64), _mm256_extracti128_si256(acc64, 1));
  s = _mm_add_epi64(s, _mm_unpackhi_epi64(s, s));
  uint64_t ssd;
  _mm_storel_epi64((__m128i *)&ssd, s);
  return (unsigned)(ssd >> (2 * (KVZ_BIT_DEPTH - 8)));
}

unsigned pixels_calc_ssd_hbd_avx2(const kvz_pixel *const ref, const kvz_pixel *const rec,
                                  const int ref_stride, const int rec_stride, const int width)
{
  switch (width) {
    case 4:  return pixels_calc_ssd_hbd_w<4>(ref, rec, ref_stride, rec_stride);
    case 8:  return pixels_calc_ssd_hbd_w<8>(ref, rec, ref_stride, rec_stride);
    case 16: return pixels_calc_ssd_hbd_w<16>(ref, rec, ref_stride, rec_stride);
    case 32: return pixels_calc_ssd_hbd_w<32>(ref, rec, ref_stride, rec_stride);
    case 64: return pixels_calc_ssd_hbd_w<64>(ref, rec, ref_stride, rec_stride);
    default: break;
  }
  uint64_t ssd = 0;
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      const int diff = ref[x + y * ref_stride] - rec[x + y * rec_stride];
      ssd += (uint32_t)(diff * diff);
    }
  }
  return (unsigned)(ssd >> (2 * (KVZ_BIT_DEPTH - 8)));
}

/*
 * Bi-prediction average. The generic code computes
 * clip((s0 + s1 + offset) >> shift) with shift = 15 - KVZ_BIT_DEPTH,
 * offset = 1 << (shift - 1) and samples s0, s1 at 14-bit precision.
 * - Two pixels: s = p << (shift - 1), so the result is (p0 + p1 + 1) >> 1,
 *   which is avg_epu16 and needs no clipping.
 * - Otherwise the sums are computed with saturating 16-bit additions. When a
 *   sum saturates to INT16_MAX, the result 32767 >> shift equals PIXEL_MAX,
 *   which is also what clipping the exact result gives; a sum saturated to
 *   INT16_MIN gives a negative result, which clips to 0 like the exact one.
 */
constexpr int kBipredShift = 15 - KVZ_BIT_DEPTH;
static_assert((32767 >> kBipredShift) == PIXEL_MAX, "saturation must clip to PIXEL_MAX");

enum { BIPRED_PX_PX, BIPRED_IM_IM, BIPRED_PX_IM };

template <int MODE>
KVZ_FORCE_INLINE __m256i bipred_hbd_combine(__m256i a, __m256i b)
{
  if constexpr (MODE == BIPRED_PX_PX) {
    return _mm256_avg_epu16(a, b);
  } else {
    if constexpr (MODE == BIPRED_PX_IM) {
      a = _mm256_slli_epi16(a, 14 - KVZ_BIT_DEPTH);
    }
    __m256i s = _mm256_adds_epi16(_mm256_adds_epi16(a, b), _mm256_set1_epi16(1 << (kBipredShift - 1)));
    s = _mm256_srai_epi16(s, kBipredShift);
    return _mm256_min_epi16(_mm256_max_epi16(s, _mm256_setzero_si256()), _mm256_set1_epi16(PIXEL_MAX));
  }
}

template <int MODE>
KVZ_FORCE_INLINE __m128i bipred_hbd_combine(__m128i a, __m128i b)
{
  if constexpr (MODE == BIPRED_PX_PX) {
    return _mm_avg_epu16(a, b);
  } else {
    if constexpr (MODE == BIPRED_PX_IM) {
      a = _mm_slli_epi16(a, 14 - KVZ_BIT_DEPTH);
    }
    __m128i s = _mm_adds_epi16(_mm_adds_epi16(a, b), _mm_set1_epi16(1 << (kBipredShift - 1)));
    s = _mm_srai_epi16(s, kBipredShift);
    return _mm_min_epi16(_mm_max_epi16(s, _mm_setzero_si128()), _mm_set1_epi16(PIXEL_MAX));
  }
}

// The generic formula for one pixel.
template <int MODE, typename T0, typename T1>
KVZ_FORCE_INLINE kvz_pixel bipred_hbd_scalar(T0 a, T1 b)
{
  const int32_t offset = 1 << (kBipredShift - 1);
  const int16_t s0 = MODE == BIPRED_IM_IM ? (int16_t)a : (int16_t)(a << (14 - KVZ_BIT_DEPTH));
  const int16_t s1 = MODE == BIPRED_PX_PX ? (int16_t)(b << (14 - KVZ_BIT_DEPTH)) : (int16_t)b;
  return kvz_fast_clip_32bit_to_pixel((s0 + s1 + offset) >> kBipredShift);
}

/*
 * src0 and src1 are w x h blocks with a stride of w. The 4-pixel wide blocks
 * are processed four rows (16 contiguous source samples) at a time.
 */
template <int MODE, typename T0, typename T1>
void bipred_average_hbd_block(kvz_pixel *dst, const T0 *src0, const T1 *src1,
                              unsigned w, unsigned h, unsigned dst_stride)
{
  unsigned y = 0;
  if (w == 4) {
    for (; y + 4 <= h; y += 4) {
      const __m256i r = bipred_hbd_combine<MODE>(hbd_load16((const kvz_pixel *)src0 + 4 * y),
                                                 hbd_load16((const kvz_pixel *)src1 + 4 * y));
      const __m128i lo = _mm256_castsi256_si128(r);
      const __m128i hi = _mm256_extracti128_si256(r, 1);
      _mm_storel_epi64((__m128i *)(dst + (y + 0) * dst_stride), lo);
      _mm_storeh_pd((double *)(dst + (y + 1) * dst_stride), _mm_castsi128_pd(lo));
      _mm_storel_epi64((__m128i *)(dst + (y + 2) * dst_stride), hi);
      _mm_storeh_pd((double *)(dst + (y + 3) * dst_stride), _mm_castsi128_pd(hi));
    }
  }
  for (; y < h; ++y) {
    const T0 *a = src0 + y * w;
    const T1 *b = src1 + y * w;
    kvz_pixel *d = dst + y * dst_stride;
    unsigned x = 0;
    for (; x + 16 <= w; x += 16) {
      _mm256_storeu_si256((__m256i *)(d + x),
                          bipred_hbd_combine<MODE>(hbd_load16((const kvz_pixel *)a + x),
                                                   hbd_load16((const kvz_pixel *)b + x)));
    }
    if (x + 8 <= w) {
      _mm_storeu_si128((__m128i *)(d + x),
                       bipred_hbd_combine<MODE>(hbd_load8((const kvz_pixel *)a + x),
                                                hbd_load8((const kvz_pixel *)b + x)));
      x += 8;
    }
    if (x + 4 <= w) {
      _mm_storel_epi64((__m128i *)(d + x),
                       bipred_hbd_combine<MODE>(hbd_load4((const kvz_pixel *)a + x),
                                                hbd_load4((const kvz_pixel *)b + x)));
      x += 4;
    }
    for (; x < w; ++x) {
      d[x] = bipred_hbd_scalar<MODE>(a[x], b[x]);
    }
  }
}

void bipred_average_hbd_avx2(lcu_t *const lcu,
                             const yuv_t *const px_L0,
                             const yuv_t *const px_L1,
                             const yuv_im_t *const im_L0,
                             const yuv_im_t *const im_L1,
                             const unsigned pu_x,
                             const unsigned pu_y,
                             const unsigned pu_w,
                             const unsigned pu_h,
                             const unsigned im_flags_L0,
                             const unsigned im_flags_L1,
                             const bool predict_luma,
                             const bool predict_chroma,
                             const uint8_t chroma_shift_w,
                             const uint8_t chroma_shift_h)
{
  if (predict_luma) {
    const unsigned pb_offset = SUB_SCU(pu_y) * LCU_WIDTH + SUB_SCU(pu_x);
    kvz_pixel *dst = lcu->rec.y + pb_offset;

    if (!(im_flags_L0 & 1) && !(im_flags_L1 & 1)) {
      bipred_average_hbd_block<BIPRED_PX_PX>(dst, px_L0->y, px_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else if ((im_flags_L0 & 1) && (im_flags_L1 & 1)) {
      bipred_average_hbd_block<BIPRED_IM_IM>(dst, im_L0->y, im_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else {
      const kvz_pixel    *src_px = (im_flags_L0 & 1) ? px_L1->y : px_L0->y;
      const kvz_pixel_im *src_im = (im_flags_L0 & 1) ? im_L0->y : im_L1->y;
      bipred_average_hbd_block<BIPRED_PX_IM>(dst, src_px, src_im, pu_w, pu_h, LCU_WIDTH);
    }
  }
  if (predict_chroma) {
    const unsigned pb_offset = ((SUB_SCU(pu_y) >> chroma_shift_h) * (LCU_WIDTH >> chroma_shift_w)) +
                               (SUB_SCU(pu_x) >> chroma_shift_w);
    const unsigned pb_w = pu_w >> chroma_shift_w;
    const unsigned pb_h = pu_h >> chroma_shift_h;
    const unsigned stride = LCU_WIDTH >> chroma_shift_w;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_hbd_block<BIPRED_PX_PX>(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, stride);
      bipred_average_hbd_block<BIPRED_PX_PX>(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, stride);
    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_hbd_block<BIPRED_IM_IM>(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, stride);
      bipred_average_hbd_block<BIPRED_IM_IM>(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, stride);
    } else {
      const kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      const kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      const kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      const kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_hbd_block<BIPRED_PX_IM>(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, stride);
      bipred_average_hbd_block<BIPRED_PX_IM>(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, stride);
    }
  }
}

/*
 * Variance as (n * sum(x^2) - sum(x)^2) / n^2. The sums are exact. The
 * numerator is computed exactly in 64 bits when it fits, otherwise in double;
 * either way the result agrees with the generic two-pass loop to within its
 * rounding error.
 */
double pixel_var_hbd_finish(uint64_t sum, uint64_t sumsq, uint32_t len)
{
  const double n = (double)len;
  if (len <= 65536) {
    // len * sumsq <= 2^16 * 2^16 * PIXEL_MAX^2 < 2^64.
    const uint64_t num = (uint64_t)len * sumsq - sum * sum;
    return (double)num / (n * n);
  }
  const double mean = (double)sum / n;
  const double var = ((double)sumsq - mean * (double)sum) / n;
  return var < 0 ? 0 : var;
}

double pixel_var_hbd_avx2(const kvz_pixel *arr, const uint32_t len)
{
  // A madd lane of the squares takes at most 2 * PIXEL_MAX^2 per vector.
  constexpr uint32_t kBlock = 16 * (kHbdSsdAdds < 1024 ? kHbdSsdAdds : 1024);
  const __m256i ones = _mm256_set1_epi16(1);
  __m256i sum64 = _mm256_setzero_si256();
  __m256i sq64 = _mm256_setzero_si256();
  uint32_t i = 0;
  while (i + 16 <= len) {
    const uint32_t end = len - i < kBlock ? i + ((len - i) & ~15u) : i + kBlock;
    __m256i sum32 = _mm256_setzero_si256();
    __m256i sq32 = _mm256_setzero_si256();
    for (; i < end; i += 16) {
      const __m256i v = hbd_load16(arr + i);
      sum32 = _mm256_add_epi32(sum32, _mm256_madd_epi16(v, ones));
      sq32 = _mm256_add_epi32(sq32, _mm256_madd_epi16(v, v));
    }
    sum64 = hbd_widen_add64(sum64, sum32);
    sq64 = hbd_widen_add64(sq64, sq32);
  }
  ALIGNED(32) uint64_t s[4], q[4];
  _mm256_store_si256((__m256i *)s, sum64);
  _mm256_store_si256((__m256i *)q, sq64);
  uint64_t sum = s[0] + s[1] + s[2] + s[3];
  uint64_t sumsq = q[0] + q[1] + q[2] + q[3];
  for (; i < len; ++i) {
    sum += arr[i];
    sumsq += (uint32_t)arr[i] * arr[i];
  }
  return pixel_var_hbd_finish(sum, sumsq, len);
}

} // namespace hbd_sad
} // namespace

#endif // COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8

#if COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8
/*
 * SATD and temporal filter strategies for 10 and 12-bit pixels (kvz_pixel is
 * uint16_t). The results are identical to the generic functions.
 *
 * The 2D Hadamard transform of an 8x8 block multiplies the magnitude by 8 in
 * each direction. The differences have at most KVZ_BIT_DEPTH + 1 bits, so the
 * first (vertical) pass fits 16-bit lanes at 10 and 12 bits (at most 32760).
 * The last butterfly stage feeds straight into the absolute sum and
 * |a + b| + |a - b| == 2 * max(|a|, |b|) for integers, so it is replaced by a
 * maximum. At 10 bits the two remaining horizontal stages also fit 16 bits
 * (at most 32736). At 12 bits they are done in 32-bit lanes.
 */
namespace {


static INLINE __m256i hbd_load_2x128(const kvz_pixel *lo, const kvz_pixel *hi)
{
  return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)lo)),
                                 _mm_loadu_si128((const __m128i *)hi), 1);
}

// Two rows of four pixels in a 128-bit vector.
static INLINE __m128i hbd_load_2x64(const kvz_pixel *a, const kvz_pixel *b)
{
  return _mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)a), _mm_loadl_epi64((const __m128i *)b));
}

static INLINE __m256i hbd_load_4x64(const kvz_pixel *a, const kvz_pixel *b,
                                    const kvz_pixel *c, const kvz_pixel *d)
{
  return _mm256_inserti128_si256(_mm256_castsi128_si256(hbd_load_2x64(a, b)), hbd_load_2x64(c, d), 1);
}

// Sums of the 32-bit lanes of each 128-bit half, in lane 0 of the half.
static INLINE __m256i hbd_hsum_halves_epi32(__m256i v)
{
  v = _mm256_add_epi32(v, _mm256_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2)));
  return _mm256_add_epi32(v, _mm256_shuffle_epi32(v, _MM_SHUFFLE(2, 3, 0, 1)));
}

static INLINE uint32_t hbd_hsum_epi32(__m256i v)
{
  v = hbd_hsum_halves_epi32(v);
  return (uint32_t)_mm_cvtsi128_si32(_mm_add_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1)));
}

// Three butterfly stages between the eight registers (rows).
static INLINE void hbd_hadamard8_ver(__m256i r[8])
{
  const __m256i a0 = _mm256_add_epi16(r[0], r[4]), a4 = _mm256_sub_epi16(r[0], r[4]);
  const __m256i a1 = _mm256_add_epi16(r[1], r[5]), a5 = _mm256_sub_epi16(r[1], r[5]);
  const __m256i a2 = _mm256_add_epi16(r[2], r[6]), a6 = _mm256_sub_epi16(r[2], r[6]);
  const __m256i a3 = _mm256_add_epi16(r[3], r[7]), a7 = _mm256_sub_epi16(r[3], r[7]);

  const __m256i b0 = _mm256_add_epi16(a0, a2), b2 = _mm256_sub_epi16(a0, a2);
  const __m256i b1 = _mm256_add_epi16(a1, a3), b3 = _mm256_sub_epi16(a1, a3);
  const __m256i b4 = _mm256_add_epi16(a4, a6), b6 = _mm256_sub_epi16(a4, a6);
  const __m256i b5 = _mm256_add_epi16(a5, a7), b7 = _mm256_sub_epi16(a5, a7);

  r[0] = _mm256_add_epi16(b0, b1); r[1] = _mm256_sub_epi16(b0, b1);
  r[2] = _mm256_add_epi16(b2, b3); r[3] = _mm256_sub_epi16(b2, b3);
  r[4] = _mm256_add_epi16(b4, b5); r[5] = _mm256_sub_epi16(b4, b5);
  r[6] = _mm256_add_epi16(b6, b7); r[7] = _mm256_sub_epi16(b6, b7);
}

#if KVZ_BIT_DEPTH <= 10
/*
 * Horizontal 8-point transform of the rows in x and y (one row in each
 * 128-bit half), fused with the absolute sum, in 16-bit lanes. The last stage
 * pairs adjacent lanes and is replaced by their maximum. Returns 32-bit lanes
 * whose sum per half is half of the absolute sum of the rows of that half.
 */
static INLINE __m256i hbd_hor8_absmax_16(__m256i x, __m256i y)
{
  const __m256i s1 = _mm256_hadd_epi16(x, y);
  const __m256i d1 = _mm256_hsub_epi16(x, y);
  const __m256i s2 = _mm256_abs_epi16(_mm256_hadd_epi16(s1, d1));
  const __m256i d2 = _mm256_abs_epi16(_mm256_hsub_epi16(s1, d1));
  // The maximum of each pair ends up in the low 16 bits of the 32-bit lane.
  const __m256i ms = _mm256_max_epi16(s2, _mm256_srli_epi32(s2, 16));
  const __m256i md = _mm256_max_epi16(d2, _mm256_srli_epi32(d2, 16));
  // Each maximum is at most 32736, so their sum fits 16 unsigned bits.
  return _mm256_and_si256(_mm256_add_epi16(ms, md), _mm256_set1_epi32(0xFFFF));
}

#else
/*
 * The same in 32-bit lanes, for one register of two rows. The madds widen the
 * row and do the first stage. Returns 32-bit lanes whose sum per half is the
 * absolute sum of the row of that half (each maximum is counted twice).
 */
static INLINE __m256i hbd_hor8_absmax_32(__m256i x)
{
  const __m256i p = _mm256_madd_epi16(x, _mm256_set1_epi16(1));
  const __m256i q = _mm256_madd_epi16(x, _mm256_set1_epi32(-65535)); // (1, -1) pairs
  const __m256i s = _mm256_abs_epi32(_mm256_hadd_epi32(p, q));
  const __m256i d = _mm256_abs_epi32(_mm256_hsub_epi32(p, q));
  const __m256i ms = _mm256_max_epi32(s, _mm256_shuffle_epi32(s, _MM_SHUFFLE(2, 3, 0, 1)));
  const __m256i md = _mm256_max_epi32(d, _mm256_shuffle_epi32(d, _MM_SHUFFLE(2, 3, 0, 1)));
  return _mm256_add_epi32(ms, md);
}
#endif

/*
 * Horizontal pass and absolute sum of n registers of two vertically
 * transformed rows each. Returns 32-bit lanes whose sum per half is the
 * absolute sum of the transform of the rows in that half.
 */
template <int n>
static INLINE __m256i hbd_hor8_abs_sum(const __m256i *r)
{
#if KVZ_BIT_DEPTH <= 10
  // The horizontal pass fits 16-bit lanes.
  __m256i acc = hbd_hor8_absmax_16(r[0], r[1]);
  for (int i = 2; i < n; i += 2) {
    acc = _mm256_add_epi32(acc, hbd_hor8_absmax_16(r[i], r[i + 1]));
  }
  return _mm256_add_epi32(acc, acc);
#else
  __m256i acc = hbd_hor8_absmax_32(r[0]);
  for (int i = 1; i < n; ++i) {
    acc = _mm256_add_epi32(acc, hbd_hor8_absmax_32(r[i]));
  }
  return acc;
#endif
}

/*
 * SATD of two 8x8 blocks: d[i] holds row i of the differences of the first
 * block in the low half and of the second block in the high half.
 */
static INLINE void hbd_satd_8x8_dual_core(__m256i d[8], unsigned *satd0, unsigned *satd1)
{
  hbd_hadamard8_ver(d);
  const __m256i sum = hbd_hsum_halves_epi32(hbd_hor8_abs_sum<8>(d));
  // (sum + 2) >> 2 as in satd_8x8_subblock_generic.
  *satd0 = ((uint32_t)_mm_cvtsi128_si32(_mm256_castsi256_si128(sum)) + 2) >> 2;
  *satd1 = ((uint32_t)_mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 1)) + 2) >> 2;
}

/* SATD of one 8x8 block: d[i] holds row i in the low half and row i + 4 in
 * the high half. */
static INLINE unsigned hbd_satd_8x8_core(const __m256i d[4])
{
  // Stages pairing rows (i, i ^ 1) and (i, i ^ 2) in both halves.
  const __m256i a0 = _mm256_add_epi16(d[0], d[1]), a1 = _mm256_sub_epi16(d[0], d[1]);
  const __m256i a2 = _mm256_add_epi16(d[2], d[3]), a3 = _mm256_sub_epi16(d[2], d[3]);
  const __m256i b0 = _mm256_add_epi16(a0, a2), b2 = _mm256_sub_epi16(a0, a2);
  const __m256i b1 = _mm256_add_epi16(a1, a3), b3 = _mm256_sub_epi16(a1, a3);
  // The stage pairing rows i and i + 4, which are in different halves.
  const __m256i p0 = _mm256_permute2x128_si256(b0, b1, 0x20);
  const __m256i q0 = _mm256_permute2x128_si256(b0, b1, 0x31);
  const __m256i p1 = _mm256_permute2x128_si256(b2, b3, 0x20);
  const __m256i q1 = _mm256_permute2x128_si256(b2, b3, 0x31);
  __m256i r[4];
  r[0] = _mm256_add_epi16(p0, q0); r[1] = _mm256_sub_epi16(p0, q0);
  r[2] = _mm256_add_epi16(p1, q1); r[3] = _mm256_sub_epi16(p1, q1);
  return (hbd_hsum_epi32(hbd_hor8_abs_sum<4>(r)) + 2) >> 2;
}

static INLINE unsigned satd_8x8_subblock_hbd_avx2(const kvz_pixel *buf1, int32_t stride1,
                                                  const kvz_pixel *buf2, int32_t stride2)
{
  __m256i d[4];
  for (int i = 0; i < 4; ++i) {
    d[i] = _mm256_sub_epi16(hbd_load_2x128(buf1 + i * stride1, buf1 + (i + 4) * stride1),
                            hbd_load_2x128(buf2 + i * stride2, buf2 + (i + 4) * stride2));
  }
  return hbd_satd_8x8_core(d);
}

/* Sum of the SATDs of the two horizontally adjacent 8x8 blocks of a 16x8 area. */
static INLINE unsigned satd_16x8_subblocks_hbd_avx2(const kvz_pixel *buf1, int32_t stride1,
                                                    const kvz_pixel *buf2, int32_t stride2)
{
  __m256i d[8];
  for (int i = 0; i < 8; ++i) {
    d[i] = _mm256_sub_epi16(_mm256_loadu_si256((const __m256i *)(buf1 + i * stride1)),
                            _mm256_loadu_si256((const __m256i *)(buf2 + i * stride2)));
  }
  unsigned s0, s1;
  hbd_satd_8x8_dual_core(d, &s0, &s1);
  return s0 + s1;
}

/* SATDs of the 8x8 blocks of two predictions against the same original block. */
static INLINE void satd_8x8_subblock_dual_hbd_avx2(const kvz_pixel *pred0, const kvz_pixel *pred1, int32_t stride,
                                                   const kvz_pixel *orig, int32_t orig_stride,
                                                   unsigned *satd0, unsigned *satd1)
{
  __m256i d[8];
  for (int i = 0; i < 8; ++i) {
    const __m256i o = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(orig + i * orig_stride)));
    d[i] = _mm256_sub_epi16(hbd_load_2x128(pred0 + i * stride, pred1 + i * stride), o);
  }
  hbd_satd_8x8_dual_core(d, satd0, satd1);
}

/*
 * 4x4 SATD of two blocks: lo holds rows 0 and 1 and hi rows 2 and 3 of the
 * differences, the first block in the low halves and the second in the high
 * halves. The horizontal and vertical butterflies are hadd/hsub pairs, at most
 * 32760 before the last stage, which is replaced by a maximum. Returns half of
 * the absolute sum of each block, which is what hadamard_4x4_generic returns.
 */
static INLINE __m256i hbd_satd_4x4_dual_core(__m256i lo, __m256i hi)
{
  const __m256i r0 = _mm256_hadd_epi16(lo, hi);
  const __m256i r1 = _mm256_hsub_epi16(lo, hi);
  const __m256i r2 = _mm256_hadd_epi16(r0, r1);
  const __m256i r3 = _mm256_hsub_epi16(r0, r1);
  const __m256i s = _mm256_abs_epi16(_mm256_hadd_epi16(r2, r3));
  const __m256i d = _mm256_abs_epi16(_mm256_hsub_epi16(r2, r3));
  const __m256i ms = _mm256_max_epi16(s, _mm256_srli_epi32(s, 16));
  const __m256i md = _mm256_max_epi16(d, _mm256_srli_epi32(d, 16));
  const __m256i m = _mm256_and_si256(_mm256_add_epi16(ms, md), _mm256_set1_epi32(0xFFFF));
  return hbd_hsum_halves_epi32(m);
}

static INLINE unsigned hbd_satd_4x4_core(__m128i lo, __m128i hi)
{
  const __m128i r0 = _mm_hadd_epi16(lo, hi);
  const __m128i r1 = _mm_hsub_epi16(lo, hi);
  const __m128i r2 = _mm_hadd_epi16(r0, r1);
  const __m128i r3 = _mm_hsub_epi16(r0, r1);
  const __m128i s = _mm_abs_epi16(_mm_hadd_epi16(r2, r3));
  const __m128i d = _mm_abs_epi16(_mm_hsub_epi16(r2, r3));
  const __m128i ms = _mm_max_epi16(s, _mm_srli_epi32(s, 16));
  const __m128i md = _mm_max_epi16(d, _mm_srli_epi32(d, 16));
  __m128i m = _mm_and_si128(_mm_add_epi16(ms, md), _mm_set1_epi32(0xFFFF));
  m = _mm_add_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(1, 0, 3, 2)));
  m = _mm_add_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(2, 3, 0, 1)));
  return (unsigned)_mm_cvtsi128_si32(m);
}

static INLINE void hbd_store_dual(__m256i sums, unsigned *out)
{
  out[0] = (unsigned)_mm_cvtsi128_si32(_mm256_castsi256_si128(sums));
  out[1] = (unsigned)_mm_cvtsi128_si32(_mm256_extracti128_si256(sums, 1));
}

static unsigned satd_4x4_subblock_hbd_avx2(const kvz_pixel *buf1, int32_t stride1,
                                           const kvz_pixel *buf2, int32_t stride2)
{
  return hbd_satd_4x4_core(_mm_sub_epi16(hbd_load_2x64(buf1, buf1 + stride1),
                                         hbd_load_2x64(buf2, buf2 + stride2)),
                           _mm_sub_epi16(hbd_load_2x64(buf1 + 2 * stride1, buf1 + 3 * stride1),
                                         hbd_load_2x64(buf2 + 2 * stride2, buf2 + 3 * stride2)));
}

/* Sum of the SATDs of two 4x4 blocks, the second at (dx, dy) from the first. */
static unsigned satd_4x4_subblock_pair_hbd_avx2(const kvz_pixel *buf1, int32_t stride1,
                                                const kvz_pixel *buf2, int32_t stride2,
                                                int dx, int dy)
{
  const int o1 = dy * stride1 + dx, o2 = dy * stride2 + dx;
  const __m256i lo = _mm256_sub_epi16(
    hbd_load_4x64(buf1, buf1 + stride1, buf1 + o1, buf1 + o1 + stride1),
    hbd_load_4x64(buf2, buf2 + stride2, buf2 + o2, buf2 + o2 + stride2));
  buf1 += 2 * stride1;
  buf2 += 2 * stride2;
  const __m256i hi = _mm256_sub_epi16(
    hbd_load_4x64(buf1, buf1 + stride1, buf1 + o1, buf1 + o1 + stride1),
    hbd_load_4x64(buf2, buf2 + stride2, buf2 + o2, buf2 + o2 + stride2));
  const __m256i s = hbd_satd_4x4_dual_core(lo, hi);
  return (unsigned)_mm_cvtsi128_si32(_mm256_castsi256_si128(s)) +
         (unsigned)_mm_cvtsi128_si32(_mm256_extracti128_si256(s, 1));
}

static void satd_4x4_subblock_quad_hbd_avx2(const kvz_pixel *preds[4], const int stride,
                                            const kvz_pixel *orig, const int orig_stride,
                                            unsigned costs[4])
{
  const __m256i o_lo = _mm256_broadcastsi128_si256(hbd_load_2x64(orig, orig + orig_stride));
  const __m256i o_hi = _mm256_broadcastsi128_si256(hbd_load_2x64(orig + 2 * orig_stride, orig + 3 * orig_stride));
  for (int k = 0; k < 4; k += 2) {
    const kvz_pixel *a = preds[k], *b = preds[k + 1];
    const __m256i lo = _mm256_sub_epi16(o_lo, hbd_load_4x64(a, a + stride, b, b + stride));
    const __m256i hi = _mm256_sub_epi16(o_hi, hbd_load_4x64(a + 2 * stride, a + 3 * stride,
                                                             b + 2 * stride, b + 3 * stride));
    hbd_store_dual(hbd_satd_4x4_dual_core(lo, hi), &costs[k]);
  }
}

static void satd_8x8_subblock_quad_hbd_avx2(const kvz_pixel **preds, const int stride,
                                            const kvz_pixel *orig, const int orig_stride,
                                            unsigned *costs)
{
  satd_8x8_subblock_dual_hbd_avx2(preds[0], preds[1], stride, orig, orig_stride, &costs[0], &costs[1]);
  satd_8x8_subblock_dual_hbd_avx2(preds[2], preds[3], stride, orig, orig_stride, &costs[2], &costs[3]);
}

static unsigned satd_4x4_hbd_avx2(const kvz_pixel *org, const kvz_pixel *cur)
{
  return hbd_satd_4x4_core(_mm_sub_epi16(_mm_loadu_si128((const __m128i *)org), _mm_loadu_si128((const __m128i *)cur)),
                           _mm_sub_epi16(_mm_loadu_si128((const __m128i *)(org + 8)),
                                         _mm_loadu_si128((const __m128i *)(cur + 8))));
}

static void satd_4x4_dual_hbd_avx2(const pred_buffer preds, const kvz_pixel *const orig,
                                   unsigned /* num_modes */, unsigned *costs_out)
{
  const __m256i o_lo = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)orig));
  const __m256i o_hi = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(orig + 8)));
  const __m256i lo = _mm256_sub_epi16(hbd_load_2x128(preds[0], preds[1]), o_lo);
  const __m256i hi = _mm256_sub_epi16(hbd_load_2x128(preds[0] + 8, preds[1] + 8), o_hi);
  hbd_store_dual(hbd_satd_4x4_dual_core(lo, hi), costs_out);
}

/* Sum of the 8x8 SATDs of an area of width x height pixels, both multiples of 8. */
static INLINE unsigned satd_8x8_blocks_hbd_avx2(int width, int height,
                                                const kvz_pixel *buf1, int stride1,
                                                const kvz_pixel *buf2, int stride2)
{
  unsigned sum = 0;
  for (int y = 0; y < height; y += 8) {
    const kvz_pixel *row1 = buf1 + y * stride1;
    const kvz_pixel *row2 = buf2 + y * stride2;
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      sum += satd_16x8_subblocks_hbd_avx2(row1 + x, stride1, row2 + x, stride2);
    }
    for (; x < width; x += 8) {
      sum += satd_8x8_subblock_hbd_avx2(row1 + x, stride1, row2 + x, stride2);
    }
  }
  return sum;
}

template <int N>
unsigned satd_nxn_hbd_avx2(const kvz_pixel *const block1, const kvz_pixel *const block2)
{
  return satd_8x8_blocks_hbd_avx2(N, N, block1, N, block2, N) >> (KVZ_BIT_DEPTH - 8);
}

template <int N>
void satd_nxn_dual_hbd_avx2(const pred_buffer preds, const kvz_pixel *const orig,
                            unsigned /* num_modes */, unsigned *costs_out)
{
  if constexpr (N == 8) {
    unsigned s0, s1;
    satd_8x8_subblock_dual_hbd_avx2(preds[0], preds[1], 8, orig, 8, &s0, &s1);
    costs_out[0] = s0 >> (KVZ_BIT_DEPTH - 8);
    costs_out[1] = s1 >> (KVZ_BIT_DEPTH - 8);
    return;
  }
  costs_out[0] = satd_8x8_blocks_hbd_avx2(N, N, preds[0], N, orig, N) >> (KVZ_BIT_DEPTH - 8);
  costs_out[1] = satd_8x8_blocks_hbd_avx2(N, N, preds[1], N, orig, N) >> (KVZ_BIT_DEPTH - 8);
}

/* Follows the SATD_ANY_SIZE macro of strategies-picture.h. */
static unsigned satd_any_size_hbd_avx2(int width, int height,
                                       const kvz_pixel *block1, int stride1,
                                       const kvz_pixel *block2, int stride2)
{
  unsigned sum = 0;
  if (width % 8 != 0) {
    // Process the first column using 4x4 blocks, two at a time.
    int y = 0;
    for (; y + 8 <= height; y += 8) {
      sum += satd_4x4_subblock_pair_hbd_avx2(&block1[y * stride1], stride1,
                                             &block2[y * stride2], stride2, 0, 4);
    }
    for (; y < height; y += 4) {
      sum += satd_4x4_subblock_hbd_avx2(&block1[y * stride1], stride1,
                                        &block2[y * stride2], stride2);
    }
    block1 += 4;
    block2 += 4;
    width -= 4;
  }
  if (height % 8 != 0) {
    // Process the first row using 4x4 blocks, two at a time.
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      sum += satd_4x4_subblock_pair_hbd_avx2(&block1[x], stride1, &block2[x], stride2, 4, 0);
    }
    for (; x < width; x += 4) {
      sum += satd_4x4_subblock_hbd_avx2(&block1[x], stride1, &block2[x], stride2);
    }
    block1 += 4 * stride1;
    block2 += 4 * stride2;
    height -= 4;
  }
  // The rest can now be processed with 8x8 blocks.
  sum += satd_8x8_blocks_hbd_avx2(width, height, block1, stride1, block2, stride2);
  return sum >> (KVZ_BIT_DEPTH - 8);
}

/* SATD of four N x N predictions against the same original block. */
template <int N>
void satd_quad_square_hbd_avx2(const kvz_pixel *const *preds, const int stride,
                               const kvz_pixel *orig, const int orig_stride,
                               unsigned *costs_out)
{
  if constexpr (N == 8) {
    unsigned c[4];
    satd_8x8_subblock_dual_hbd_avx2(preds[0], preds[1], stride, orig, orig_stride, &c[0], &c[1]);
    satd_8x8_subblock_dual_hbd_avx2(preds[2], preds[3], stride, orig, orig_stride, &c[2], &c[3]);
    for (int k = 0; k < 4; ++k) costs_out[k] = c[k] >> (KVZ_BIT_DEPTH - 8);
  } else {
    // Two horizontally adjacent blocks of a prediction at a time.
    for (int k = 0; k < 4; ++k) {
      costs_out[k] = satd_8x8_blocks_hbd_avx2(N, N, preds[k], stride, orig, orig_stride) >> (KVZ_BIT_DEPTH - 8);
    }
  }
}

/* Follows satd_any_size_quad_generic step by step, so the results are
 * identical for every input the generic accepts. */
static void satd_any_size_quad_hbd_avx2(int width, int height,
                                        const kvz_pixel **preds, const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned /* num_modes */, unsigned *costs_out, int8_t * /* valid */)
{
  // The generic quad SATD used to count the 4-high first row twice and to skip
  // the last 4 rows when the height is not a multiple of 8; satd_any_size
  // handles those sizes correctly.
  if (height % 8 != 0) {
    for (int i = 0; i < 4; i++) {
      costs_out[i] = satd_any_size_hbd_avx2(width, height, preds[i], stride, orig, orig_stride);
    }
    return;
  }
  if (width == height) {
    switch (width) {
      case  8: satd_quad_square_hbd_avx2< 8>(preds, stride, orig, orig_stride, costs_out); return;
      case 16: satd_quad_square_hbd_avx2<16>(preds, stride, orig, orig_stride, costs_out); return;
      case 32: satd_quad_square_hbd_avx2<32>(preds, stride, orig, orig_stride, costs_out); return;
      case 64: satd_quad_square_hbd_avx2<64>(preds, stride, orig, orig_stride, costs_out); return;
      default: break;
    }
  }
  unsigned sums[4] = { 0 };
  const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };
  const kvz_pixel *orig_ptr = orig;
  costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0;
  const int width_mod_8 = width % 8;
  if (width_mod_8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    // Process the first column using 4x4 blocks.
    for (int y = 0; y < height; y += 4) {
      satd_4x4_subblock_quad_hbd_avx2(pred_ptrs_tmp, stride, &orig[y * orig_stride], orig_stride, sums);
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs_tmp[blk] += 4 * stride;
        costs_out[blk] += sums[blk];
      }
    }
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] += 4;
    }
    width -= 4;
  }
  if (height % 8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    // Process the first row using 4x4 blocks.
    for (int x = 0; x < width; x += 4) {
      satd_4x4_subblock_quad_hbd_avx2(pred_ptrs_tmp, stride, &orig_ptr[x], orig_stride, sums);
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs_tmp[blk] += 4;
        costs_out[blk] += sums[blk];
      }
    }
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] += 4 * stride;
    }
    height -= 4;
  }
  // The rest can now be processed with 8x8 blocks.
  for (int y = height % 8; y < height; y += 8) {
    orig_ptr = &orig[y * orig_stride + width_mod_8];
    for (int blk = 0; blk < 4; ++blk) {
      pred_ptrs[blk] = &preds[blk][y * stride + width_mod_8];
    }
    for (int x = width_mod_8; x < width; x += 8) {
      satd_8x8_subblock_quad_hbd_avx2(pred_ptrs, stride, orig_ptr, orig_stride, sums);
      orig_ptr += 8;
      for (int blk = 0; blk < 4; ++blk) {
        pred_ptrs[blk] += 8;
        costs_out[blk] += sums[blk];
      }
    }
  }
  for (int blk = 0; blk < 4; ++blk) {
    costs_out[blk] >>= (KVZ_BIT_DEPTH - 8);
  }
}

//////////////////////////////////////////////////////////////////////////
// Temporal filter

/*
 * Quarter-pel bilinear interpolation of 16 pixels. The weights sum to 16, so
 * the result before the shift is at most 16 * 4095 + 8 = 65528: the products
 * and sums are exact in unsigned 16-bit arithmetic.
 */
static INLINE __m256i hbd_bilinear16(__m256i a, __m256i b, __m256i c, __m256i d,
                                     __m256i w00, __m256i w10, __m256i w01, __m256i w11)
{
  __m256i p = _mm256_add_epi16(_mm256_mullo_epi16(a, w00), _mm256_mullo_epi16(b, w10));
  p = _mm256_add_epi16(p, _mm256_mullo_epi16(c, w01));
  p = _mm256_add_epi16(p, _mm256_mullo_epi16(d, w11));
  return _mm256_srli_epi16(_mm256_add_epi16(p, _mm256_set1_epi16(8)), 4);
}

/*
 * Sum of squared differences of a size x size block (8 or 16) and the
 * bilinear interpolation of ref with the fraction (fx, fy), or ref itself
 * when both are 0. The squares of 13-bit differences fit 32-bit lanes, and
 * the sum wraps like the unsigned sum of the generic code.
 */
static INLINE uint32_t hbd_bilinear_ssd_block(const kvz_pixel *org, int org_stride,
                                              const kvz_pixel *ref, int ref_stride,
                                              int size, int fx, int fy)
{
  const __m256i w00 = _mm256_set1_epi16((int16_t)((4 - fx) * (4 - fy)));
  const __m256i w10 = _mm256_set1_epi16((int16_t)(fx * (4 - fy)));
  const __m256i w01 = _mm256_set1_epi16((int16_t)((4 - fx) * fy));
  const __m256i w11 = _mm256_set1_epi16((int16_t)(fx * fy));
  const bool interp = fx || fy;
  __m256i acc = _mm256_setzero_si256();
  if (size == 16) {
    for (int y = 0; y < 16; y++) {
      const kvz_pixel *r0 = &ref[y * ref_stride];
      const kvz_pixel *r1 = r0 + ref_stride;
      __m256i p = _mm256_loadu_si256((const __m256i *)r0);
      if (interp) {
        p = hbd_bilinear16(p, _mm256_loadu_si256((const __m256i *)(r0 + 1)),
                           _mm256_loadu_si256((const __m256i *)r1),
                           _mm256_loadu_si256((const __m256i *)(r1 + 1)), w00, w10, w01, w11);
      }
      const __m256i diff = _mm256_sub_epi16(_mm256_loadu_si256((const __m256i *)&org[y * org_stride]), p);
      acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
    }
  } else {
    // 8x8, two rows per iteration.
    for (int y = 0; y < 8; y += 2) {
      const kvz_pixel *r0 = &ref[y * ref_stride];
      const kvz_pixel *r1 = r0 + ref_stride;
      const kvz_pixel *r2 = r1 + ref_stride;
      __m256i p = hbd_load_2x128(r0, r1);
      if (interp) {
        p = hbd_bilinear16(p, hbd_load_2x128(r0 + 1, r1 + 1), hbd_load_2x128(r1, r2),
                           hbd_load_2x128(r1 + 1, r2 + 1), w00, w10, w01, w11);
      }
      const __m256i diff = _mm256_sub_epi16(hbd_load_2x128(&org[y * org_stride], &org[(y + 1) * org_stride]), p);
      acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
    }
  }
  return hbd_hsum_epi32(acc);
}

/* The generic bilinear SSD, for block sizes other than 8 and 16. */
static unsigned hbd_bilinear_ssd_scalar(const kvz_pixel *org, int org_stride,
                                        const kvz_pixel *ref, int ref_stride,
                                        int size, int fx, int fy)
{
  const int w00 = (4 - fx) * (4 - fy), w10 = fx * (4 - fy), w01 = (4 - fx) * fy, w11 = fx * fy;
  unsigned ssd = 0;
  for (int y = 0; y < size; y++) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    const kvz_pixel *r1 = r0 + ref_stride;
    const kvz_pixel *o = &org[y * org_stride];
    for (int x = 0; x < size; x++) {
      const int pred = (w00 * r0[x] + w10 * r0[x + 1] + w01 * r1[x] + w11 * r1[x + 1] + 8) >> 4;
      const int diff = o[x] - pred;
      ssd += diff * diff;
    }
  }
  return ssd;
}

static unsigned bilinear_ssd_hbd_avx2(const kvz_pixel *org, int org_stride,
                                      const kvz_pixel *ref, int ref_stride,
                                      int size, int fx, int fy)
{
  if (size != 8 && size != 16) {
    return hbd_bilinear_ssd_scalar(org, org_stride, ref, ref_stride, size, fx, fy);
  }
  // fx = fy = 0 gives ref itself, with or without the interpolation.
  return hbd_bilinear_ssd_block(org, org_stride, ref, ref_stride, size, fx, fy);
}

static void bilinear_ssd_multi_hbd_avx2(const kvz_pixel *org, int org_stride,
                                        const kvz_pixel *ref, int ref_stride,
                                        int size, const int *mv_x, const int *mv_y,
                                        int num, unsigned *out)
{
  if (size != 8 && size != 16) {
    kvz_bilinear_ssd_multi_generic(org, org_stride, ref, ref_stride, size, mv_x, mv_y, num, out);
    return;
  }
  for (int i = 0; i < num; i++) {
    const int fx = mv_x[i] & 3, fy = mv_y[i] & 3;
    const kvz_pixel *r = ref + (mv_y[i] >> 2) * ref_stride + (mv_x[i] >> 2);
    const uint32_t ssd = hbd_bilinear_ssd_block(org, org_stride, r, ref_stride, size, fx, fy);
    // Integer positions are scaled to 8 bits like kvz_pixels_calc_ssd.
    out[i] = (fx || fy) ? ssd : ssd >> (2 * (KVZ_BIT_DEPTH - 8));
  }
}

/*
 * The weights of tf_bilinear_pred sum to at most 64, which exceeds 16 bits
 * with 12-bit pixels, so the interpolation uses madd into 32-bit lanes.
 * Returns the 16 interpolated pixels of (a, b, c, d) in the same order.
 */
static INLINE __m256i hbd_tf_bilinear16(__m256i a, __m256i b, __m256i c, __m256i d,
                                        __m256i w0, __m256i w1, __m256i rnd, __m128i sh)
{
  __m256i lo = _mm256_add_epi32(_mm256_madd_epi16(_mm256_unpacklo_epi16(a, b), w0),
                                _mm256_madd_epi16(_mm256_unpacklo_epi16(c, d), w1));
  __m256i hi = _mm256_add_epi32(_mm256_madd_epi16(_mm256_unpackhi_epi16(a, b), w0),
                                _mm256_madd_epi16(_mm256_unpackhi_epi16(c, d), w1));
  lo = _mm256_srl_epi32(_mm256_add_epi32(lo, rnd), sh);
  hi = _mm256_srl_epi32(_mm256_add_epi32(hi, rnd), sh);
  return _mm256_packus_epi32(lo, hi);
}

static void tf_bilinear_pred_hbd_avx2(kvz_pixel *dst, int dst_stride,
                                      const kvz_pixel *ref, int ref_stride,
                                      int width, int height,
                                      int w00, int w10, int w01, int w11, int shift)
{
  if (width != 8 && width != 4) {
    kvz_tf_bilinear_pred_generic(dst, dst_stride, ref, ref_stride, width, height,
                                 w00, w10, w01, w11, shift);
    return;
  }
  // (w00, w10) and (w01, w11) pairs for the interleaved pixels.
  const __m256i w0 = _mm256_set1_epi32((int32_t)(((uint32_t)w10 << 16) | (uint16_t)w00));
  const __m256i w1 = _mm256_set1_epi32((int32_t)(((uint32_t)w11 << 16) | (uint16_t)w01));
  const __m256i rnd = _mm256_set1_epi32(1 << (shift - 1));
  const __m128i sh = _mm_cvtsi32_si128(shift);
  const int rows = width == 8 ? 2 : 4;
  int y = 0;
  for (; y + rows <= height; y += rows) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    if (width == 8) {
      const kvz_pixel *r1 = r0 + ref_stride;
      const kvz_pixel *r2 = r1 + ref_stride;
      const __m256i p = hbd_tf_bilinear16(hbd_load_2x128(r0, r1), hbd_load_2x128(r0 + 1, r1 + 1),
                                          hbd_load_2x128(r1, r2), hbd_load_2x128(r1 + 1, r2 + 1),
                                          w0, w1, rnd, sh);
      _mm_storeu_si128((__m128i *)&dst[y * dst_stride], _mm256_castsi256_si128(p));
      _mm_storeu_si128((__m128i *)&dst[(y + 1) * dst_stride], _mm256_extracti128_si256(p, 1));
    } else {
      const kvz_pixel *rr[5];
      for (int i = 0; i < 5; i++) rr[i] = r0 + i * ref_stride;
      const __m256i p = hbd_tf_bilinear16(hbd_load_4x64(rr[0], rr[1], rr[2], rr[3]),
                                          hbd_load_4x64(rr[0] + 1, rr[1] + 1, rr[2] + 1, rr[3] + 1),
                                          hbd_load_4x64(rr[1], rr[2], rr[3], rr[4]),
                                          hbd_load_4x64(rr[1] + 1, rr[2] + 1, rr[3] + 1, rr[4] + 1),
                                          w0, w1, rnd, sh);
      const __m128i lo = _mm256_castsi256_si128(p);
      const __m128i hi = _mm256_extracti128_si256(p, 1);
      _mm_storel_epi64((__m128i *)&dst[y * dst_stride], lo);
      _mm_storel_epi64((__m128i *)&dst[(y + 1) * dst_stride], _mm_unpackhi_epi64(lo, lo));
      _mm_storel_epi64((__m128i *)&dst[(y + 2) * dst_stride], hi);
      _mm_storel_epi64((__m128i *)&dst[(y + 3) * dst_stride], _mm_unpackhi_epi64(hi, hi));
    }
  }
  if (y < height) {
    kvz_tf_bilinear_pred_generic(&dst[y * dst_stride], dst_stride, &ref[y * ref_stride], ref_stride,
                                 width, height - y, w00, w10, w01, w11, shift);
  }
}

/* Sum of the four signed 32-bit lanes, in 64 bits. */
static INLINE int64_t hbd_hsum_epi32_i64(__m128i v)
{
  int32_t t[4];
  _mm_storeu_si128((__m128i *)t, v);
  return (int64_t)t[0] + t[1] + t[2] + t[3];
}

/*
 * 8x8 blocks only. Per 32-bit lane the largest sum, that of the squared
 * gradients, is at most 7 rows * 4 * 8190^2 < 2^31 at 12 bits.
 */
static void tf_block_stats_hbd_avx2(const kvz_pixel *org, int org_stride,
                                    const kvz_pixel *pred, int pred_stride,
                                    int width, int height, int64_t out[5])
{
  if (width != 8 || height != 8) {
    kvz_tf_block_stats_generic(org, org_stride, pred, pred_stride, width, height, out);
    return;
  }
  const __m128i ones = _mm_set1_epi16(1);
  // Lanes 0..6: the pixels that have a right neighbour.
  const __m128i mask7 = _mm_setr_epi16(-1, -1, -1, -1, -1, -1, -1, 0);
  __m128i o[8], d[8];
  for (int y = 0; y < 8; y++) {
    o[y] = _mm_loadu_si128((const __m128i *)&org[y * org_stride]);
    d[y] = _mm_sub_epi16(o[y], _mm_loadu_si128((const __m128i *)&pred[y * pred_stride]));
  }
  __m128i sum = _mm_setzero_si128(), sum2 = _mm_setzero_si128(), sse = _mm_setzero_si128();
  __m128i var_d = _mm_setzero_si128(), diff = _mm_setzero_si128();
  for (int y = 0; y < 8; y++) {
    sum = _mm_add_epi32(sum, _mm_madd_epi16(o[y], ones));
    sum2 = _mm_add_epi32(sum2, _mm_madd_epi16(o[y], o[y]));
    sse = _mm_add_epi32(sse, _mm_madd_epi16(d[y], d[y]));
    if (y < 7) {
      const __m128i dm = _mm_and_si128(d[y], mask7);
      var_d = _mm_add_epi32(var_d, _mm_madd_epi16(dm, dm));
      const __m128i dh = _mm_and_si128(_mm_sub_epi16(_mm_srli_si128(d[y], 2), d[y]), mask7);
      const __m128i dv = _mm_and_si128(_mm_sub_epi16(d[y + 1], d[y]), mask7);
      diff = _mm_add_epi32(diff, _mm_add_epi32(_mm_madd_epi16(dh, dh), _mm_madd_epi16(dv, dv)));
    }
  }
  out[0] = hbd_hsum_epi32_i64(sum);
  out[1] = hbd_hsum_epi32_i64(sum2);
  out[2] = hbd_hsum_epi32_i64(sse);
  out[3] = hbd_hsum_epi32_i64(var_d);
  out[4] = hbd_hsum_epi32_i64(diff);
}

/* Filter 8 pixels given as 8 16-bit lanes of src and of each prediction. The
 * operations and their order are those of kvz_tf_filter_block_generic. */
static INLINE __m128i hbd_tf_filter_8px(__m128i src8, const __m128i *pred8, int num_refs,
                                        const float *weights, const float *const *luts,
                                        __m256i max_val)
{
  const __m256i o = _mm256_cvtepu16_epi32(src8);
  __m256 sum = _mm256_cvtepi32_ps(o);
  __m256 wsum = _mm256_set1_ps(1.0f);
  for (int r = 0; r < num_refs; r++) {
    const __m256i p = _mm256_cvtepu16_epi32(pred8[r]);
    const __m256 lut = _mm256_i32gather_ps(luts[r], _mm256_sub_epi32(p, o), 4);
    const __m256 wgt = _mm256_mul_ps(_mm256_set1_ps(weights[r]), lut);
    sum = _mm256_add_ps(sum, _mm256_mul_ps(wgt, _mm256_cvtepi32_ps(p)));
    wsum = _mm256_add_ps(wsum, wgt);
  }
  __m256i res = _mm256_cvttps_epi32(_mm256_add_ps(_mm256_div_ps(sum, wsum), _mm256_set1_ps(0.5f)));
  res = _mm256_min_epi32(_mm256_max_epi32(res, _mm256_setzero_si256()), max_val);
  return _mm_packus_epi32(_mm256_castsi256_si128(res), _mm256_extracti128_si256(res, 1));
}

static void tf_filter_block_hbd_avx2(kvz_pixel *dst, int dst_stride,
                                     const kvz_pixel *src, int src_stride,
                                     const kvz_pixel *const *preds, int pred_stride,
                                     int num_refs, const float *weights,
                                     const float *const *luts,
                                     int width, int height, int max_val)
{
  if (num_refs > 8 || !((width % 8 == 0) || (width == 4 && height % 2 == 0))) {
    kvz_tf_filter_block_generic(dst, dst_stride, src, src_stride, preds, pred_stride,
                                num_refs, weights, luts, width, height, max_val);
    return;
  }
  const __m256i maxv = _mm256_set1_epi32(max_val);
  __m128i p8[8];
  if (width == 4) {
    for (int y = 0; y < height; y += 2) {
      for (int r = 0; r < num_refs; r++) {
        p8[r] = hbd_load_2x64(&preds[r][y * pred_stride], &preds[r][(y + 1) * pred_stride]);
      }
      const __m128i res = hbd_tf_filter_8px(hbd_load_2x64(&src[y * src_stride], &src[(y + 1) * src_stride]),
                                            p8, num_refs, weights, luts, maxv);
      _mm_storel_epi64((__m128i *)&dst[y * dst_stride], res);
      _mm_storel_epi64((__m128i *)&dst[(y + 1) * dst_stride], _mm_unpackhi_epi64(res, res));
    }
    return;
  }
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x += 8) {
      for (int r = 0; r < num_refs; r++) {
        p8[r] = _mm_loadu_si128((const __m128i *)&preds[r][y * pred_stride + x]);
      }
      const __m128i res = hbd_tf_filter_8px(_mm_loadu_si128((const __m128i *)&src[y * src_stride + x]),
                                            p8, num_refs, weights, luts, maxv);
      _mm_storeu_si128((__m128i *)&dst[y * dst_stride + x], res);
    }
  }
}

} // namespace
#endif // COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8

#if COMPILE_INTEL_AVX2
#if KVZ_BIT_DEPTH == 8

/**
 * \brief Calculate Sum of Absolute Differences (SAD)
 *
 * Calculate Sum of Absolute Differences (SAD) between two rectangular regions
 * located in arbitrary points in the picture.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param stride  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
extern "C" uint32_t kvz_reg_sad_avx2(const uint8_t * const data1, const uint8_t * const data2,
                          const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  if (width == 0)
    return 0;
  if (width < 4) {
    // Too narrow for the vector paths. A plain loop is faster than their setup.
    uint32_t sum = 0;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        int d = (int)data1[y * stride1 + x] - (int)data2[y * stride2 + x];
        sum += (uint32_t)(d < 0 ? -d : d);
      }
    }
    return sum;
  }
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
  if (width == 32)
    return reg_sad_w32(data1, data2, height, stride1, stride2);
  if (width == 64)
    return reg_sad_w64(data1, data2, height, stride1, stride2);
  else
    return reg_sad_arbitrary(data1, data2, width, height, stride1, stride2);
}

namespace {

/**
* \brief Calculate SAD for 8x8 bytes in continuous memory.
*/
inline __m256i inline_8bit_sad_8x8_avx2(const __m256i *const a, const __m256i *const b)
{
  // The SAD inputs are caller-provided buffers (e.g. lcu->rec, only 16-byte
  // aligned), so use unaligned loads.
  __m256i sum0, sum1;
  sum0 = _mm256_sad_epu8(_mm256_loadu_si256(a + 0), _mm256_loadu_si256(b + 0));
  sum1 = _mm256_sad_epu8(_mm256_loadu_si256(a + 1), _mm256_loadu_si256(b + 1));

  return _mm256_add_epi32(sum0, sum1);
}


/**
* \brief Calculate SAD for 16x16 bytes in continuous memory.
*/
inline __m256i inline_8bit_sad_16x16_avx2(const __m256i *const a, const __m256i *const b)
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
inline uint32_t m256i_horizontal_sum(const __m256i sum)
{
  // Add the high 128 bits to low 128 bits.
  __m128i mm128_result = _mm_add_epi32(_mm256_castsi256_si128(sum), _mm256_extractf128_si256(sum, 1));
  // Add the high 64 bits  to low 64 bits.
  uint32_t result[4];
  _mm_storeu_si128((__m128i*)result, mm128_result);
  return result[0] + result[2];
}


unsigned sad_8bit_8x8_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_8x8_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


unsigned sad_8bit_16x16_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  const __m256i *const a = (const __m256i *)buf1;
  const __m256i *const b = (const __m256i *)buf2;
  __m256i sum = inline_8bit_sad_16x16_avx2(a, b);

  return m256i_horizontal_sum(sum);
}


unsigned sad_8bit_32x32_avx2(const uint8_t *buf1, const uint8_t *buf2)
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


unsigned sad_8bit_64x64_avx2(const uint8_t * buf1, const uint8_t * buf2)
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

/*
 * 4x4 SATD cores on 16-byte vectors holding a 4x4 block row by row. The
 * horizontal and vertical Hadamard butterflies are hadd/hsub pairs; the sum of
 * the 16 absolute values (at most 65280) is read as an unsigned 16-bit lane.
 */
inline unsigned satd_4x4_core_avx2(__m128i org_lo, __m128i org_hi, __m128i cur_lo, __m128i cur_hi)
{
  const __m128i diff_lo = _mm_sub_epi16(_mm_cvtepu8_epi16(cur_lo), _mm_cvtepu8_epi16(org_lo));
  const __m128i diff_hi = _mm_sub_epi16(_mm_cvtepu8_epi16(cur_hi), _mm_cvtepu8_epi16(org_hi));

  //Hor
  __m128i row0 = _mm_hadd_epi16(diff_lo, diff_hi);
  __m128i row1 = _mm_hsub_epi16(diff_lo, diff_hi);

  __m128i row2 = _mm_hadd_epi16(row0, row1);
  __m128i row3 = _mm_hsub_epi16(row0, row1);

  //Ver
  row0 = _mm_hadd_epi16(row2, row3);
  row1 = _mm_hsub_epi16(row2, row3);

  row2 = _mm_hadd_epi16(row0, row1);
  row3 = _mm_hsub_epi16(row0, row1);

  //Abs and sum
  row2 = _mm_abs_epi16(row2);
  row3 = _mm_abs_epi16(row3);

  row3 = _mm_add_epi16(row2, row3);

  row3 = _mm_add_epi16(row3, _mm_shuffle_epi32(row3, _MM_SHUFFLE(1, 0, 3, 2) ));
  row3 = _mm_add_epi16(row3, _mm_shuffle_epi32(row3, _MM_SHUFFLE(0, 1, 0, 1) ));
  row3 = _mm_add_epi16(row3, _mm_shufflelo_epi16(row3, _MM_SHUFFLE(0, 1, 0, 1) ));

  unsigned sum = _mm_extract_epi16(row3, 0);
  return (sum + 1) >> 1;
}

inline void satd_4x4_dual_core_avx2(__m128i orig_lo, __m128i orig_hi, __m128i p0_lo, __m128i p0_hi,
                                           __m128i p1_lo, __m128i p1_hi, unsigned *satds_out)
{
  const __m256i original_lo = _mm256_broadcastsi128_si256(_mm_cvtepu8_epi16(orig_lo));
  const __m256i original_hi = _mm256_broadcastsi128_si256(_mm_cvtepu8_epi16(orig_hi));
  __m256i pred_lo = _mm256_cvtepu8_epi16(p0_lo);
  pred_lo = _mm256_inserti128_si256(pred_lo, _mm_cvtepu8_epi16(p1_lo), 1);
  __m256i pred_hi = _mm256_cvtepu8_epi16(p0_hi);
  pred_hi = _mm256_inserti128_si256(pred_hi, _mm_cvtepu8_epi16(p1_hi), 1);

  const __m256i diff_lo = _mm256_sub_epi16(pred_lo, original_lo);
  const __m256i diff_hi = _mm256_sub_epi16(pred_hi, original_hi);

  //Hor
  __m256i row0 = _mm256_hadd_epi16(diff_lo, diff_hi);
  __m256i row1 = _mm256_hsub_epi16(diff_lo, diff_hi);

  __m256i row2 = _mm256_hadd_epi16(row0, row1);
  __m256i row3 = _mm256_hsub_epi16(row0, row1);

  //Ver
  row0 = _mm256_hadd_epi16(row2, row3);
  row1 = _mm256_hsub_epi16(row2, row3);

  row2 = _mm256_hadd_epi16(row0, row1);
  row3 = _mm256_hsub_epi16(row0, row1);

  //Abs and sum
  row2 = _mm256_abs_epi16(row2);
  row3 = _mm256_abs_epi16(row3);

  row3 = _mm256_add_epi16(row2, row3);

  row3 = _mm256_add_epi16(row3, _mm256_shuffle_epi32(row3, _MM_SHUFFLE(1, 0, 3, 2) ));
  row3 = _mm256_add_epi16(row3, _mm256_shuffle_epi32(row3, _MM_SHUFFLE(0, 1, 0, 1) ));
  row3 = _mm256_add_epi16(row3, _mm256_shufflelo_epi16(row3, _MM_SHUFFLE(0, 1, 0, 1) ));

  unsigned sum1 = _mm_extract_epi16(_mm256_castsi256_si128(row3), 0);
  unsigned sum2 = _mm_extract_epi16(_mm256_extracti128_si256(row3, 1), 0);
  satds_out[0] = (sum1 + 1) >> 1;
  satds_out[1] = (sum2 + 1) >> 1;
}

// Packs rows y and y + 1 of a 4-wide strided block into the low 8 bytes.
inline __m128i load_2rows4_avx2(const uint8_t *buf, int32_t stride)
{
  uint32_t r0, r1;
  memcpy(&r0, buf, 4);
  memcpy(&r1, buf + stride, 4);
  return _mm_cvtsi64_si128((long long)((uint64_t)r0 | ((uint64_t)r1 << 32)));
}
#define LOAD_4X4_HALVES(buf, stride) load_2rows4_avx2((buf), (stride)), load_2rows4_avx2((buf) + 2 * (stride), (stride))

unsigned satd_4x4_8bit_avx2(const uint8_t *org, const uint8_t *cur)
{
  return satd_4x4_core_avx2(_mm_loadl_epi64((const __m128i *)org), _mm_loadl_epi64((const __m128i *)(org + 8)),
                            _mm_loadl_epi64((const __m128i *)cur), _mm_loadl_epi64((const __m128i *)(cur + 8)));
}

void satd_8bit_4x4_dual_avx2(
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out) 
{
  satd_4x4_dual_core_avx2(_mm_loadl_epi64((const __m128i *)orig), _mm_loadl_epi64((const __m128i *)(orig + 8)),
                          _mm_loadl_epi64((const __m128i *)preds[0]), _mm_loadl_epi64((const __m128i *)(preds[0] + 8)),
                          _mm_loadl_epi64((const __m128i *)preds[1]), _mm_loadl_epi64((const __m128i *)(preds[1] + 8)),
                          satds_out);
}

inline void hor_transform_row_avx2(__m128i* row){
  
  __m128i mask_pos = _mm_set1_epi16(1);
  __m128i mask_neg = _mm_set1_epi16(-1);
  __m128i sign_mask = _mm_unpacklo_epi64(mask_pos, mask_neg);
  __m128i temp = _mm_shuffle_epi32(*row, _MM_SHUFFLE(1, 0, 3, 2));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);

  sign_mask = _mm_unpacklo_epi32(mask_pos, mask_neg);
  temp = _mm_shuffle_epi32(*row, _MM_SHUFFLE(2, 3, 0, 1));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);

  sign_mask = _mm_unpacklo_epi16(mask_pos, mask_neg);
  temp = _mm_shufflelo_epi16(*row, _MM_SHUFFLE(2,3,0,1));
  temp = _mm_shufflehi_epi16(temp, _MM_SHUFFLE(2,3,0,1));
  *row = _mm_sign_epi16(*row, sign_mask);
  *row = _mm_add_epi16(*row, temp);
}

inline void hor_transform_row_dual_avx2(__m256i* row){
  
  __m256i mask_pos = _mm256_set1_epi16(1);
  __m256i mask_neg = _mm256_set1_epi16(-1);
  __m256i sign_mask = _mm256_unpacklo_epi64(mask_pos, mask_neg);
  __m256i temp = _mm256_shuffle_epi32(*row, _MM_SHUFFLE(1, 0, 3, 2));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);

  sign_mask = _mm256_unpacklo_epi32(mask_pos, mask_neg);
  temp = _mm256_shuffle_epi32(*row, _MM_SHUFFLE(2, 3, 0, 1));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);

  sign_mask = _mm256_unpacklo_epi16(mask_pos, mask_neg);
  temp = _mm256_shufflelo_epi16(*row, _MM_SHUFFLE(2,3,0,1));
  temp = _mm256_shufflehi_epi16(temp, _MM_SHUFFLE(2,3,0,1));
  *row = _mm256_sign_epi16(*row, sign_mask);
  *row = _mm256_add_epi16(*row, temp);
}

inline void add_sub_avx2(__m128i *out, __m128i *in, unsigned out_idx0, unsigned out_idx1, unsigned in_idx0, unsigned in_idx1)
{
  out[out_idx0] = _mm_add_epi16(in[in_idx0], in[in_idx1]);
  out[out_idx1] = _mm_sub_epi16(in[in_idx0], in[in_idx1]);
}

inline void add_sub_dual_avx2(__m256i *out, __m256i *in, unsigned out_idx0, unsigned out_idx1, unsigned in_idx0, unsigned in_idx1)
{
  out[out_idx0] = _mm256_add_epi16(in[in_idx0], in[in_idx1]);
  out[out_idx1] = _mm256_sub_epi16(in[in_idx0], in[in_idx1]);
}


/*
 * Last vertical butterfly stage fused with the absolute sum. The stage pairs
 * rows (0, 4), (1, 5), (2, 6) and (3, 7), and |a + b| + |a - b| ==
 * 2 * max(|a|, |b|) for integers, so the stage, half of the absolute values and
 * half of the accumulation are replaced by maxima. Returns half of the
 * absolute sum. The terms are at most 8160 for 8-bit input, so the sum of the
 * four maxima (at most 32640) fits a 16-bit lane.
 */
inline unsigned ver_transform_absmax_sum_avx2(__m128i (*rows)[8])
{
  __m128i temp0[8], temp1[8];
  add_sub_avx2(temp0, (*rows), 0, 1, 0, 1);
  add_sub_avx2(temp0, (*rows), 2, 3, 2, 3);
  add_sub_avx2(temp0, (*rows), 4, 5, 4, 5);
  add_sub_avx2(temp0, (*rows), 6, 7, 6, 7);

  add_sub_avx2(temp1, temp0, 0, 1, 0, 2);
  add_sub_avx2(temp1, temp0, 2, 3, 1, 3);
  add_sub_avx2(temp1, temp0, 4, 5, 4, 6);
  add_sub_avx2(temp1, temp0, 6, 7, 5, 7);

  const __m128i m0 = _mm_max_epi16(_mm_abs_epi16(temp1[0]), _mm_abs_epi16(temp1[4]));
  const __m128i m1 = _mm_max_epi16(_mm_abs_epi16(temp1[1]), _mm_abs_epi16(temp1[5]));
  const __m128i m2 = _mm_max_epi16(_mm_abs_epi16(temp1[2]), _mm_abs_epi16(temp1[6]));
  const __m128i m3 = _mm_max_epi16(_mm_abs_epi16(temp1[3]), _mm_abs_epi16(temp1[7]));
  const __m128i sum16 = _mm_add_epi16(_mm_add_epi16(m0, m1), _mm_add_epi16(m2, m3));
  __m128i sum = _mm_madd_epi16(sum16, _mm_set1_epi16(1));
  sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
  sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));
  return (unsigned)_mm_cvtsi128_si32(sum);
}

inline void ver_transform_absmax_sum_dual_avx2(__m256i (*rows)[8], unsigned *half0, unsigned *half1)
{
  __m256i temp0[8], temp1[8];
  add_sub_dual_avx2(temp0, (*rows), 0, 1, 0, 1);
  add_sub_dual_avx2(temp0, (*rows), 2, 3, 2, 3);
  add_sub_dual_avx2(temp0, (*rows), 4, 5, 4, 5);
  add_sub_dual_avx2(temp0, (*rows), 6, 7, 6, 7);

  add_sub_dual_avx2(temp1, temp0, 0, 1, 0, 2);
  add_sub_dual_avx2(temp1, temp0, 2, 3, 1, 3);
  add_sub_dual_avx2(temp1, temp0, 4, 5, 4, 6);
  add_sub_dual_avx2(temp1, temp0, 6, 7, 5, 7);

  const __m256i m0 = _mm256_max_epi16(_mm256_abs_epi16(temp1[0]), _mm256_abs_epi16(temp1[4]));
  const __m256i m1 = _mm256_max_epi16(_mm256_abs_epi16(temp1[1]), _mm256_abs_epi16(temp1[5]));
  const __m256i m2 = _mm256_max_epi16(_mm256_abs_epi16(temp1[2]), _mm256_abs_epi16(temp1[6]));
  const __m256i m3 = _mm256_max_epi16(_mm256_abs_epi16(temp1[3]), _mm256_abs_epi16(temp1[7]));
  const __m256i sum16 = _mm256_add_epi16(_mm256_add_epi16(m0, m1), _mm256_add_epi16(m2, m3));
  __m256i sum = _mm256_madd_epi16(sum16, _mm256_set1_epi16(1));
  sum = _mm256_add_epi32(sum, _mm256_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
  sum = _mm256_add_epi32(sum, _mm256_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));
  *half0 = (unsigned)_mm_cvtsi128_si32(_mm256_castsi256_si128(sum));
  *half1 = (unsigned)_mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 1));
}

inline __m128i diff_row_avx2(const uint8_t *buf1, const uint8_t *buf2)
{
  __m128i buf1_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)buf1));
  __m128i buf2_row = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i*)buf2));
  return _mm_sub_epi16(buf1_row, buf2_row);
}

inline __m256i diff_row_dual_avx2(const uint8_t *buf1, const uint8_t *buf2, const uint8_t *orig)
{
  __m128i temp1 = _mm_loadl_epi64((__m128i*)buf1);
  __m128i temp2 = _mm_loadl_epi64((__m128i*)buf2);
  __m128i temp3 = _mm_loadl_epi64((__m128i*)orig);
  __m256i buf1_row = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(temp1, temp2));
  __m256i buf2_row = _mm256_cvtepu8_epi16(_mm_broadcastq_epi64(temp3));

  return _mm256_sub_epi16(buf1_row, buf2_row);
}

inline void diff_blocks_avx2(__m128i (*row_diff)[8],
                                                           const uint8_t * buf1, unsigned stride1,
                                                           const uint8_t * orig, unsigned stride_orig)
{
  (*row_diff)[0] = diff_row_avx2(buf1 + 0 * stride1, orig + 0 * stride_orig);
  (*row_diff)[1] = diff_row_avx2(buf1 + 1 * stride1, orig + 1 * stride_orig);
  (*row_diff)[2] = diff_row_avx2(buf1 + 2 * stride1, orig + 2 * stride_orig);
  (*row_diff)[3] = diff_row_avx2(buf1 + 3 * stride1, orig + 3 * stride_orig);
  (*row_diff)[4] = diff_row_avx2(buf1 + 4 * stride1, orig + 4 * stride_orig);
  (*row_diff)[5] = diff_row_avx2(buf1 + 5 * stride1, orig + 5 * stride_orig);
  (*row_diff)[6] = diff_row_avx2(buf1 + 6 * stride1, orig + 6 * stride_orig);
  (*row_diff)[7] = diff_row_avx2(buf1 + 7 * stride1, orig + 7 * stride_orig);

}

inline void diff_blocks_dual_avx2(__m256i (*row_diff)[8],
                                                           const uint8_t * buf1, unsigned stride1,
                                                           const uint8_t * buf2, unsigned stride2,
                                                           const uint8_t * orig, unsigned stride_orig)
{
  (*row_diff)[0] = diff_row_dual_avx2(buf1 + 0 * stride1, buf2 + 0 * stride2, orig + 0 * stride_orig);
  (*row_diff)[1] = diff_row_dual_avx2(buf1 + 1 * stride1, buf2 + 1 * stride2, orig + 1 * stride_orig);
  (*row_diff)[2] = diff_row_dual_avx2(buf1 + 2 * stride1, buf2 + 2 * stride2, orig + 2 * stride_orig);
  (*row_diff)[3] = diff_row_dual_avx2(buf1 + 3 * stride1, buf2 + 3 * stride2, orig + 3 * stride_orig);
  (*row_diff)[4] = diff_row_dual_avx2(buf1 + 4 * stride1, buf2 + 4 * stride2, orig + 4 * stride_orig);
  (*row_diff)[5] = diff_row_dual_avx2(buf1 + 5 * stride1, buf2 + 5 * stride2, orig + 5 * stride_orig);
  (*row_diff)[6] = diff_row_dual_avx2(buf1 + 6 * stride1, buf2 + 6 * stride2, orig + 6 * stride_orig);
  (*row_diff)[7] = diff_row_dual_avx2(buf1 + 7 * stride1, buf2 + 7 * stride2, orig + 7 * stride_orig);

}

inline void hor_transform_block_avx2(__m128i (*row_diff)[8])
{
  hor_transform_row_avx2((*row_diff) + 0);
  hor_transform_row_avx2((*row_diff) + 1);
  hor_transform_row_avx2((*row_diff) + 2);
  hor_transform_row_avx2((*row_diff) + 3);
  hor_transform_row_avx2((*row_diff) + 4);
  hor_transform_row_avx2((*row_diff) + 5);
  hor_transform_row_avx2((*row_diff) + 6);
  hor_transform_row_avx2((*row_diff) + 7);
}

inline void hor_transform_block_dual_avx2(__m256i (*row_diff)[8])
{
  hor_transform_row_dual_avx2((*row_diff) + 0);
  hor_transform_row_dual_avx2((*row_diff) + 1);
  hor_transform_row_dual_avx2((*row_diff) + 2);
  hor_transform_row_dual_avx2((*row_diff) + 3);
  hor_transform_row_dual_avx2((*row_diff) + 4);
  hor_transform_row_dual_avx2((*row_diff) + 5);
  hor_transform_row_dual_avx2((*row_diff) + 6);
  hor_transform_row_dual_avx2((*row_diff) + 7);
}

void kvz_satd_8bit_8x8_general_dual_avx2(const uint8_t * buf1, unsigned stride1,
                                                const uint8_t * buf2, unsigned stride2,
                                                const uint8_t * orig, unsigned stride_orig,
                                                unsigned *sum0, unsigned *sum1)
{
  __m256i temp[8];

  diff_blocks_dual_avx2(&temp, buf1, stride1, buf2, stride2, orig, stride_orig);
  hor_transform_block_dual_avx2(&temp);
  unsigned half0, half1;
  ver_transform_absmax_sum_dual_avx2(&temp, &half0, &half1);

  // (sum + 2) >> 2 with sum = 2 * half.
  *sum0 = (2 * half0 + 2) >> 2;
  *sum1 = (2 * half1 + 2) >> 2;
}

/**
* \brief  Calculate SATD between two 4x4 blocks inside bigger arrays.
*/
unsigned kvz_satd_4x4_subblock_8bit_avx2(const uint8_t * buf1,
                                                const int32_t     stride1,
                                                const uint8_t * buf2,
                                                const int32_t     stride2)
{
  return satd_4x4_core_avx2(LOAD_4X4_HALVES(buf2, stride2), LOAD_4X4_HALVES(buf1, stride1));
}

void kvz_satd_4x4_subblock_quad_avx2(const uint8_t *preds[4],
                                       const int stride,
                                       const uint8_t *orig,
                                       const int orig_stride,
                                       unsigned costs[4])
{
  const __m128i o_lo = load_2rows4_avx2(orig, orig_stride);
  const __m128i o_hi = load_2rows4_avx2(orig + 2 * orig_stride, orig_stride);
  satd_4x4_dual_core_avx2(o_lo, o_hi, LOAD_4X4_HALVES(preds[0], stride), LOAD_4X4_HALVES(preds[1], stride), &costs[0]);
  satd_4x4_dual_core_avx2(o_lo, o_hi, LOAD_4X4_HALVES(preds[2], stride), LOAD_4X4_HALVES(preds[3], stride), &costs[2]);
}

unsigned satd_8x8_subblock_8bit_avx2(const uint8_t * buf1, unsigned stride1, const uint8_t * buf2, unsigned stride2)
{
  __m128i temp[8];

  diff_blocks_avx2(&temp, buf1, stride1, buf2, stride2);
  hor_transform_block_avx2(&temp);

  // (sum + 2) >> 2 with sum = 2 * half.
  return (2 * ver_transform_absmax_sum_avx2(&temp) + 2) >> 2;
}

void satd_8x8_subblock_quad_avx2(const uint8_t **preds,
  const int stride,
  const uint8_t *orig,
  const int orig_stride,
  unsigned *costs)
{
  kvz_satd_8bit_8x8_general_dual_avx2(preds[0], stride, preds[1], stride, orig, orig_stride, &costs[0], &costs[1]);
  kvz_satd_8bit_8x8_general_dual_avx2(preds[2], stride, preds[3], stride, orig, orig_stride, &costs[2], &costs[3]);
}

SATD_NxN(8bit_avx2,  8)
SATD_NxN(8bit_avx2, 16)
SATD_NxN(8bit_avx2, 32)
SATD_NxN(8bit_avx2, 64)
SATD_ANY_SIZE(8bit_avx2)

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NXN_DUAL_AVX2(n) \
void satd_8bit_ ## n ## x ## n ## _dual_avx2( \
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out)  \
{ \
  unsigned x, y; \
  satds_out[0] = 0; \
  satds_out[1] = 0; \
  unsigned sum1 = 0; \
  unsigned sum2 = 0; \
  for (y = 0; y < (n); y += 8) { \
  unsigned row = y * (n); \
  for (x = 0; x < (n); x += 8) { \
  /* Plain pointers: the 32x32 rows of pred_buffer are too small for 64x64. */ \
  kvz_satd_8bit_8x8_general_dual_avx2((const kvz_pixel *)preds[0] + row + x, (n), (const kvz_pixel *)preds[1] + row + x, (n), &orig[row + x], (n), &sum1, &sum2); \
  satds_out[0] += sum1; \
  satds_out[1] += sum2; \
    } \
    } \
  satds_out[0] >>= (KVZ_BIT_DEPTH-8); \
  satds_out[1] >>= (KVZ_BIT_DEPTH-8); \
}

void satd_8bit_8x8_dual_avx2(
  const pred_buffer preds, const uint8_t * const orig, unsigned num_modes, unsigned *satds_out) 
{ 
  unsigned x, y; 
  satds_out[0] = 0;
  satds_out[1] = 0;
  unsigned sum1 = 0;
  unsigned sum2 = 0;
  for (y = 0; y < (8); y += 8) { 
  unsigned row = y * (8); 
  for (x = 0; x < (8); x += 8) { 
  kvz_satd_8bit_8x8_general_dual_avx2(&preds[0][row + x], (8), &preds[1][row + x], (8), &orig[row + x], (8), &sum1, &sum2); 
  satds_out[0] += sum1;
  satds_out[1] += sum2;
      } 
      } 
  satds_out[0] >>= (KVZ_BIT_DEPTH-8);
  satds_out[1] >>= (KVZ_BIT_DEPTH-8);
}

//SATD_NXN_DUAL_AVX2(8) //Use the non-macro version
SATD_NXN_DUAL_AVX2(16)
SATD_NXN_DUAL_AVX2(32)
SATD_NXN_DUAL_AVX2(64)

/*
 * SATD of four N x N predictions against the same original block, for the
 * square sizes. N is a template parameter, so the block loops are constant
 * and the 8x8 dual kernel is inlined into them.
 */
template <int N>
void satd_quad_square_avx2(const kvz_pixel *const *preds, const int stride,
                           const kvz_pixel *orig, const int orig_stride,
                           unsigned *costs_out)
{
  unsigned c0 = 0, c1 = 0, c2 = 0, c3 = 0;
  for (int y = 0; y < N; y += 8) {
    for (int x = 0; x < N; x += 8) {
      const int p = y * stride + x;
      const kvz_pixel *o = orig + y * orig_stride + x;
      unsigned s0, s1, s2, s3;
      kvz_satd_8bit_8x8_general_dual_avx2(preds[0] + p, stride, preds[1] + p, stride, o, orig_stride, &s0, &s1);
      kvz_satd_8bit_8x8_general_dual_avx2(preds[2] + p, stride, preds[3] + p, stride, o, orig_stride, &s2, &s3);
      c0 += s0; c1 += s1; c2 += s2; c3 += s3;
    }
  }
  costs_out[0] = c0 >> (KVZ_BIT_DEPTH - 8);
  costs_out[1] = c1 >> (KVZ_BIT_DEPTH - 8);
  costs_out[2] = c2 >> (KVZ_BIT_DEPTH - 8);
  costs_out[3] = c3 >> (KVZ_BIT_DEPTH - 8);
}

#define SATD_ANY_SIZE_MULTI_AVX2(suffix, num_parallel_blocks) \
  static cost_pixel_any_size_multi_func satd_any_size_## suffix; \
  static void satd_any_size_ ## suffix ( \
      int width, int height, \
      const kvz_pixel **preds, \
      const int stride, \
      const kvz_pixel *orig, \
      const int orig_stride, \
      unsigned num_modes, \
      unsigned *costs_out, \
      int8_t *valid) \
  { \
    /* Heights that are not multiples of 8 as satd_any_size (see the HBD version). */ \
    if (height % 8 != 0) { \
      for (int i = 0; i < num_parallel_blocks; i++) { \
        costs_out[i] = satd_any_size_8bit_avx2(width, height, preds[i], stride, orig, orig_stride); \
      } \
      return; \
    } \
    if (num_parallel_blocks == 4 && width == height) { \
      switch (width) { \
        case  8: satd_quad_square_avx2< 8>(preds, stride, orig, orig_stride, costs_out); return; \
        case 16: satd_quad_square_avx2<16>(preds, stride, orig, orig_stride, costs_out); return; \
        case 32: satd_quad_square_avx2<32>(preds, stride, orig, orig_stride, costs_out); return; \
        case 64: satd_quad_square_avx2<64>(preds, stride, orig, orig_stride, costs_out); return; \
        default: break; \
      } \
    } \
    unsigned sums[num_parallel_blocks] = { 0 }; \
    const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };\
    const kvz_pixel *orig_ptr = orig; \
    costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0; \
    const int width_mod_8 = width % 8; \
    if (width_mod_8 != 0) { \
      const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };\
      /* Process the first column using 4x4 blocks. */ \
      for (int y = 0; y < height; y += 4) { \
        kvz_satd_4x4_subblock_ ## suffix(pred_ptrs_tmp, stride, &orig[y*orig_stride], orig_stride, sums); \
        for(int blk = 0; blk < num_parallel_blocks; ++blk){ pred_ptrs_tmp[blk] += 4*stride; costs_out[blk] += sums[blk]; }\
      } \
      for(int blk = 0; blk < num_parallel_blocks; ++blk){\
        pred_ptrs[blk] += 4; \
      }\
      width -= 4; \
    } \
    if (height % 8 != 0) { \
      const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };\
      /* Process the first row using 4x4 blocks. */ \
      for (int x = 0; x < width; x += 4 ) { \
        kvz_satd_4x4_subblock_ ## suffix(pred_ptrs_tmp, stride, &orig_ptr[x], orig_stride, sums); \
        for(int blk = 0; blk < num_parallel_blocks; ++blk){ pred_ptrs_tmp[blk] += 4;  costs_out[blk] += sums[blk]; }\
      } \
      for(int blk = 0; blk < num_parallel_blocks; ++blk){\
        pred_ptrs[blk] += 4 * stride; \
      }\
      height -= 4; \
    } \
    /* The rest can now be processed with 8x8 blocks. */ \
    for (int y = height % 8; y < height; y += 8) { \
      orig_ptr = &orig[y * orig_stride + width_mod_8]; \
      pred_ptrs[0] = &preds[0][y * stride + width_mod_8]; \
      pred_ptrs[1] = &preds[1][y * stride + width_mod_8]; \
      pred_ptrs[2] = &preds[2][y * stride + width_mod_8]; \
      pred_ptrs[3] = &preds[3][y * stride + width_mod_8]; \
      for (int x = width_mod_8; x < width; x += 8) { \
        satd_8x8_subblock_ ## suffix(pred_ptrs, stride, orig_ptr, orig_stride, sums); \
        orig_ptr += 8; \
        pred_ptrs[0] += 8; \
        pred_ptrs[1] += 8; \
        pred_ptrs[2] += 8; \
        pred_ptrs[3] += 8; \
        costs_out[0] += sums[0]; \
        costs_out[1] += sums[1]; \
        costs_out[2] += sums[2]; \
        costs_out[3] += sums[3]; \
      } \
    } \
    for(int i = 0; i < num_parallel_blocks; ++i){\
      costs_out[i] = costs_out[i] >> (KVZ_BIT_DEPTH - 8);\
    } \
    return; \
  }

SATD_ANY_SIZE_MULTI_AVX2(quad_avx2, 4)


unsigned pixels_calc_ssd_avx2(const uint8_t *const ref, const uint8_t *const rec,
                 const int ref_stride, const int rec_stride,
                 const int width)
{
  __m256i ssd_part;
  __m256i diff = _mm256_setzero_si256();
  __m128i sum;

  __m256i ref_epi16;
  __m256i rec_epi16;

  __m128i ref_row0, ref_row1, ref_row2, ref_row3;
  __m128i rec_row0, rec_row1, rec_row2, rec_row3;

  int ssd;

  switch (width) {

  case 4:

    ref_row0 = _mm_cvtsi32_si128(*(int32_t*)&(ref[0 * ref_stride]));
    ref_row1 = _mm_cvtsi32_si128(*(int32_t*)&(ref[1 * ref_stride]));
    ref_row2 = _mm_cvtsi32_si128(*(int32_t*)&(ref[2 * ref_stride]));
    ref_row3 = _mm_cvtsi32_si128(*(int32_t*)&(ref[3 * ref_stride]));

    ref_row0 = _mm_unpacklo_epi32(ref_row0, ref_row1);
    ref_row1 = _mm_unpacklo_epi32(ref_row2, ref_row3);
    ref_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(ref_row0, ref_row1) );

    rec_row0 = _mm_cvtsi32_si128(*(int32_t*)&(rec[0 * rec_stride]));
    rec_row1 = _mm_cvtsi32_si128(*(int32_t*)&(rec[1 * rec_stride]));
    rec_row2 = _mm_cvtsi32_si128(*(int32_t*)&(rec[2 * rec_stride]));
    rec_row3 = _mm_cvtsi32_si128(*(int32_t*)&(rec[3 * rec_stride]));

    rec_row0 = _mm_unpacklo_epi32(rec_row0, rec_row1);
    rec_row1 = _mm_unpacklo_epi32(rec_row2, rec_row3);
    rec_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(rec_row0, rec_row1) );

    diff = _mm256_sub_epi16(ref_epi16, rec_epi16);
    ssd_part =  _mm256_madd_epi16(diff, diff);

    sum = _mm_add_epi32(_mm256_castsi256_si128(ssd_part), _mm256_extracti128_si256(ssd_part, 1));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));

    ssd = _mm_cvtsi128_si32(sum);

    return ssd >> (2*(KVZ_BIT_DEPTH-8));
    break;

  default:

    ssd_part = _mm256_setzero_si256();
    for (int y = 0; y < width; y += 8) {
      for (int x = 0; x < width; x += 8) {
        for (int i = 0; i < 8; i += 2) {
          ref_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)&(ref[x + (y + i) * ref_stride])), _mm_loadl_epi64((__m128i*)&(ref[x + (y + i + 1) * ref_stride]))));
          rec_epi16 = _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)&(rec[x + (y + i) * rec_stride])), _mm_loadl_epi64((__m128i*)&(rec[x + (y + i + 1) * rec_stride]))));
          diff = _mm256_sub_epi16(ref_epi16, rec_epi16);
          ssd_part = _mm256_add_epi32(ssd_part, _mm256_madd_epi16(diff, diff));
        }
      }
    }

    sum = _mm_add_epi32(_mm256_castsi256_si128(ssd_part), _mm256_extracti128_si256(ssd_part, 1));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(1, 0, 3, 2)));
    sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, _MM_SHUFFLE(0, 1, 0, 1)));

    ssd = _mm_cvtsi128_si32(sum);

    return ssd >> (2*(KVZ_BIT_DEPTH-8));
    break;
  }
}

inline void scatter_ymm_4x8_8bit(kvz_pixel * dst, __m256i ymm, unsigned dst_stride)
{
  __m128i ymm_lo = _mm256_castsi256_si128(ymm);
  __m128i ymm_hi = _mm256_extracti128_si256(ymm, 1);
  *(uint32_t *)dst = _mm_cvtsi128_si32(ymm_lo); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 1); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 2); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_lo, 3); dst += dst_stride;
  *(uint32_t *)dst = _mm_cvtsi128_si32(ymm_hi); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 1); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 2); dst += dst_stride;
  *(uint32_t *)dst = _mm_extract_epi32(ymm_hi, 3);
}

inline void scatter_ymm_8x4_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m256d ymm_as_m256d = _mm256_castsi256_pd(ymm);
  __m128d ymm_lo = _mm256_castpd256_pd128(ymm_as_m256d);
  __m128d ymm_hi = _mm256_extractf128_pd(ymm_as_m256d, 1);
  _mm_storel_pd((double*)dst, ymm_lo); dst += dst_stride;
  _mm_storeh_pd((double*)dst, ymm_lo); dst += dst_stride;
  _mm_storel_pd((double*)dst, ymm_hi); dst += dst_stride;
  _mm_storeh_pd((double*)dst, ymm_hi);
}

inline void scatter_ymm_16x2_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m128i ymm_lo = _mm256_castsi256_si128(ymm);
  __m128i ymm_hi = _mm256_extracti128_si256(ymm, 1);
  _mm_storeu_si128((__m128i *)dst, ymm_lo); dst += dst_stride;
  _mm_storeu_si128((__m128i *)dst, ymm_hi);
}

inline void scatter_ymm_12x2_8bit(kvz_pixel *dst, __m256i ymm, unsigned dst_stride)
{
  __m256i mask_a = _mm256_setr_epi32(-1, -1, -1, 0, 0, 0, 0, 0);
  __m256i mask_b = _mm256_setr_epi32(0, 0, 0, -1, -1, -1, 0, 0);
  _mm256_maskstore_epi32((int32_t*)dst, mask_a, ymm); dst += dst_stride - 3 * 4;
  _mm256_maskstore_epi32((int32_t*)dst, mask_b, ymm);
}

template <unsigned W>
void bipred_average_px_px_template_avx2(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_h,
  unsigned dst_stride)
{
  constexpr unsigned pu_w = W;
  constexpr bool has_pow2_width = (W & (W - 1)) == 0;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {

      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_L0 = _mm256_loadu_si256((__m256i*)&px_L0[i]);
      __m256i sample_L1 = _mm256_loadu_si256((__m256i*)&px_L1[i]);
      __m256i avg       = _mm256_avg_epu8(sample_L0, sample_L1);

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i *)&dst[y * dst_stride + x], avg); break;
        default:
          for (int py = 0; py < pu_h; ++py) {
            for (int px = 0; px < pu_w; ++px) {
              dst[py * dst_stride + px] = (px_L0[py * pu_w + px] + px_L1[py * pu_w + px] + 1) >> 1;
            }
          }
          return;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 are not used to simplify the loop
      __m256i mask      = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m256i sample_L0 = _mm256_maskload_epi64((const long long*)&px_L0[i], mask);
      __m256i sample_L1 = _mm256_maskload_epi64((const long long*)&px_L1[i], mask);
      __m256i avg       = _mm256_avg_epu8(sample_L0, sample_L1);

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          for (int py = 0; py < pu_h; ++py) {
            for (int px = 0; px < pu_w; ++px) {
              dst[py * dst_stride + px] = (px_L0[py * pu_w + px] + px_L1[py * pu_w + px] + 1) >> 1;
            }
          }
          return;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      __m128i sample_L0, sample_L1, avg;
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          sample_L0 = _mm_loadu_si128((__m128i*)&px_L0[i]);
          sample_L1 = _mm_loadu_si128((__m128i*)&px_L1[i]);
          avg       = _mm_avg_epu8(sample_L0, sample_L1);
          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m128i mask      = _mm_setr_epi32(-1, -1, -1, 0);
          __m128i sample_L0 = _mm_maskload_epi32((const int*)(&px_L0[i]), mask);
          __m128i sample_L1 = _mm_maskload_epi32((const int*)(&px_L1[i]), mask);
          __m128i avg       = _mm_avg_epu8(sample_L0, sample_L1);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        for (int py = 0; py < pu_h; ++py) {
          for (int px = 0; px < pu_w; ++px) {
            dst[py * dst_stride + px] = (px_L0[py * pu_w + px] + px_L1[py * pu_w + px] + 1) >> 1;
          }
        }
        return;
    }
  }
}

inline void bipred_average_px_px_avx2(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_px_px_template_avx2<4>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case  8: bipred_average_px_px_template_avx2<8>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 16: bipred_average_px_px_template_avx2<16>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 32: bipred_average_px_px_template_avx2<32>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 64: bipred_average_px_px_template_avx2<64>(dst, px_L0, px_L1, pu_h, dst_stride); break;

      case  6: bipred_average_px_px_template_avx2<6>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 12: bipred_average_px_px_template_avx2<12>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 24: bipred_average_px_px_template_avx2<24>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      case 48: bipred_average_px_px_template_avx2<48>(dst, px_L0, px_L1, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
    int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
    int32_t offset = 1 << (shift - 1);

    for (int i = 0; i < pu_w * pu_h; ++i)
    {
      int y = i / pu_w;
      int x = i % pu_w;
      int16_t sample_L0 = px_L0[i] << (14 - KVZ_BIT_DEPTH);
      int16_t sample_L1 = px_L1[i] << (14 - KVZ_BIT_DEPTH);
      int32_t rounded = (sample_L0 + sample_L1 + offset) >> shift;
      dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
    }
  }
}

template <unsigned W>
void bipred_average_im_im_template_avx2(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_h,
  unsigned dst_stride)
{
  constexpr unsigned pu_w = W;
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t scalar_offset = 1 << (shift - 1);
  __m256i offset = _mm256_set1_epi32(scalar_offset);

  constexpr bool has_pow2_width = (W & (W - 1)) == 0;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {
      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_L0_a_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
      __m256i sample_L1_a_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);
      __m256i sample_L0_b_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i + 16]);
      __m256i sample_L1_b_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i + 16]);

      __m256i sample_L0_L1_a_lo = _mm256_unpacklo_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_a_hi = _mm256_unpackhi_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_b_lo = _mm256_unpacklo_epi16(sample_L0_b_16bit, sample_L1_b_16bit);
      __m256i sample_L0_L1_b_hi = _mm256_unpackhi_epi16(sample_L0_b_16bit, sample_L1_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_L0_L1_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_L0_L1_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i*)&dst[y * dst_stride + x], avg); break;
        default:
          for (int py = 0; py < pu_h; ++py) {
            for (int px = 0; px < pu_w; ++px) {
              int idx = py * pu_w + px;
              int16_t sample_L0 = im_L0[idx];
              int16_t sample_L1 = im_L1[idx];
              int32_t rounded = (sample_L0 + sample_L1 + scalar_offset) >> shift;
              dst[py * dst_stride + px] = kvz_fast_clip_32bit_to_pixel(rounded);
            }
          }
          return;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 are not used to simplify the loop
      __m256i mask              = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m256i sample_L0_a_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
      __m256i sample_L1_a_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);
      __m256i sample_L0_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im_L0[i + 16]));
      __m256i sample_L1_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im_L1[i + 16]));

      __m256i sample_L0_L1_a_lo = _mm256_unpacklo_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_a_hi = _mm256_unpackhi_epi16(sample_L0_a_16bit, sample_L1_a_16bit);
      __m256i sample_L0_L1_b_lo = _mm256_unpacklo_epi16(sample_L0_b_16bit, sample_L1_b_16bit);
      __m256i sample_L0_L1_b_hi = _mm256_unpackhi_epi16(sample_L0_b_16bit, sample_L1_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_L0_L1_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_L0_L1_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          for (int py = 0; py < pu_h; ++py) {
            for (int px = 0; px < pu_w; ++px) {
              int idx = py * pu_w + px;
              int16_t sample_L0 = im_L0[idx];
              int16_t sample_L1 = im_L1[idx];
              int32_t rounded = (sample_L0 + sample_L1 + scalar_offset) >> shift;
              dst[py * dst_stride + px] = kvz_fast_clip_32bit_to_pixel(rounded);
            }
          }
          return;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          __m256i sample_L0_16bit = _mm256_loadu_si256((__m256i*)&im_L0[i]);
          __m256i sample_L1_16bit = _mm256_loadu_si256((__m256i*)&im_L1[i]);

          __m256i sample_L0_L1_lo = _mm256_unpacklo_epi16(sample_L0_16bit, sample_L1_16bit);
          __m256i sample_L0_L1_hi = _mm256_unpackhi_epi16(sample_L0_16bit, sample_L1_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_L0_L1_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_L0_L1_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m256i mask            = _mm256_setr_epi64x(-1, -1, -1, 0);
          __m256i sample_L0_16bit = _mm256_maskload_epi64((const long long*)(&im_L0[i]), mask);
          __m256i sample_L1_16bit = _mm256_maskload_epi64((const long long*)(&im_L1[i]), mask);

          __m256i sample_L0_L1_lo = _mm256_unpacklo_epi16(sample_L0_16bit, sample_L1_16bit);
          __m256i sample_L0_L1_hi = _mm256_unpackhi_epi16(sample_L0_16bit, sample_L1_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_a_lo = _mm256_madd_epi16(sample_L0_L1_lo, all_ones);
          __m256i avg_a_hi = _mm256_madd_epi16(sample_L0_L1_hi, all_ones);

          avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
          avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);

          avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
          avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        for (int py = 0; py < pu_h; ++py) {
          for (int px = 0; px < pu_w; ++px) {
            int idx = py * pu_w + px;
            int16_t sample_L0 = im_L0[idx];
            int16_t sample_L1 = im_L1[idx];
            int32_t rounded = (sample_L0 + sample_L1 + scalar_offset) >> shift;
            dst[py * dst_stride + px] = kvz_fast_clip_32bit_to_pixel(rounded);
          }
        }
        return;
    }
  }
}

void bipred_average_im_im_avx2(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_im_im_template_avx2<4>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case  8: bipred_average_im_im_template_avx2<8>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 16: bipred_average_im_im_template_avx2<16>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 32: bipred_average_im_im_template_avx2<32>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 64: bipred_average_im_im_template_avx2<64>(dst, im_L0, im_L1, pu_h, dst_stride); break;

      case  6: bipred_average_im_im_template_avx2<6>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 12: bipred_average_im_im_template_avx2<12>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 24: bipred_average_im_im_template_avx2<24>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      case 48: bipred_average_im_im_template_avx2<48>(dst, im_L0, im_L1, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
    int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
    int32_t offset = 1 << (shift - 1);

    for (int i = 0; i < pu_w * pu_h; ++i)
    {
      int y = i / pu_w;
      int x = i % pu_w;
      int16_t sample_L0 = im_L0[i];
      int16_t sample_L1 = im_L1[i];
      int32_t rounded = (sample_L0 + sample_L1 + offset) >> shift;
      dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
    }
  }
}

template <unsigned W>
void bipred_average_px_im_template_avx2(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_h,
  unsigned dst_stride)
{
  constexpr unsigned pu_w = W;
  int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
  int32_t scalar_offset = 1 << (shift - 1);
  __m256i offset = _mm256_set1_epi32(scalar_offset);

  constexpr bool has_pow2_width = (W & (W - 1)) == 0;
  bool area_mod_32 = (pu_w * pu_h) % 32;
  assert(!(pu_w == 4 && pu_h == 4) && "Branch for 4x4 not yet implemented.");
  assert(!(pu_w == 2 && pu_h == 8) && "Branch for 2x8 not yet implemented.");

  if (has_pow2_width && area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 32) {

      int y = i / pu_w;
      int x = i % pu_w;

      __m256i sample_px_a_16bit = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i*)&px[i]));
      __m256i sample_px_b_16bit = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i*)&px[i + 16]));
      sample_px_a_16bit         = _mm256_slli_epi16(sample_px_a_16bit, 14 - KVZ_BIT_DEPTH);
      sample_px_b_16bit         = _mm256_slli_epi16(sample_px_b_16bit, 14 - KVZ_BIT_DEPTH);
      __m256i sample_im_a_16bit = _mm256_loadu_si256((__m256i*)&im[i]);
      __m256i sample_im_b_16bit = _mm256_loadu_si256((__m256i*)&im[i + 16]);

      __m256i sample_px_im_a_lo = _mm256_unpacklo_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_a_hi = _mm256_unpackhi_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_b_lo = _mm256_unpacklo_epi16(sample_px_b_16bit, sample_im_b_16bit);
      __m256i sample_px_im_b_hi = _mm256_unpackhi_epi16(sample_px_b_16bit, sample_im_b_16bit);

      __m256i all_ones  = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_px_im_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_px_im_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_px_im_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_px_im_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case  4: scatter_ymm_4x8_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case  8: scatter_ymm_8x4_8bit( &dst[y * dst_stride + x], avg, dst_stride); break;
        case 16: scatter_ymm_16x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 32: // Same as case 64
        case 64: _mm256_storeu_si256((__m256i*)&dst[y * dst_stride + x], avg); break;
        default:
          for (int p_y = 0; p_y < pu_h; ++p_y) {
            for (int p_x = 0; p_x < pu_w; ++p_x) {
              int idx = p_y * pu_w + p_x;
              int16_t sample_px = px[idx] << (14 - KVZ_BIT_DEPTH);
              int16_t sample_im = im[idx];
              int32_t rounded = (sample_px + sample_im + scalar_offset) >> shift;
              dst[p_y * dst_stride + p_x] = kvz_fast_clip_32bit_to_pixel(rounded);
            }
          }
          return;
      }
    }
  } else if (area_mod_32 == 0) {
    for (int i = 0; i < pu_w * pu_h; i += 24) {

      int y = i / pu_w;
      int x = i % pu_w;

      // Last 64 bits of the 256 / 32 bits of the 128 are not used to simplify the loop
      __m256i mask              = _mm256_setr_epi64x(-1, -1, -1, 0);
      __m128i sample_px_a_8bit  = _mm_loadu_si128((__m128i*)&px[i]);
      __m128i sample_px_b_8bit  = _mm_loadl_epi64((__m128i*)&px[i + 16]);
      __m256i sample_px_a_16bit = _mm256_cvtepu8_epi16(sample_px_a_8bit);
      __m256i sample_px_b_16bit = _mm256_cvtepu8_epi16(sample_px_b_8bit);
      sample_px_a_16bit         = _mm256_slli_epi16(sample_px_a_16bit, 14 - KVZ_BIT_DEPTH);
      sample_px_b_16bit         = _mm256_slli_epi16(sample_px_b_16bit, 14 - KVZ_BIT_DEPTH);
      __m256i sample_im_a_16bit = _mm256_loadu_si256((__m256i*)&im[i]);
      __m256i sample_im_b_16bit = _mm256_castsi128_si256(_mm_loadu_si128((__m128i*)&im[i + 16]));

      __m256i sample_px_im_a_lo = _mm256_unpacklo_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_a_hi = _mm256_unpackhi_epi16(sample_px_a_16bit, sample_im_a_16bit);
      __m256i sample_px_im_b_lo = _mm256_unpacklo_epi16(sample_px_b_16bit, sample_im_b_16bit);
      __m256i sample_px_im_b_hi = _mm256_unpackhi_epi16(sample_px_b_16bit, sample_im_b_16bit);

      __m256i all_ones = _mm256_set1_epi16(1);
      __m256i avg_a_lo = _mm256_madd_epi16(sample_px_im_a_lo, all_ones);
      __m256i avg_a_hi = _mm256_madd_epi16(sample_px_im_a_hi, all_ones);
      __m256i avg_b_lo = _mm256_madd_epi16(sample_px_im_b_lo, all_ones);
      __m256i avg_b_hi = _mm256_madd_epi16(sample_px_im_b_hi, all_ones);

      avg_a_lo = _mm256_add_epi32(avg_a_lo, offset);
      avg_a_hi = _mm256_add_epi32(avg_a_hi, offset);
      avg_b_lo = _mm256_add_epi32(avg_b_lo, offset);
      avg_b_hi = _mm256_add_epi32(avg_b_hi, offset);

      avg_a_lo = _mm256_srai_epi32(avg_a_lo, shift);
      avg_a_hi = _mm256_srai_epi32(avg_a_hi, shift);
      avg_b_lo = _mm256_srai_epi32(avg_b_lo, shift);
      avg_b_hi = _mm256_srai_epi32(avg_b_hi, shift);

      __m256i avg_01  = _mm256_packus_epi32(avg_a_lo, avg_a_hi);
      __m256i avg_23  = _mm256_packus_epi32(avg_b_lo, avg_b_hi);
      __m256i avg0213 = _mm256_packus_epi16(avg_01, avg_23);
      __m256i avg     = _mm256_permute4x64_epi64(avg0213, _MM_SHUFFLE(3, 1, 2, 0));

      switch (pu_w) {
        case 12: scatter_ymm_12x2_8bit(&dst[y * dst_stride + x], avg, dst_stride); break;
        case 24: // Same as case 48
        case 48: _mm256_maskstore_epi64((long long*)&dst[y * dst_stride + x], mask, avg); break;
        default:
          for (int p_y = 0; p_y < pu_h; ++p_y) {
            for (int p_x = 0; p_x < pu_w; ++p_x) {
              int idx = p_y * pu_w + p_x;
              int16_t sample_px = px[idx] << (14 - KVZ_BIT_DEPTH);
              int16_t sample_im = im[idx];
              int32_t rounded = (sample_px + sample_im + scalar_offset) >> shift;
              dst[p_y * dst_stride + p_x] = kvz_fast_clip_32bit_to_pixel(rounded);
            }
          }
          return;
      }
    }
  } else {
    // 8x2, 8x6, 6x8 blocks (and maybe 2x8 in the future)
    switch (pu_w) {
      case 8: // 8x2, 8x6
        for (int i = 0; i < pu_w * pu_h; i += 16) {

          int y = i / pu_w;

          __m128i sample_px_8bit  = _mm_loadu_si128((__m128i*)&px[i]);
          __m256i sample_px_16bit = _mm256_cvtepu8_epi16(sample_px_8bit);
          sample_px_16bit         = _mm256_slli_epi16(sample_px_16bit, 14 - KVZ_BIT_DEPTH);
          __m256i sample_im_16bit = _mm256_loadu_si256((__m256i*)&im[i]);

          __m256i sample_px_im_lo = _mm256_unpacklo_epi16(sample_px_16bit, sample_im_16bit);
          __m256i sample_px_im_hi = _mm256_unpackhi_epi16(sample_px_16bit, sample_im_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_px_im_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_px_im_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          _mm_storel_epi64((__m128i*)&dst[y * dst_stride], avg);
          _mm_storeh_pd((double*)&dst[(y + 1) * dst_stride], _mm_castsi128_pd(avg));
        }
        break;
      case 6: // 6x8
        for (int i = 0; i < pu_w * pu_h; i += 12) {

          int y = i / pu_w;

          __m128i mask128         = _mm_setr_epi32(-1, -1, -1, 0);
          __m128i sample_px_8bit  = _mm_maskload_epi32((const int*)(&px[i]), mask128);

          __m256i mask            = _mm256_setr_epi64x(-1, -1, -1, 0);
          __m256i sample_px_16bit = _mm256_cvtepu8_epi16(sample_px_8bit);
          sample_px_16bit         = _mm256_slli_epi16(sample_px_16bit, 14 - KVZ_BIT_DEPTH);
          __m256i sample_im_16bit = _mm256_maskload_epi64((const long long*)(&im[i]), mask);

          __m256i sample_px_im_lo = _mm256_unpacklo_epi16(sample_px_16bit, sample_im_16bit);
          __m256i sample_px_im_hi = _mm256_unpackhi_epi16(sample_px_16bit, sample_im_16bit);

          __m256i all_ones = _mm256_set1_epi16(1);
          __m256i avg_lo   = _mm256_madd_epi16(sample_px_im_lo, all_ones);
          __m256i avg_hi   = _mm256_madd_epi16(sample_px_im_hi, all_ones);

          avg_lo = _mm256_add_epi32(avg_lo, offset);
          avg_hi = _mm256_add_epi32(avg_hi, offset);

          avg_lo = _mm256_srai_epi32(avg_lo, shift);
          avg_hi = _mm256_srai_epi32(avg_hi, shift);

          __m256i avg256 = _mm256_packus_epi32(avg_lo, avg_hi);
          avg256         = _mm256_packus_epi16(avg256, avg256);
          avg256         = _mm256_permute4x64_epi64(avg256, _MM_SHUFFLE(3, 1, 2, 0));
          __m128i avg    = _mm256_castsi256_si128(avg256);

          uint32_t elements_0123 = _mm_cvtsi128_si32(avg);
          uint16_t elements_45   = _mm_extract_epi16(avg, 2);
          uint16_t elements_67   = _mm_extract_epi16(avg, 3);
          uint32_t elements_89ab = _mm_extract_epi32(avg, 2);
          *(uint32_t*)&dst[(y + 0) * dst_stride + 0] = elements_0123;
          *(uint16_t*)&dst[(y + 0) * dst_stride + 4] = elements_45;
          *(uint16_t*)&dst[(y + 1) * dst_stride + 0] = elements_67;
          *(uint32_t*)&dst[(y + 1) * dst_stride + 2] = elements_89ab;
        }
        break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  }
}

void bipred_average_px_im_avx2(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  // Use scalar code for yet unoptimized block sizes (4x4, 2x8)
  if (!(pu_w == 4 && pu_h == 4) && pu_w > 2) {
    switch (pu_w) {
      case  4: bipred_average_px_im_template_avx2<4>(dst, px, im, pu_h, dst_stride); break;
      case  8: bipred_average_px_im_template_avx2<8>(dst, px, im, pu_h, dst_stride); break;
      case 16: bipred_average_px_im_template_avx2<16>(dst, px, im, pu_h, dst_stride); break;
      case 32: bipred_average_px_im_template_avx2<32>(dst, px, im, pu_h, dst_stride); break;
      case 64: bipred_average_px_im_template_avx2<64>(dst, px, im, pu_h, dst_stride); break;

      case  6: bipred_average_px_im_template_avx2<6>(dst, px, im, pu_h, dst_stride); break;
      case 12: bipred_average_px_im_template_avx2<12>(dst, px, im, pu_h, dst_stride); break;
      case 24: bipred_average_px_im_template_avx2<24>(dst, px, im, pu_h, dst_stride); break;
      case 48: bipred_average_px_im_template_avx2<48>(dst, px, im, pu_h, dst_stride); break;
      default:
        assert(0 && "Unexpected block width.");
        break;
    }
  } else {
    int32_t shift = 15 - KVZ_BIT_DEPTH; // TODO: defines
    int32_t offset = 1 << (shift - 1);

    for (int i = 0; i < pu_w * pu_h; ++i)
    {
      int y = i / pu_w;
      int x = i % pu_w;
      int16_t sample_px = px[i] << (14 - KVZ_BIT_DEPTH);
      int16_t sample_im = im[i];
      int32_t rounded = (sample_px + sample_im + offset) >> shift;
      dst[y * dst_stride + x] = kvz_fast_clip_32bit_to_pixel(rounded);
    }
  }
}

void bipred_average_avx2(lcu_t *const lcu,
  const yuv_t *const px_L0,
  const yuv_t *const px_L1,
  const yuv_im_t *const im_L0,
  const yuv_im_t *const im_L1,
  const unsigned pu_x,
  const unsigned pu_y,
  const unsigned pu_w,
  const unsigned pu_h,
  const unsigned im_flags_L0,
  const unsigned im_flags_L1,
  const bool predict_luma,
  const bool predict_chroma,
  const uint8_t chroma_shift_w,
  const uint8_t chroma_shift_h) {

  //After reconstruction, merge the predictors by taking an average of each pixel
  if (predict_luma) {
    unsigned pb_offset = SUB_SCU(pu_y) * LCU_WIDTH + SUB_SCU(pu_x);

    if (!(im_flags_L0 & 1) && !(im_flags_L1 & 1)) {
      bipred_average_px_px_avx2(lcu->rec.y + pb_offset, px_L0->y, px_L1->y, pu_w, pu_h, LCU_WIDTH);

    } else if ((im_flags_L0 & 1) && (im_flags_L1 & 1)) {
      bipred_average_im_im_avx2(lcu->rec.y + pb_offset, im_L0->y, im_L1->y, pu_w, pu_h, LCU_WIDTH);

    } else {
      kvz_pixel *src_px    = (im_flags_L0 & 1) ? px_L1->y : px_L0->y;
      kvz_pixel_im *src_im = (im_flags_L0 & 1) ? im_L0->y : im_L1->y;
      bipred_average_px_im_avx2(lcu->rec.y + pb_offset, src_px, src_im, pu_w, pu_h, LCU_WIDTH);
    }
  }
  if (predict_chroma) {
    unsigned pb_offset = ((SUB_SCU(pu_y) >> chroma_shift_h) * (LCU_WIDTH >> chroma_shift_w)) + (SUB_SCU(pu_x) >> chroma_shift_w);
    unsigned pb_w = pu_w >> chroma_shift_w;
    unsigned pb_h = pu_h >> chroma_shift_h;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_px_px_avx2(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_px_avx2(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);

    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_im_im_avx2(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_im_im_avx2(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);

    } else {
      kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_px_im_avx2(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_im_avx2(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    }
  }
}

optimized_sad_func_ptr_t get_optimized_sad_avx2(int32_t width)
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
  if (width == 32)
    return reg_sad_w32;
  if (width == 64)
    return reg_sad_w64;
  else
    return NULL;
}

uint32_t ver_sad_avx2(const uint8_t *pic_data, const uint8_t *ref_data,
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

uint32_t hor_sad_avx2(const uint8_t *pic_data, const uint8_t *ref_data,
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
    return hor_sad_avx2_w32 (pic_data, ref_data, height,
                             pic_stride, ref_stride, left, right);
  else
    return hor_sad_sse41_arbitrary(pic_data, ref_data, width, height,
                                   pic_stride, ref_stride, left, right);
}

double pixel_var_avx2_largebuf(const uint8_t *buf, const uint32_t len)
{
  const float len_f  = (float)len;
  const __m256i zero = _mm256_setzero_si256();

  int64_t sum;
  size_t i;
  __m256i sums = zero;
  for (i = 0; i + 31 < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));
    __m256i curr_sum = _mm256_sad_epu8(curr, zero);
            sums = _mm256_add_epi64(sums, curr_sum);
  }
  __m128i sum_lo = _mm256_castsi256_si128  (sums);
  __m128i sum_hi = _mm256_extracti128_si256(sums,   1);
  __m128i sum_3  = _mm_add_epi64           (sum_lo, sum_hi);
  __m128i sum_4  = _mm_shuffle_epi32       (sum_3,  _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sum_5  = _mm_add_epi64           (sum_3,  sum_4);

  _mm_storel_epi64((__m128i *)&sum, sum_5);

  // Remaining len mod 32 pixels
  for (; i < len; ++i) {
    sum += buf[i];
  }

  float   mean_f = (float)sum / len_f;
  __m256  mean   = _mm256_set1_ps(mean_f);
  __m256  accum  = _mm256_setzero_ps();

  for (i = 0; i + 31 < len; i += 32) {
    __m128i curr0    = _mm_loadl_epi64((const __m128i *)(buf + i +  0));
    __m128i curr1    = _mm_loadl_epi64((const __m128i *)(buf + i +  8));
    __m128i curr2    = _mm_loadl_epi64((const __m128i *)(buf + i + 16));
    __m128i curr3    = _mm_loadl_epi64((const __m128i *)(buf + i + 24));

    __m256i curr0_32 = _mm256_cvtepu8_epi32(curr0);
    __m256i curr1_32 = _mm256_cvtepu8_epi32(curr1);
    __m256i curr2_32 = _mm256_cvtepu8_epi32(curr2);
    __m256i curr3_32 = _mm256_cvtepu8_epi32(curr3);

    __m256  curr0_f  = _mm256_cvtepi32_ps  (curr0_32);
    __m256  curr1_f  = _mm256_cvtepi32_ps  (curr1_32);
    __m256  curr2_f  = _mm256_cvtepi32_ps  (curr2_32);
    __m256  curr3_f  = _mm256_cvtepi32_ps  (curr3_32);

    __m256  curr0_sd = _mm256_sub_ps       (curr0_f,  mean);
    __m256  curr1_sd = _mm256_sub_ps       (curr1_f,  mean);
    __m256  curr2_sd = _mm256_sub_ps       (curr2_f,  mean);
    __m256  curr3_sd = _mm256_sub_ps       (curr3_f,  mean);

    __m256  curr0_v  = _mm256_mul_ps       (curr0_sd, curr0_sd);
    __m256  curr1_v  = _mm256_mul_ps       (curr1_sd, curr1_sd);
    __m256  curr2_v  = _mm256_mul_ps       (curr2_sd, curr2_sd);
    __m256  curr3_v  = _mm256_mul_ps       (curr3_sd, curr3_sd);

    __m256  curr01   = _mm256_add_ps       (curr0_v,  curr1_v);
    __m256  curr23   = _mm256_add_ps       (curr2_v,  curr3_v);
    __m256  curr     = _mm256_add_ps       (curr01,   curr23);
            accum    = _mm256_add_ps       (accum,    curr);
  }
  __m256d accum_d  = _mm256_castps_pd     (accum);
  __m256d accum2_d = _mm256_permute4x64_pd(accum_d, _MM_SHUFFLE(1, 0, 3, 2));
  __m256  accum2   = _mm256_castpd_ps     (accum2_d);

  __m256  accum3   = _mm256_add_ps        (accum,  accum2);
  __m256  accum4   = _mm256_permute_ps    (accum3, _MM_SHUFFLE(1, 0, 3, 2));
  __m256  accum5   = _mm256_add_ps        (accum3, accum4);
  __m256  accum6   = _mm256_permute_ps    (accum5, _MM_SHUFFLE(2, 3, 0, 1));
  __m256  accum7   = _mm256_add_ps        (accum5, accum6);

  __m128  accum8   = _mm256_castps256_ps128(accum7);
  float   var_sum  = _mm_cvtss_f32         (accum8);

  // Remaining len mod 32 pixels
  for (; i < len; ++i) {
    float diff = buf[i] - mean_f;
    var_sum += diff * diff;
  }

  return  var_sum / len_f;
}

#ifdef INACCURATE_VARIANCE_CALCULATION

// Assumes that u is a power of two
inline uint32_t ilog2(uint32_t u)
{
  return _tzcnt_u32(u);
}

// A B C D | E F G H (8x32b)
//        ==>
// A+B C+D | E+F G+H (4x64b)
__m256i hsum_epi32_to_epi64(const __m256i v)
{
  const __m256i zero    = _mm256_setzero_si256();
        __m256i v_shufd = _mm256_shuffle_epi32(v, _MM_SHUFFLE(3, 3, 1, 1));
        __m256i sums_32 = _mm256_add_epi32    (v, v_shufd);
        __m256i sums_64 = _mm256_blend_epi32  (sums_32, zero, 0xaa);
  return        sums_64;
}

double pixel_var_avx2(const uint8_t *buf, const uint32_t len)
{
  assert(sizeof(*buf) == 1);
  assert((len & 31) == 0);

  // Uses Q8.7 numbers to measure mean and deviation, so variances are Q16.14
  const uint64_t sum_maxwid     = ilog2(len) + (8 * sizeof(*buf));
  const __m128i normalize_sum   = _mm_cvtsi32_si128(sum_maxwid - 15); // Normalize mean to [0, 32767], so signed 16-bit subtraction never overflows
  const __m128i debias_sum      = _mm_cvtsi32_si128(1 << (sum_maxwid - 16));
  const float varsum_to_f       = 1.0f / (float)(1 << (14 + ilog2(len)));

  const bool power_of_two = (len & (len - 1)) == 0;
  if (sum_maxwid > 32 || sum_maxwid < 15 || !power_of_two) {
    return pixel_var_avx2_largebuf(buf, len);
  }

  const __m256i zero      = _mm256_setzero_si256();
  const __m256i himask_15 = _mm256_set1_epi16(0x7f00);

  uint64_t vars;
  size_t i;
  __m256i sums = zero;
  for (i = 0; i < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));
    __m256i curr_sum = _mm256_sad_epu8(curr, zero);
            sums = _mm256_add_epi64(sums, curr_sum);
  }
  __m128i sum_lo = _mm256_castsi256_si128  (sums);
  __m128i sum_hi = _mm256_extracti128_si256(sums,   1);
  __m128i sum_3  = _mm_add_epi64           (sum_lo, sum_hi);
  __m128i sum_4  = _mm_shuffle_epi32       (sum_3,  _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sum_5  = _mm_add_epi64           (sum_3,  sum_4);
  __m128i sum_5n = _mm_srl_epi32           (sum_5,  normalize_sum);
          sum_5n = _mm_add_epi32           (sum_5n, debias_sum);

  __m256i sum_n  = _mm256_broadcastw_epi16 (sum_5n);

  __m256i accum = zero;
  for (i = 0; i < len; i += 32) {
    __m256i curr = _mm256_loadu_si256((const __m256i *)(buf + i));

    __m256i curr0    = _mm256_slli_epi16  (curr,  7);
    __m256i curr1    = _mm256_srli_epi16  (curr,  1);
            curr0    = _mm256_and_si256   (curr0, himask_15);
            curr1    = _mm256_and_si256   (curr1, himask_15);

    __m256i dev0     = _mm256_sub_epi16   (curr0, sum_n);
    __m256i dev1     = _mm256_sub_epi16   (curr1, sum_n);

    __m256i vars0    = _mm256_madd_epi16  (dev0,  dev0);
    __m256i vars1    = _mm256_madd_epi16  (dev1,  dev1);

    __m256i varsum   = _mm256_add_epi32   (vars0, vars1);
            varsum   = hsum_epi32_to_epi64(varsum);
            accum    = _mm256_add_epi64   (accum, varsum);
  }
  __m256i accum2 = _mm256_permute4x64_epi64(accum,  _MM_SHUFFLE(1, 0, 3, 2));
  __m256i accum3 = _mm256_add_epi64        (accum,  accum2);
  __m256i accum4 = _mm256_permute4x64_epi64(accum3, _MM_SHUFFLE(2, 3, 1, 0));
  __m256i v_tot  = _mm256_add_epi64        (accum3, accum4);
  __m128i vt128  = _mm256_castsi256_si128  (v_tot);

  _mm_storel_epi64((__m128i *)&vars, vt128);

  return (float)vars * varsum_to_f;
}

#else // INACCURATE_VARIANCE_CALCULATION

double pixel_var_avx2(const uint8_t *buf, const uint32_t len)
{
  return pixel_var_avx2_largebuf(buf, len);
}

#endif // !INACCURATE_VARIANCE_CALCULATION

static unsigned bilinear_ssd_avx2(const kvz_pixel *org, int org_stride,
                                  const kvz_pixel *ref, int ref_stride,
                                  int size, int fx, int fy)
{
  const __m256i w00 = _mm256_set1_epi16((int16_t)((4 - fx) * (4 - fy)));
  const __m256i w10 = _mm256_set1_epi16((int16_t)(fx * (4 - fy)));
  const __m256i w01 = _mm256_set1_epi16((int16_t)((4 - fx) * fy));
  const __m256i w11 = _mm256_set1_epi16((int16_t)(fx * fy));
  const __m256i rnd = _mm256_set1_epi16(8);
  __m256i acc = _mm256_setzero_si256();

  if (size == 16) {
    for (int y = 0; y < 16; y++) {
      const kvz_pixel *r0 = &ref[y * ref_stride];
      const kvz_pixel *r1 = r0 + ref_stride;
      __m256i a = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)r0));
      __m256i b = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(r0 + 1)));
      __m256i c = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)r1));
      __m256i d = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(r1 + 1)));
      __m256i o = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&org[y * org_stride]));
      __m256i p = _mm256_add_epi16(_mm256_mullo_epi16(a, w00), _mm256_mullo_epi16(b, w10));
      p = _mm256_add_epi16(p, _mm256_mullo_epi16(c, w01));
      p = _mm256_add_epi16(p, _mm256_mullo_epi16(d, w11));
      p = _mm256_srli_epi16(_mm256_add_epi16(p, rnd), 4);
      __m256i diff = _mm256_sub_epi16(o, p);
      acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
    }
  } else {
    // 8x8: two rows per iteration.
    for (int y = 0; y < size; y += 2) {
      const kvz_pixel *r0 = &ref[y * ref_stride];
      const kvz_pixel *r1 = r0 + ref_stride;
      const kvz_pixel *r2 = r1 + ref_stride;
      const kvz_pixel *o0 = &org[y * org_stride];
#define LOAD2(p, q) _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)(p)), \
                                                              _mm_loadl_epi64((const __m128i *)(q))))
      __m256i a = LOAD2(r0, r1);
      __m256i b = LOAD2(r0 + 1, r1 + 1);
      __m256i c = LOAD2(r1, r2);
      __m256i d = LOAD2(r1 + 1, r2 + 1);
      __m256i o = LOAD2(o0, o0 + org_stride);
#undef LOAD2
      __m256i p = _mm256_add_epi16(_mm256_mullo_epi16(a, w00), _mm256_mullo_epi16(b, w10));
      p = _mm256_add_epi16(p, _mm256_mullo_epi16(c, w01));
      p = _mm256_add_epi16(p, _mm256_mullo_epi16(d, w11));
      p = _mm256_srli_epi16(_mm256_add_epi16(p, rnd), 4);
      __m256i diff = _mm256_sub_epi16(o, p);
      acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
    }
  }
  __m128i s = _mm_add_epi32(_mm256_castsi256_si128(acc), _mm256_extracti128_si256(acc, 1));
  s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(1, 0, 3, 2)));
  s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(2, 3, 0, 1)));
  return (unsigned)_mm_cvtsi128_si32(s);
}


static INLINE __m128i load_4x2_u8(const kvz_pixel *a, const kvz_pixel *b)
{
  int32_t x, y;
  memcpy(&x, a, 4);
  memcpy(&y, b, 4);
  return _mm_unpacklo_epi32(_mm_cvtsi32_si128(x), _mm_cvtsi32_si128(y));
}


static void tf_bilinear_pred_avx2(kvz_pixel *dst, int dst_stride,
                                  const kvz_pixel *ref, int ref_stride,
                                  int width, int height,
                                  int w00, int w10, int w01, int w11, int shift)
{
  if (width != 8 && width != 4) {
    kvz_tf_bilinear_pred_generic(dst, dst_stride, ref, ref_stride, width, height,
                                 w00, w10, w01, w11, shift);
    return;
  }
  const __m256i v00 = _mm256_set1_epi16((int16_t)w00);
  const __m256i v10 = _mm256_set1_epi16((int16_t)w10);
  const __m256i v01 = _mm256_set1_epi16((int16_t)w01);
  const __m256i v11 = _mm256_set1_epi16((int16_t)w11);
  const __m256i rnd = _mm256_set1_epi16((int16_t)(1 << (shift - 1)));
  const __m128i sh = _mm_cvtsi32_si128(shift);
  const int rows = width == 8 ? 2 : 4;
  int y = 0;
  for (; y + rows <= height; y += rows) {
    const kvz_pixel *r0 = &ref[y * ref_stride];
    __m256i a, b, c, d;
    if (width == 8) {
      const kvz_pixel *r1 = r0 + ref_stride;
      const kvz_pixel *r2 = r1 + ref_stride;
#define LOAD2(p, q) _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)(p)), \
                                                              _mm_loadl_epi64((const __m128i *)(q))))
      a = LOAD2(r0, r1);
      b = LOAD2(r0 + 1, r1 + 1);
      c = LOAD2(r1, r2);
      d = LOAD2(r1 + 1, r2 + 1);
#undef LOAD2
    } else {
      const kvz_pixel *rr[5];
      for (int i = 0; i < 5; i++) rr[i] = r0 + i * ref_stride;
#define LOAD4(o) _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(load_4x2_u8(rr[0 + (o)], rr[1 + (o)]), \
                                                           load_4x2_u8(rr[2 + (o)], rr[3 + (o)])))
#define LOAD4X(o) _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(load_4x2_u8(rr[0 + (o)] + 1, rr[1 + (o)] + 1), \
                                                            load_4x2_u8(rr[2 + (o)] + 1, rr[3 + (o)] + 1)))
      a = LOAD4(0);
      b = LOAD4X(0);
      c = LOAD4(1);
      d = LOAD4X(1);
#undef LOAD4
#undef LOAD4X
    }
    __m256i p = _mm256_add_epi16(_mm256_mullo_epi16(a, v00), _mm256_mullo_epi16(b, v10));
    p = _mm256_add_epi16(p, _mm256_mullo_epi16(c, v01));
    p = _mm256_add_epi16(p, _mm256_mullo_epi16(d, v11));
    p = _mm256_srl_epi16(_mm256_add_epi16(p, rnd), sh);
    const __m256i packed = _mm256_packus_epi16(p, p);
    const __m128i lo = _mm256_castsi256_si128(packed);
    const __m128i hi = _mm256_extracti128_si256(packed, 1);
    if (width == 8) {
      _mm_storel_epi64((__m128i *)&dst[y * dst_stride], lo);
      _mm_storel_epi64((__m128i *)&dst[(y + 1) * dst_stride], hi);
    } else {
      const int32_t v0 = _mm_cvtsi128_si32(lo), v1 = _mm_extract_epi32(lo, 1);
      const int32_t v2 = _mm_cvtsi128_si32(hi), v3 = _mm_extract_epi32(hi, 1);
      memcpy(&dst[y * dst_stride], &v0, 4);
      memcpy(&dst[(y + 1) * dst_stride], &v1, 4);
      memcpy(&dst[(y + 2) * dst_stride], &v2, 4);
      memcpy(&dst[(y + 3) * dst_stride], &v3, 4);
    }
  }
  if (y < height) {
    kvz_tf_bilinear_pred_generic(&dst[y * dst_stride], dst_stride, &ref[y * ref_stride], ref_stride,
                                 width, height - y, w00, w10, w01, w11, shift);
  }
}


static INLINE int32_t hsum_epi32_128(__m128i v)
{
  v = _mm_add_epi32(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2)));
  v = _mm_add_epi32(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(2, 3, 0, 1)));
  return _mm_cvtsi128_si32(v);
}


static void tf_block_stats_avx2(const kvz_pixel *org, int org_stride,
                                const kvz_pixel *pred, int pred_stride,
                                int width, int height, int64_t out[5])
{
  if (width != 8 || height != 8) {
    kvz_tf_block_stats_generic(org, org_stride, pred, pred_stride, width, height, out);
    return;
  }
  const __m128i ones = _mm_set1_epi16(1);
  // Lanes 0..6: the pixels that have a right neighbour.
  const __m128i mask7 = _mm_setr_epi16(-1, -1, -1, -1, -1, -1, -1, 0);
  __m128i o[8], d[8];
  for (int y = 0; y < 8; y++) {
    o[y] = _mm_cvtepu8_epi16(_mm_loadl_epi64((const __m128i *)&org[y * org_stride]));
    const __m128i p = _mm_cvtepu8_epi16(_mm_loadl_epi64((const __m128i *)&pred[y * pred_stride]));
    d[y] = _mm_sub_epi16(o[y], p);
  }
  __m128i sum = _mm_setzero_si128(), sum2 = _mm_setzero_si128(), sse = _mm_setzero_si128();
  __m128i var_d = _mm_setzero_si128(), diff = _mm_setzero_si128();
  for (int y = 0; y < 8; y++) {
    sum = _mm_add_epi32(sum, _mm_madd_epi16(o[y], ones));
    sum2 = _mm_add_epi32(sum2, _mm_madd_epi16(o[y], o[y]));
    sse = _mm_add_epi32(sse, _mm_madd_epi16(d[y], d[y]));
    if (y < 7) {
      const __m128i dm = _mm_and_si128(d[y], mask7);
      var_d = _mm_add_epi32(var_d, _mm_madd_epi16(dm, dm));
      const __m128i dh = _mm_and_si128(_mm_sub_epi16(_mm_srli_si128(d[y], 2), d[y]), mask7);
      const __m128i dv = _mm_and_si128(_mm_sub_epi16(d[y + 1], d[y]), mask7);
      diff = _mm_add_epi32(diff, _mm_add_epi32(_mm_madd_epi16(dh, dh), _mm_madd_epi16(dv, dv)));
    }
  }
  out[0] = hsum_epi32_128(sum);
  out[1] = hsum_epi32_128(sum2);
  out[2] = hsum_epi32_128(sse);
  out[3] = hsum_epi32_128(var_d);
  out[4] = hsum_epi32_128(diff);
}


/** \brief Filter 8 pixels given as 8 bytes of src and of each prediction. */
static INLINE __m128i tf_filter_8px(__m128i src8, const __m128i *pred8, int num_refs,
                                    const float *weights, const float *const *luts,
                                    __m256i max_val)
{
  const __m256i o = _mm256_cvtepu8_epi32(src8);
  __m256 sum = _mm256_cvtepi32_ps(o);
  __m256 wsum = _mm256_set1_ps(1.0f);
  for (int r = 0; r < num_refs; r++) {
    const __m256i p = _mm256_cvtepu8_epi32(pred8[r]);
    const __m256 lut = _mm256_i32gather_ps(luts[r], _mm256_sub_epi32(p, o), 4);
    const __m256 wgt = _mm256_mul_ps(_mm256_set1_ps(weights[r]), lut);
    sum = _mm256_add_ps(sum, _mm256_mul_ps(wgt, _mm256_cvtepi32_ps(p)));
    wsum = _mm256_add_ps(wsum, wgt);
  }
  __m256i res = _mm256_cvttps_epi32(_mm256_add_ps(_mm256_div_ps(sum, wsum), _mm256_set1_ps(0.5f)));
  res = _mm256_min_epi32(_mm256_max_epi32(res, _mm256_setzero_si256()), max_val);
  const __m256i w16 = _mm256_packus_epi32(res, res);
  const __m256i b8 = _mm256_packus_epi16(w16, w16);
  return _mm_unpacklo_epi32(_mm256_castsi256_si128(b8), _mm256_extracti128_si256(b8, 1));
}


static void tf_filter_block_avx2(kvz_pixel *dst, int dst_stride,
                                 const kvz_pixel *src, int src_stride,
                                 const kvz_pixel *const *preds, int pred_stride,
                                 int num_refs, const float *weights,
                                 const float *const *luts,
                                 int width, int height, int max_val)
{
  if (num_refs > 8 || !((width % 8 == 0) || (width == 4 && height % 2 == 0))) {
    kvz_tf_filter_block_generic(dst, dst_stride, src, src_stride, preds, pred_stride,
                                num_refs, weights, luts, width, height, max_val);
    return;
  }
  const __m256i maxv = _mm256_set1_epi32(max_val);
  __m128i p8[8];
  if (width == 4) {
    for (int y = 0; y < height; y += 2) {
      for (int r = 0; r < num_refs; r++) {
        p8[r] = load_4x2_u8(&preds[r][y * pred_stride], &preds[r][(y + 1) * pred_stride]);
      }
      const __m128i res = tf_filter_8px(load_4x2_u8(&src[y * src_stride], &src[(y + 1) * src_stride]),
                                        p8, num_refs, weights, luts, maxv);
      const int32_t v0 = _mm_cvtsi128_si32(res), v1 = _mm_extract_epi32(res, 1);
      memcpy(&dst[y * dst_stride], &v0, 4);
      memcpy(&dst[(y + 1) * dst_stride], &v1, 4);
    }
    return;
  }
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x += 8) {
      for (int r = 0; r < num_refs; r++) {
        p8[r] = _mm_loadl_epi64((const __m128i *)&preds[r][y * pred_stride + x]);
      }
      const __m128i res = tf_filter_8px(_mm_loadl_epi64((const __m128i *)&src[y * src_stride + x]),
                                        p8, num_refs, weights, luts, maxv);
      _mm_storel_epi64((__m128i *)&dst[y * dst_stride + x], res);
    }
  }
}


static void bilinear_ssd_multi_avx2(const kvz_pixel *org, int org_stride,
                                    const kvz_pixel *ref, int ref_stride,
                                    int size, const int *mv_x, const int *mv_y,
                                    int num, unsigned *out)
{
  if (size != 8) {
    kvz_bilinear_ssd_multi_generic(org, org_stride, ref, ref_stride, size, mv_x, mv_y, num, out);
    return;
  }
#define LOAD2(p, q) _mm256_cvtepu8_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)(p)), \
                                                            _mm_loadl_epi64((const __m128i *)(q))))
  // The original block is loaded once, two rows per register.
  __m256i o[4];
  for (int i = 0; i < 4; i++) {
    o[i] = LOAD2(&org[2 * i * org_stride], &org[(2 * i + 1) * org_stride]);
  }
  const __m256i rnd = _mm256_set1_epi16(8);
  for (int n = 0; n < num; n++) {
    const int fx = mv_x[n] & 3, fy = mv_y[n] & 3;
    const kvz_pixel *r = ref + (mv_y[n] >> 2) * ref_stride + (mv_x[n] >> 2);
    __m256i acc = _mm256_setzero_si256();
    if (!fx && !fy) {
      for (int i = 0; i < 4; i++) {
        const kvz_pixel *r0 = &r[2 * i * ref_stride];
        const __m256i diff = _mm256_sub_epi16(o[i], LOAD2(r0, r0 + ref_stride));
        acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
      }
    } else {
      const __m256i w00 = _mm256_set1_epi16((int16_t)((4 - fx) * (4 - fy)));
      const __m256i w10 = _mm256_set1_epi16((int16_t)(fx * (4 - fy)));
      const __m256i w01 = _mm256_set1_epi16((int16_t)((4 - fx) * fy));
      const __m256i w11 = _mm256_set1_epi16((int16_t)(fx * fy));
      for (int i = 0; i < 4; i++) {
        const kvz_pixel *r0 = &r[2 * i * ref_stride];
        const kvz_pixel *r1 = r0 + ref_stride;
        const kvz_pixel *r2 = r1 + ref_stride;
        const __m256i a = LOAD2(r0, r1);
        const __m256i b = LOAD2(r0 + 1, r1 + 1);
        const __m256i c = LOAD2(r1, r2);
        const __m256i d = LOAD2(r1 + 1, r2 + 1);
        __m256i p = _mm256_add_epi16(_mm256_mullo_epi16(a, w00), _mm256_mullo_epi16(b, w10));
        p = _mm256_add_epi16(p, _mm256_mullo_epi16(c, w01));
        p = _mm256_add_epi16(p, _mm256_mullo_epi16(d, w11));
        p = _mm256_srli_epi16(_mm256_add_epi16(p, rnd), 4);
        const __m256i diff = _mm256_sub_epi16(o[i], p);
        acc = _mm256_add_epi32(acc, _mm256_madd_epi16(diff, diff));
      }
    }
    out[n] = (unsigned)hsum_epi32_128(_mm_add_epi32(_mm256_castsi256_si128(acc),
                                                    _mm256_extracti128_si256(acc, 1)));
  }
#undef LOAD2
}

} // namespace

#endif // KVZ_BIT_DEPTH == 8
#endif //COMPILE_INTEL_AVX2

int kvz_strategy_register_picture_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8
  if (bitdepth == KVZ_BIT_DEPTH) {
    success &= kvz_strategyselector_register(opaque, "reg_sad", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::reg_sad_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "sad_4x4", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_avx2<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_avx2<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_avx2<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_avx2<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_avx2<64>));
    success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_avx2<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_avx2<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_avx2<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_avx2<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_avx2<64>));
    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::pixels_calc_ssd_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "bipred_average", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::bipred_average_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::get_optimized_sad_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "ver_sad", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::ver_sad_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "hor_sad", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::hor_sad_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "pixel_var", "avx2", 40, KVZ_STRATEGY_FPTR(hbd_sad::pixel_var_hbd_avx2));

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "avx2", 40, KVZ_STRATEGY_FPTR(satd_4x4_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_avx2<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_avx2<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_avx2<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_hbd_avx2<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_4x4_dual_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_avx2<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_avx2<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_avx2<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_hbd_avx2<64>));
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "avx2", 40, KVZ_STRATEGY_FPTR(satd_any_size_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "avx2", 40, KVZ_STRATEGY_FPTR(satd_any_size_quad_hbd_avx2));

    success &= kvz_strategyselector_register(opaque, "bilinear_ssd", "avx2", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "bilinear_ssd_multi", "avx2", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_multi_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_bilinear_pred", "avx2", 40, KVZ_STRATEGY_FPTR(tf_bilinear_pred_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_block_stats", "avx2", 40, KVZ_STRATEGY_FPTR(tf_block_stats_hbd_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_filter_block", "avx2", 40, KVZ_STRATEGY_FPTR(tf_filter_block_hbd_avx2));
  }
#endif // COMPILE_INTEL_AVX2 && KVZ_BIT_DEPTH > 8
#if COMPILE_INTEL_AVX2
#if KVZ_BIT_DEPTH == 8
  // We don't actually use SAD for intra right now, other than 4x4 for
  // transform skip, but we might again one day and this is some of the
  // simplest code to look at for anyone interested in doing more
  // optimizations, so it's worth it to keep this maintained.
  if (bitdepth == 8){

    success &= kvz_strategyselector_register(opaque, "reg_sad", "avx2", 40, KVZ_STRATEGY_FPTR(kvz_reg_sad_avx2));
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "avx2", 40, KVZ_STRATEGY_FPTR(sad_8bit_8x8_avx2));
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "avx2", 40, KVZ_STRATEGY_FPTR(sad_8bit_16x16_avx2));
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "avx2", 40, KVZ_STRATEGY_FPTR(sad_8bit_32x32_avx2));
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "avx2", 40, KVZ_STRATEGY_FPTR(sad_8bit_64x64_avx2));

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "avx2", 40, KVZ_STRATEGY_FPTR(satd_4x4_8bit_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8x8_8bit_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "avx2", 40, KVZ_STRATEGY_FPTR(satd_16x16_8bit_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "avx2", 40, KVZ_STRATEGY_FPTR(satd_32x32_8bit_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "avx2", 40, KVZ_STRATEGY_FPTR(satd_64x64_8bit_avx2));

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8bit_4x4_dual_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8bit_8x8_dual_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8bit_16x16_dual_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8bit_32x32_dual_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "avx2", 40, KVZ_STRATEGY_FPTR(satd_8bit_64x64_dual_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "avx2", 40, KVZ_STRATEGY_FPTR(satd_any_size_8bit_avx2));
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "avx2", 40, KVZ_STRATEGY_FPTR(satd_any_size_quad_avx2));

    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "avx2", 40, KVZ_STRATEGY_FPTR(pixels_calc_ssd_avx2));
    success &= kvz_strategyselector_register(opaque, "bipred_average", "avx2", 40, KVZ_STRATEGY_FPTR(bipred_average_avx2));
    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "avx2", 40, KVZ_STRATEGY_FPTR(get_optimized_sad_avx2));
    success &= kvz_strategyselector_register(opaque, "ver_sad", "avx2", 40, KVZ_STRATEGY_FPTR(ver_sad_avx2));
    success &= kvz_strategyselector_register(opaque, "hor_sad", "avx2", 40, KVZ_STRATEGY_FPTR(hor_sad_avx2));

    success &= kvz_strategyselector_register(opaque, "pixel_var", "avx2", 40, KVZ_STRATEGY_FPTR(pixel_var_avx2));
    success &= kvz_strategyselector_register(opaque, "bilinear_ssd", "avx2", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_avx2));
    success &= kvz_strategyselector_register(opaque, "bilinear_ssd_multi", "avx2", 40, KVZ_STRATEGY_FPTR(bilinear_ssd_multi_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_bilinear_pred", "avx2", 40, KVZ_STRATEGY_FPTR(tf_bilinear_pred_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_block_stats", "avx2", 40, KVZ_STRATEGY_FPTR(tf_block_stats_avx2));
    success &= kvz_strategyselector_register(opaque, "tf_filter_block", "avx2", 40, KVZ_STRATEGY_FPTR(tf_filter_block_avx2));

  }
#endif // KVZ_BIT_DEPTH == 8
#endif
  return success;
}
