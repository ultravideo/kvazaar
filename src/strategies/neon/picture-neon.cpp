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

#include "strategies/strategies-cpp.h"

extern "C" {
#include "strategies/neon/picture-neon.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "strategies/neon/neon-common.h"
#include "strategies/optimized_sad_func_ptr_t.h"
#include "strategies/strategies-picture.h"
#include "strategyselector.h"
}

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH > 8

// In a namespace of their own to keep the names apart from the other
// high bit depth functions of this file.
namespace {
namespace hbd_sad {

/*
 * NEON versions of the SAD, SSD, variance and bi-prediction functions for
 * 10 and 12-bit pixels (kvz_pixel = uint16_t). Every function except
 * pixel_var computes exactly the same result as the generic version.
 *
 * SAD sums are accumulated with vaba into unsigned 16-bit lanes, which take
 * kHbdSadAdds absolute differences before they are flushed into 32-bit lanes.
 */
constexpr int kHbdSadAdds = 65535 / PIXEL_MAX;
static_assert(kHbdSadAdds >= 16, "SAD accumulation needs at least 16 additions per lane");

// Two 4-pixel rows in one vector.
KVZ_FORCE_INLINE uint16x8_t hbd_load_2x4(const kvz_pixel *p, unsigned stride)
{
  return vcombine_u16(vld1_u16(p), vld1_u16(p + stride));
}

// One 4-pixel row in the low half of a vector, zeros in the high half.
KVZ_FORCE_INLINE uint16x8_t hbd_load_1x4(const kvz_pixel *p)
{
  return vcombine_u16(vld1_u16(p), vdup_n_u16(0));
}

/*
 * SAD of a W-wide block. 8-pixel parts are loaded one row at a time, a
 * 4-pixel part as a row pair. Narrow blocks use separate accumulators for
 * even and odd rows to shorten the vaba dependency chains. Only the pixels
 * of the block are read.
 */
template <int W>
uint32_t reg_sad_hbd_w(const kvz_pixel *const data1, const kvz_pixel *const data2,
                       const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  static_assert(W > 0 && W % 4 == 0 && W <= 64, "unsupported width");
  constexpr int k8 = W / 8;                    // 8-pixel parts per row
  constexpr bool k4 = (W & 4) != 0;            // 4-pixel part at x = W - 4
  constexpr int kSplit = k8 <= 2 ? 2 : 1;      // accumulators per 8-pixel column
  constexpr int kAcc = k8 > 0 ? k8 * kSplit : 1;
  // Additions to a lane per row pair: 2 / kSplit for the 8-pixel parts, 1
  // for the 4-pixel part.
  constexpr int kRows = 2 * (kHbdSadAdds / (2 / kSplit));

  uint32x4_t total = vdupq_n_u32(0);
  int y = 0;
  while (y < height) {
    const int y_end = height - y < kRows ? height : y + kRows;
    uint16x8_t acc[kAcc];
    uint16x8_t acc4 = vdupq_n_u16(0);
    KVZ_UNROLL
    for (int i = 0; i < kAcc; ++i) {
      acc[i] = vdupq_n_u16(0);
    }
    for (; y + 2 <= y_end; y += 2) {
      const kvz_pixel *a = data1 + y * stride1;
      const kvz_pixel *b = data2 + y * stride2;
      KVZ_UNROLL
      for (int i = 0; i < k8; ++i) {
        acc[i * kSplit] = vabaq_u16(acc[i * kSplit], vld1q_u16(a + 8 * i), vld1q_u16(b + 8 * i));
        acc[i * kSplit + kSplit - 1] = vabaq_u16(acc[i * kSplit + kSplit - 1],
                                                 vld1q_u16(a + stride1 + 8 * i),
                                                 vld1q_u16(b + stride2 + 8 * i));
      }
      if constexpr (k4) {
        acc4 = vabaq_u16(acc4, hbd_load_2x4(a + W - 4, stride1), hbd_load_2x4(b + W - 4, stride2));
      }
    }
    if (y < y_end) {
      // Last row of an odd height.
      const kvz_pixel *a = data1 + y * stride1;
      const kvz_pixel *b = data2 + y * stride2;
      KVZ_UNROLL
      for (int i = 0; i < k8; ++i) {
        acc[i * kSplit] = vabaq_u16(acc[i * kSplit], vld1q_u16(a + 8 * i), vld1q_u16(b + 8 * i));
      }
      if constexpr (k4) {
        acc4 = vabaq_u16(acc4, hbd_load_1x4(a + W - 4), hbd_load_1x4(b + W - 4));
      }
      ++y;
    }
    if constexpr (k8 > 0) {
      KVZ_UNROLL
      for (int i = 0; i < kAcc; ++i) {
        total = vpadalq_u16(total, acc[i]);
      }
    }
    if constexpr (k4) {
      total = vpadalq_u16(total, acc4);
    }
  }
  return vaddvq_u32(total);
}

uint32_t reg_sad_hbd_w0(const kvz_pixel *const, const kvz_pixel *const,
                        const int32_t, const uint32_t, const uint32_t)
{
  return 0;
}

/*
 * SAD of a block of any width. With kHor, every row of data1 is compared
 * with the single pixel data2[y * stride2] (hor_sad). The absolute
 * differences are widened to 32 bits right away.
 */
template <bool kHor>
uint32_t sad_hbd_any(const kvz_pixel *data1, const kvz_pixel *data2,
                     int width, int height, unsigned stride1, unsigned stride2)
{
  uint32x4_t acc = vdupq_n_u32(0);
  uint32_t tail = 0;
  for (int y = 0; y < height; ++y) {
    const kvz_pixel *a = data1 + y * stride1;
    const kvz_pixel *b = data2 + y * stride2;
    const uint16x8_t ref8 = vdupq_n_u16(kHor ? b[0] : 0);
    int x = 0;
    for (; x + 8 <= width; x += 8) {
      acc = vpadalq_u16(acc, vabdq_u16(vld1q_u16(a + x), kHor ? ref8 : vld1q_u16(b + x)));
    }
    if (x + 4 <= width) {
      acc = vaddw_u16(acc, vabd_u16(vld1_u16(a + x), kHor ? vget_low_u16(ref8) : vld1_u16(b + x)));
      x += 4;
    }
    for (; x < width; ++x) {
      tail += abs((int)a[x] - (int)(kHor ? b[0] : b[x]));
    }
  }
  return vaddvq_u32(acc) + tail;
}

uint32_t reg_sad_hbd_neon(const kvz_pixel *const data1, const kvz_pixel *const data2,
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

optimized_sad_func_ptr_t get_optimized_sad_hbd_neon(int32_t width)
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
unsigned sad_nxn_hbd_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return reg_sad_hbd_w<N>(block1, block2, N, N, N) >> (KVZ_BIT_DEPTH - 8);
}

template <int N>
void sad_nxn_dual_hbd_neon(const pred_buffer preds, const kvz_pixel *const orig,
                           unsigned /* num_modes */, unsigned *costs_out)
{
  costs_out[0] = reg_sad_hbd_w<N>(preds[0], orig, N, N, N) >> (KVZ_BIT_DEPTH - 8);
  costs_out[1] = reg_sad_hbd_w<N>(preds[1], orig, N, N, N) >> (KVZ_BIT_DEPTH - 8);
}

// The reference row is the same for every row of the block.
uint32_t ver_sad_hbd_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                          int32_t width, int32_t height, uint32_t stride)
{
  return reg_sad_hbd_neon(pic_data, ref_data, width, height, stride, 0);
}

uint32_t hor_sad_hbd_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                          int32_t width, int32_t height, uint32_t pic_stride,
                          uint32_t ref_stride, uint32_t left, uint32_t right)
{
  uint32_t result = 0;
  if (left) {
    result += sad_hbd_any<true>(pic_data, ref_data + left, left, height, pic_stride, ref_stride);
    result += reg_sad_hbd_neon(pic_data + left, ref_data + left, width - left,
                               height, pic_stride, ref_stride);
  } else if (right) {
    result += reg_sad_hbd_neon(pic_data, ref_data, width - right,
                               height, pic_stride, ref_stride);
    result += sad_hbd_any<true>(pic_data + width - right, ref_data + width - right - 1,
                                right, height, pic_stride, ref_stride);
  } else {
    result += reg_sad_hbd_neon(pic_data, ref_data, width, height, pic_stride, ref_stride);
  }
  return result;
}

/*
 * SSD of a W x W block. The squares of the low and high half of a vector go
 * to the same 32-bit lanes, so a vector adds at most 2 * PIXEL_MAX^2 to a
 * lane; the lanes are moved to 64-bit lanes before they could overflow.
 */
constexpr uint32_t kHbdSsdAdds = 0xFFFFFFFFu / (2u * PIXEL_MAX * PIXEL_MAX);

KVZ_FORCE_INLINE uint32x4_t hbd_sqdiff_acc(uint32x4_t acc, uint16x8_t a, uint16x8_t b)
{
  const uint16x8_t d = vabdq_u16(a, b);
  acc = vmlal_u16(acc, vget_low_u16(d), vget_low_u16(d));
  return vmlal_high_u16(acc, d, d);
}

template <int W>
unsigned pixels_calc_ssd_hbd_w(const kvz_pixel *const ref, const kvz_pixel *const rec,
                               const int ref_stride, const int rec_stride)
{
  // Rows per vector and vectors per row.
  constexpr int kRowsPerVec = W >= 8 ? 1 : 2;
  constexpr int kVecsPerRow = W >= 8 ? W / 8 : 1;
  constexpr int kRowsPerFlush = (int)(kHbdSsdAdds / kVecsPerRow) * kRowsPerVec;
  static_assert(kRowsPerFlush >= kRowsPerVec, "SSD accumulation overflows");

  uint64x2_t acc64 = vdupq_n_u64(0);
  for (int y0 = 0; y0 < W; y0 += kRowsPerFlush) {
    const int y1 = W - y0 < kRowsPerFlush ? W : y0 + kRowsPerFlush;
    uint32x4_t acc = vdupq_n_u32(0);
    for (int y = y0; y < y1; y += kRowsPerVec) {
      const kvz_pixel *a = ref + y * ref_stride;
      const kvz_pixel *b = rec + y * rec_stride;
      if constexpr (W == 4) {
        acc = hbd_sqdiff_acc(acc, hbd_load_2x4(a, ref_stride), hbd_load_2x4(b, rec_stride));
      } else {
        KVZ_UNROLL
        for (int i = 0; i < kVecsPerRow; ++i) {
          acc = hbd_sqdiff_acc(acc, vld1q_u16(a + 8 * i), vld1q_u16(b + 8 * i));
        }
      }
    }
    acc64 = vpadalq_u32(acc64, acc);
  }
  const uint64_t ssd = vaddvq_u64(acc64);
  return (unsigned)(ssd >> (2 * (KVZ_BIT_DEPTH - 8)));
}

unsigned pixels_calc_ssd_hbd_neon(const kvz_pixel *const ref, const kvz_pixel *const rec,
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
 *   which is vrhadd and needs no clipping.
 * - Otherwise the sum is computed with a saturating 16-bit addition and
 *   rounded with vrshr, which does not overflow. When the sum saturates to
 *   INT16_MAX, the result is at least 32767 >> shift = PIXEL_MAX, which is
 *   also what clipping the exact result gives; a sum saturated to INT16_MIN
 *   gives a negative result, which clips to 0 like the exact one.
 */
constexpr int kBipredShift = 15 - KVZ_BIT_DEPTH;
static_assert((32767 >> kBipredShift) == PIXEL_MAX, "saturation must clip to PIXEL_MAX");

enum { BIPRED_PX_PX, BIPRED_IM_IM, BIPRED_PX_IM };

template <int MODE>
KVZ_FORCE_INLINE uint16x8_t bipred_hbd_combine(uint16x8_t a, uint16x8_t b)
{
  if constexpr (MODE == BIPRED_PX_PX) {
    return vrhaddq_u16(a, b);
  } else {
    if constexpr (MODE == BIPRED_PX_IM) {
      a = vshlq_n_u16(a, 14 - KVZ_BIT_DEPTH);
    }
    int16x8_t s = vqaddq_s16(vreinterpretq_s16_u16(a), vreinterpretq_s16_u16(b));
    s = vrshrq_n_s16(s, kBipredShift);
    s = vminq_s16(vmaxq_s16(s, vdupq_n_s16(0)), vdupq_n_s16(PIXEL_MAX));
    return vreinterpretq_u16_s16(s);
  }
}

template <int MODE>
KVZ_FORCE_INLINE uint16x4_t bipred_hbd_combine(uint16x4_t a, uint16x4_t b)
{
  if constexpr (MODE == BIPRED_PX_PX) {
    return vrhadd_u16(a, b);
  } else {
    if constexpr (MODE == BIPRED_PX_IM) {
      a = vshl_n_u16(a, 14 - KVZ_BIT_DEPTH);
    }
    int16x4_t s = vqadd_s16(vreinterpret_s16_u16(a), vreinterpret_s16_u16(b));
    s = vrshr_n_s16(s, kBipredShift);
    s = vmin_s16(vmax_s16(s, vdup_n_s16(0)), vdup_n_s16(PIXEL_MAX));
    return vreinterpret_u16_s16(s);
  }
}

// The generic formula for one pixel.
template <int MODE, typename T0, typename T1>
KVZ_FORCE_INLINE kvz_pixel bipred_hbd_scalar(T0 a, T1 b)
{
  const int32_t offset = 1 << (kBipredShift - 1);
  const int16_t s0 = MODE == BIPRED_IM_IM ? (int16_t)a : (int16_t)(a << (14 - KVZ_BIT_DEPTH));
  const int16_t s1 = MODE == BIPRED_PX_PX ? (int16_t)(b << (14 - KVZ_BIT_DEPTH)) : (int16_t)b;
  return (kvz_pixel)CLIP(0, PIXEL_MAX, (s0 + s1 + offset) >> kBipredShift);
}

/*
 * src0 and src1 are w x h blocks with a stride of w. The 4-pixel wide blocks
 * are processed two rows (8 contiguous source samples) at a time.
 */
template <int MODE, typename T0, typename T1>
void bipred_average_hbd_block(kvz_pixel *dst, const T0 *src0, const T1 *src1,
                              unsigned w, unsigned h, unsigned dst_stride)
{
  const uint16_t *s0 = (const uint16_t *)src0;
  const uint16_t *s1 = (const uint16_t *)src1;
  unsigned y = 0;
  if (w == 4) {
    for (; y + 2 <= h; y += 2) {
      const uint16x8_t r = bipred_hbd_combine<MODE>(vld1q_u16(s0 + 4 * y), vld1q_u16(s1 + 4 * y));
      vst1_u16(dst + y * dst_stride, vget_low_u16(r));
      vst1_u16(dst + (y + 1) * dst_stride, vget_high_u16(r));
    }
  }
  for (; y < h; ++y) {
    const uint16_t *a = s0 + y * w;
    const uint16_t *b = s1 + y * w;
    kvz_pixel *d = dst + y * dst_stride;
    unsigned x = 0;
    for (; x + 8 <= w; x += 8) {
      vst1q_u16(d + x, bipred_hbd_combine<MODE>(vld1q_u16(a + x), vld1q_u16(b + x)));
    }
    if (x + 4 <= w) {
      vst1_u16(d + x, bipred_hbd_combine<MODE>(vld1_u16(a + x), vld1_u16(b + x)));
      x += 4;
    }
    for (; x < w; ++x) {
      d[x] = bipred_hbd_scalar<MODE>(src0[y * w + x], src1[y * w + x]);
    }
  }
}

void bipred_average_hbd_neon(lcu_t *const lcu,
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
double pixel_var_hbd_neon(const kvz_pixel *arr, const uint32_t len)
{
  // A 32-bit lane of the squares takes at most 2 * PIXEL_MAX^2 per vector.
  constexpr uint32_t kBlock = 8 * (kHbdSsdAdds < 1024 ? kHbdSsdAdds : 1024);
  uint64x2_t sum64 = vdupq_n_u64(0);
  uint64x2_t sq64 = vdupq_n_u64(0);
  uint32_t i = 0;
  while (i + 8 <= len) {
    const uint32_t end = len - i < kBlock ? i + ((len - i) & ~7u) : i + kBlock;
    uint32x4_t sum32 = vdupq_n_u32(0);
    uint32x4_t sq32 = vdupq_n_u32(0);
    for (; i < end; i += 8) {
      const uint16x8_t v = vld1q_u16(arr + i);
      sum32 = vpadalq_u16(sum32, v);
      sq32 = vmlal_u16(sq32, vget_low_u16(v), vget_low_u16(v));
      sq32 = vmlal_high_u16(sq32, v, v);
    }
    sum64 = vpadalq_u32(sum64, sum32);
    sq64 = vpadalq_u32(sq64, sq32);
  }
  uint64_t sum = vaddvq_u64(sum64);
  uint64_t sumsq = vaddvq_u64(sq64);
  for (; i < len; ++i) {
    sum += arr[i];
    sumsq += (uint32_t)arr[i] * arr[i];
  }

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

} // namespace hbd_sad
} // namespace

#endif // COMPILE_ARM_NEON && KVZ_BIT_DEPTH > 8

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

namespace {

/*
 * 8-bit NEON versions of the functions of picture-generic.c. Every function
 * computes exactly the same result as the generic version, including the
 * block-splitting order of the any-size SATD functions.
 */

//////////////////////////////////////////////////////////////////////////
// SAD

/* Loads a 4-pixel row pair (rows y and y + 1) into one 64-bit vector. */
inline uint8x8_t load_row_pair4(const kvz_pixel *p, unsigned stride)
{
  uint32_t r0, r1;
  memcpy(&r0, p, 4);
  memcpy(&r1, p + stride, 4);
  return vreinterpret_u8_u64(vcreate_u64((uint64_t)r0 | ((uint64_t)r1 << 32)));
}

/*
 * SAD of a W-wide block of `height` rows. |a - b| is accumulated with vabal
 * into 16-bit lanes. A row adds at most ceil(W / 8) * 255 to a lane, so the
 * lanes are flushed into 32-bit lanes every kRows rows, which keeps the sums
 * exact for any height. Only the bytes of the block are read: a 4-wide strip
 * (W = 4 or 12) is loaded as row pairs.
 */
template <int W>
unsigned sad_rows(const kvz_pixel *p1, unsigned s1, const kvz_pixel *p2, unsigned s2, int height)
{
  static_assert(W == 4 || W == 8 || W == 12 || W == 24 || (W % 16 == 0 && W <= 64), "unsupported width");
  constexpr int kPerLane = (W + 7) / 8;
  constexpr int kRows = (65535 / (255 * kPerLane)) & ~1;
  constexpr int kStrip = W % 8 ? W - 4 : -1;          // x of the 4-wide strip
  constexpr int kEight = W % 16 >= 8 ? W - 8 - (kStrip >= 0 ? 4 : 0) : -1; // x of a lone 8-wide part
  constexpr int kSixteen = W - (W % 16);              // columns covered by 16-wide parts

  uint32x4_t total = vdupq_n_u32(0);
  unsigned tail = 0;
  for (int y0 = 0; y0 < height; y0 += kRows) {
    const int y1 = height - y0 < kRows ? height : y0 + kRows;
    uint16x8_t acc = vdupq_n_u16(0);
    if constexpr (kStrip >= 0) {
      int y = y0;
      for (; y + 2 <= y1; y += 2) {
        acc = vabal_u8(acc, load_row_pair4(p1 + y * s1 + kStrip, s1), load_row_pair4(p2 + y * s2 + kStrip, s2));
      }
      for (; y < y1; ++y) {
        for (int x = kStrip; x < W; ++x) {
          tail += abs((int)p1[y * s1 + x] - (int)p2[y * s2 + x]);
        }
      }
    }
    for (int y = y0; y < y1; ++y) {
      const kvz_pixel *r1 = p1 + y * s1;
      const kvz_pixel *r2 = p2 + y * s2;
      if constexpr (kEight >= 0) {
        acc = vabal_u8(acc, vld1_u8(r1 + kEight), vld1_u8(r2 + kEight));
      }
      for (int x = 0; x < kSixteen; x += 16) {
        const uint8x16_t a = vld1q_u8(r1 + x);
        const uint8x16_t b = vld1q_u8(r2 + x);
        acc = vabal_u8(acc, vget_low_u8(a), vget_low_u8(b));
        acc = vabal_high_u8(acc, a, b);
      }
    }
    total = vpadalq_u16(total, acc);
  }
  return vaddvq_u32(total) + tail;
}

/**
 * Sum of absolute differences of a width x height block. Common widths go
 * to the width-specialized kernel; others use the generic row loop, which
 * accumulates in 32-bit lanes.
 */
unsigned reg_sad_neon(const kvz_pixel *const data1, const kvz_pixel *const data2,
                      const int width, const int height,
                      const unsigned stride1, const unsigned stride2)
{
  switch (width) {
    case 4:  return sad_rows<4>(data1, stride1, data2, stride2, height);
    case 8:  return sad_rows<8>(data1, stride1, data2, stride2, height);
    case 12: return sad_rows<12>(data1, stride1, data2, stride2, height);
    case 16: return sad_rows<16>(data1, stride1, data2, stride2, height);
    case 24: return sad_rows<24>(data1, stride1, data2, stride2, height);
    case 32: return sad_rows<32>(data1, stride1, data2, stride2, height);
    case 48: return sad_rows<48>(data1, stride1, data2, stride2, height);
    case 64: return sad_rows<64>(data1, stride1, data2, stride2, height);
    default: break;
  }

  if (width < 8) {
    /* Too narrow for a vector to pay for its setup and reduction. */
    unsigned sum = 0;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        sum += abs((int)data1[y * stride1 + x] - (int)data2[y * stride2 + x]);
      }
    }
    return sum;
  }

  uint32x4_t acc = vdupq_n_u32(0);
  unsigned tail = 0;

  for (int y = 0; y < height; ++y) {
    const kvz_pixel *p1 = data1 + y * stride1;
    const kvz_pixel *p2 = data2 + y * stride2;
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      uint8x16_t d = vabdq_u8(vld1q_u8(p1 + x), vld1q_u8(p2 + x));
      acc = vpadalq_u16(acc, vpaddlq_u8(d));
    }
    if (x + 8 <= width) {
      uint8x8_t d = vabd_u8(vld1_u8(p1 + x), vld1_u8(p2 + x));
      acc = vaddw_u16(acc, vpaddl_u8(d));
      x += 8;
    }
    for (; x < width; ++x) {
      tail += abs((int)p1[x] - (int)p2[x]);
    }
  }
  return vaddvq_u32(acc) + tail;
}

template <int N>
unsigned sad_nxn_neon(const kvz_pixel *block1, const kvz_pixel *block2)
{
  return sad_rows<N>(block1, N, block2, N, N);
}

template <int N>
void sad_nxn_dual_neon(const pred_buffer preds, const kvz_pixel *const orig, unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = sad_rows<N>(preds[0], N, orig, N, N);
  costs_out[1] = sad_rows<N>(preds[1], N, orig, N, N);
}

/* Fixed-width SAD entry points for get_optimized_sad. */
template <int W>
uint32_t reg_sad_w_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                        const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return sad_rows<W>(data1, stride1, data2, stride2, height);
}

uint32_t reg_sad_w0_neon(const kvz_pixel * const data1, const kvz_pixel * const data2,
                         const int32_t height, const uint32_t stride1, const uint32_t stride2)
{
  return 0;
}

optimized_sad_func_ptr_t get_optimized_sad_neon(int32_t width)
{
  switch (width) {
    case 0:  return reg_sad_w0_neon;
    case 4:  return reg_sad_w_neon<4>;
    case 8:  return reg_sad_w_neon<8>;
    case 12: return reg_sad_w_neon<12>;
    case 16: return reg_sad_w_neon<16>;
    case 24: return reg_sad_w_neon<24>;
    case 32: return reg_sad_w_neon<32>;
    case 64: return reg_sad_w_neon<64>;
    default: return NULL;
  }
}

/* Vertically interpolated SAD: the reference row is shared by every row of
 * the picture block, so the second stride is zero. */
static uint32_t ver_sad_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                             int block_width, int block_height, unsigned pic_stride)
{
  return reg_sad_neon(pic_data, ref_data, block_width, block_height, pic_stride, 0);
}

/* SAD of a block against one reference sample per row, broadcast across the
 * row. */
static unsigned hor_sad_neon(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                             int block_width, int block_height, unsigned pic_stride, unsigned ref_stride)
{
  if (block_width < 8) {
    /* Narrow frame-border strips: a plain loop beats per-row vector setup. */
    unsigned sum = 0;
    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *p = pic_data + y * pic_stride;
      const int ref = ref_data[y * ref_stride];
      for (int x = 0; x < block_width; ++x) {
        sum += abs((int)p[x] - ref);
      }
    }
    return sum;
  }

  uint32x4_t acc = vdupq_n_u32(0);
  unsigned tail = 0;

  for (int y = 0; y < block_height; ++y) {
    const kvz_pixel *p = pic_data + y * pic_stride;
    const kvz_pixel ref = ref_data[y * ref_stride];
    const uint8x16_t vref16 = vdupq_n_u8(ref);
    const uint8x8_t vref8 = vdup_n_u8(ref);
    int x = 0;
    for (; x + 16 <= block_width; x += 16) {
      uint8x16_t d = vabdq_u8(vld1q_u8(p + x), vref16);
      acc = vpadalq_u16(acc, vpaddlq_u8(d));
    }
    if (x + 8 <= block_width) {
      uint8x8_t d = vabd_u8(vld1_u8(p + x), vref8);
      acc = vaddw_u16(acc, vpaddl_u8(d));
      x += 8;
    }
    for (; x < block_width; ++x) {
      tail += abs((int)p[x] - (int)ref);
    }
  }
  return vaddvq_u32(acc) + tail;
}

static uint32_t hor_sad_neon_wrapper(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                     int32_t width, int32_t height, uint32_t pic_stride,
                                     uint32_t ref_stride, uint32_t left, uint32_t right)
{
  uint32_t result = 0;
  if (left) {
    result += hor_sad_neon(pic_data, ref_data + left, left,
                           height, pic_stride, ref_stride);

    result += reg_sad_neon(pic_data + left, ref_data + left, width - left,
                           height, pic_stride, ref_stride);
  } else if (right) {
    result += reg_sad_neon(pic_data, ref_data, width - right,
                           height, pic_stride, ref_stride);

    result += hor_sad_neon(pic_data + width - right,
                           ref_data + width - right - 1,
                           right, height, pic_stride, ref_stride);
  } else {
    result += reg_sad_neon(pic_data, ref_data, width,
                           height, pic_stride, ref_stride);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////
// SATD

/* Loads 8 pixels of a and b and returns a - b as signed 16-bit lanes. */
static INLINE int16x8_t load_diff8(const kvz_pixel *a, const kvz_pixel *b)
{
  return vreinterpretq_s16_u16(vsubl_u8(vld1_u8(a), vld1_u8(b)));
}

/* Loads 4 pixels of a and b and returns a - b in the low four 16-bit lanes. */
static INLINE int16x4_t load_diff4(const kvz_pixel *a, const kvz_pixel *b)
{
  uint32_t va, vb;
  memcpy(&va, a, 4);
  memcpy(&vb, b, 4);
  uint8x8_t ua = vcreate_u8((uint64_t)va);
  uint8x8_t ub = vcreate_u8((uint64_t)vb);
  return vget_low_s16(vreinterpretq_s16_u16(vsubl_u8(ua, ub)));
}

/*
 * One pass of the 4-point Hadamard transform over four rows. Each vector holds
 * one row, so the transform runs down the columns. The result matches the
 * first two stages of hadamard_4x4_generic.
 */
static INLINE void hadamard4_pass(int16x4_t r[4])
{
  int16x4_t s03p = vadd_s16(r[0], r[3]);
  int16x4_t s03m = vsub_s16(r[0], r[3]);
  int16x4_t s12p = vadd_s16(r[1], r[2]);
  int16x4_t s12m = vsub_s16(r[1], r[2]);
  r[0] = vadd_s16(s03p, s12p);
  r[1] = vadd_s16(s03m, s12m);
  r[2] = vsub_s16(s03p, s12p);
  r[3] = vsub_s16(s03m, s12m);
}

static INLINE void transpose4x4_s16(int16x4_t r[4])
{
  int16x4x2_t t01 = vtrn_s16(r[0], r[1]);
  int16x4x2_t t23 = vtrn_s16(r[2], r[3]);
  int32x2x2_t c02 = vtrn_s32(vreinterpret_s32_s16(t01.val[0]), vreinterpret_s32_s16(t23.val[0]));
  int32x2x2_t c13 = vtrn_s32(vreinterpret_s32_s16(t01.val[1]), vreinterpret_s32_s16(t23.val[1]));
  r[0] = vreinterpret_s16_s32(c02.val[0]);
  r[1] = vreinterpret_s16_s32(c13.val[0]);
  r[2] = vreinterpret_s16_s32(c02.val[1]);
  r[3] = vreinterpret_s16_s32(c13.val[1]);
}

/*
 * One pass of the 8-point Hadamard transform over eight rows, the same
 * butterfly stages as the vertical pass of satd_8x8_subblock_generic.
 */
static INLINE void hadamard8_pass(int16x8_t m[8])
{
  int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  m[0] = vaddq_s16(b0, b1); m[1] = vsubq_s16(b0, b1);
  m[2] = vaddq_s16(b2, b3); m[3] = vsubq_s16(b2, b3);
  m[4] = vaddq_s16(b4, b5); m[5] = vsubq_s16(b4, b5);
  m[6] = vaddq_s16(b6, b7); m[7] = vsubq_s16(b6, b7);
}

/* 8x8 transpose of 16-bit lanes, built from 16-, 32- and 64-bit lane swaps. */
static INLINE void transpose8x8_s16(int16x8_t r[8])
{
  int16x8x2_t t01 = vtrnq_s16(r[0], r[1]);
  int16x8x2_t t23 = vtrnq_s16(r[2], r[3]);
  int16x8x2_t t45 = vtrnq_s16(r[4], r[5]);
  int16x8x2_t t67 = vtrnq_s16(r[6], r[7]);

  int32x4x2_t u02 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[0]), vreinterpretq_s32_s16(t23.val[0]));
  int32x4x2_t u13 = vtrnq_s32(vreinterpretq_s32_s16(t01.val[1]), vreinterpretq_s32_s16(t23.val[1]));
  int32x4x2_t u46 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[0]), vreinterpretq_s32_s16(t67.val[0]));
  int32x4x2_t u57 = vtrnq_s32(vreinterpretq_s32_s16(t45.val[1]), vreinterpretq_s32_s16(t67.val[1]));

  int16x8_t v00 = vreinterpretq_s16_s32(u02.val[0]);
  int16x8_t v01 = vreinterpretq_s16_s32(u02.val[1]);
  int16x8_t v10 = vreinterpretq_s16_s32(u13.val[0]);
  int16x8_t v11 = vreinterpretq_s16_s32(u13.val[1]);
  int16x8_t v20 = vreinterpretq_s16_s32(u46.val[0]);
  int16x8_t v21 = vreinterpretq_s16_s32(u46.val[1]);
  int16x8_t v30 = vreinterpretq_s16_s32(u57.val[0]);
  int16x8_t v31 = vreinterpretq_s16_s32(u57.val[1]);

  r[0] = vcombine_s16(vget_low_s16(v00), vget_low_s16(v20));
  r[4] = vcombine_s16(vget_high_s16(v00), vget_high_s16(v20));
  r[2] = vcombine_s16(vget_low_s16(v01), vget_low_s16(v21));
  r[6] = vcombine_s16(vget_high_s16(v01), vget_high_s16(v21));
  r[1] = vcombine_s16(vget_low_s16(v10), vget_low_s16(v30));
  r[5] = vcombine_s16(vget_high_s16(v10), vget_high_s16(v30));
  r[3] = vcombine_s16(vget_low_s16(v11), vget_low_s16(v31));
  r[7] = vcombine_s16(vget_high_s16(v11), vget_high_s16(v31));
}

/*
 * The last butterfly stage of the second Hadamard pass feeds straight into the
 * absolute sum, and |a + b| + |a - b| == 2 * max(|a|, |b|) for integers. So the
 * last stage and half of the absolute values are replaced by maxima. The sums
 * are exact: the 4x4 second-pass terms are at most 2040 and the 8x8 terms at
 * most 8160 (four of them, 32640, still fit a 16-bit lane).
 */
static INLINE uint32_t hadamard4_pass_absmax_sum(const int16x4_t r[4])
{
  const int16x4_t s03p = vadd_s16(r[0], r[3]);
  const int16x4_t s03m = vsub_s16(r[0], r[3]);
  const int16x4_t s12p = vadd_s16(r[1], r[2]);
  const int16x4_t s12m = vsub_s16(r[1], r[2]);
  const int16x4_t mp = vmax_s16(vabs_s16(s03p), vabs_s16(s12p));
  const int16x4_t mm = vmax_s16(vabs_s16(s03m), vabs_s16(s12m));
  return (uint32_t)vaddlv_s16(vadd_s16(mp, mm)); // half of the 4x4 absolute sum
}

static INLINE uint32_t hadamard8_pass_absmax_sum(const int16x8_t m[8])
{
  const int16x8_t a0 = vaddq_s16(m[0], m[4]), a4 = vsubq_s16(m[0], m[4]);
  const int16x8_t a1 = vaddq_s16(m[1], m[5]), a5 = vsubq_s16(m[1], m[5]);
  const int16x8_t a2 = vaddq_s16(m[2], m[6]), a6 = vsubq_s16(m[2], m[6]);
  const int16x8_t a3 = vaddq_s16(m[3], m[7]), a7 = vsubq_s16(m[3], m[7]);

  const int16x8_t b0 = vaddq_s16(a0, a2), b2 = vsubq_s16(a0, a2);
  const int16x8_t b1 = vaddq_s16(a1, a3), b3 = vsubq_s16(a1, a3);
  const int16x8_t b4 = vaddq_s16(a4, a6), b6 = vsubq_s16(a4, a6);
  const int16x8_t b5 = vaddq_s16(a5, a7), b7 = vsubq_s16(a5, a7);

  const int16x8_t m01 = vmaxq_s16(vabsq_s16(b0), vabsq_s16(b1));
  const int16x8_t m23 = vmaxq_s16(vabsq_s16(b2), vabsq_s16(b3));
  const int16x8_t m45 = vmaxq_s16(vabsq_s16(b4), vabsq_s16(b5));
  const int16x8_t m67 = vmaxq_s16(vabsq_s16(b6), vabsq_s16(b7));
  const uint16x8_t sum = vreinterpretq_u16_s16(vaddq_s16(vaddq_s16(m01, m23), vaddq_s16(m45, m67)));
  return vaddlvq_u16(sum); // half of the 8x8 absolute sum
}

/* SATD of one 4x4 block. Matches hadamard_4x4_generic, including the final
 * halving with rounding. */
static unsigned satd_4x4_subblock_neon(const kvz_pixel *buf1, int32_t stride1,
                                       const kvz_pixel *buf2, int32_t stride2)
{
  int16x4_t r[4];
  for (int i = 0; i < 4; ++i) {
    r[i] = load_diff4(buf1 + i * stride1, buf2 + i * stride2);
  }
  hadamard4_pass(r);
  transpose4x4_s16(r);
  // Generic returns (sum + 1) >> 1 of the absolute sum, which is 2 * half.
  return hadamard4_pass_absmax_sum(r);
}

/* SATD of one 8x8 block, matching satd_8x8_subblock_generic. */
static unsigned satd_8x8_subblock_neon(const kvz_pixel *org, int32_t stride_org,
                                       const kvz_pixel *cur, int32_t stride_cur)
{
  int16x8_t m[8];
  for (int i = 0; i < 8; ++i) {
    m[i] = load_diff8(org + i * stride_org, cur + i * stride_cur);
  }
  hadamard8_pass(m);
  transpose8x8_s16(m);
  // Generic returns (sum + 2) >> 2 of the absolute sum, with sum = 2 * half.
  return (2 * hadamard8_pass_absmax_sum(m) + 2) >> 2;
}

static unsigned satd_4x4_neon(const kvz_pixel *piOrg, const kvz_pixel *piCur)
{
  return satd_4x4_subblock_neon(piOrg, 4, piCur, 4);
}

static void satd_4x4_subblock_quad_neon(const kvz_pixel *preds[4], const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned costs[4])
{
  for (int k = 0; k < 4; ++k) {
    costs[k] = satd_4x4_subblock_neon(orig, orig_stride, preds[k], stride);
  }
}

static void satd_8x8_subblock_quad_neon(const kvz_pixel **preds, const int stride,
                                        const kvz_pixel *orig, const int orig_stride,
                                        unsigned *costs)
{
  for (int k = 0; k < 4; ++k) {
    costs[k] = satd_8x8_subblock_neon(orig, orig_stride, preds[k], stride);
  }
}

template <int N>
unsigned satd_nxn_neon(const kvz_pixel * const block1, const kvz_pixel * const block2)
{
  unsigned sum = 0;
  for (int y = 0; y < N; y += 8) {
    for (int x = 0; x < N; x += 8) {
      sum += satd_8x8_subblock_neon(&block1[y * N + x], N, &block2[y * N + x], N);
    }
  }
  return sum;
}

template <int N>
void satd_nxn_dual_neon(const pred_buffer preds, const kvz_pixel * const orig,
                        unsigned num_modes, unsigned *costs_out)
{
  for (int k = 0; k < 2; ++k) {
    costs_out[k] = satd_nxn_neon<N>(preds[k], orig);
  }
}

static void satd_4x4_dual_neon(const pred_buffer preds, const kvz_pixel * const orig,
                               unsigned num_modes, unsigned *costs_out)
{
  costs_out[0] = satd_4x4_neon(orig, preds[0]);
  costs_out[1] = satd_4x4_neon(orig, preds[1]);
}

static unsigned satd_any_size_neon(int width, int height,
                                   const kvz_pixel *block1, int stride1,
                                   const kvz_pixel *block2, int stride2)
{
  unsigned sum = 0;
  if (width % 8 != 0) {
    /* Process the first column using 4x4 blocks. */
    for (int y = 0; y < height; y += 4) {
      sum += satd_4x4_subblock_neon(&block1[y * stride1], stride1,
                                    &block2[y * stride2], stride2);
    }
    block1 += 4;
    block2 += 4;
    width -= 4;
  }
  if (height % 8 != 0) {
    /* Process the first row using 4x4 blocks. */
    for (int x = 0; x < width; x += 4) {
      sum += satd_4x4_subblock_neon(&block1[x], stride1,
                                    &block2[x], stride2);
    }
    block1 += 4 * stride1;
    block2 += 4 * stride2;
    height -= 4;
  }
  /* The rest can now be processed with 8x8 blocks. */
  for (int y = 0; y < height; y += 8) {
    const kvz_pixel *row1 = &block1[y * stride1];
    const kvz_pixel *row2 = &block2[y * stride2];
    for (int x = 0; x < width; x += 8) {
      sum += satd_8x8_subblock_neon(&row1[x], stride1, &row2[x], stride2);
    }
  }
  return sum;
}

/* Quad version of satd_any_size. It follows satd_any_size_quad_generic step by
 * step, so the results are identical for every input the generic accepts. */
static void satd_any_size_quad_neon(int width, int height,
                                    const kvz_pixel **preds, const int stride,
                                    const kvz_pixel *orig, const int orig_stride,
                                    unsigned num_modes, unsigned *costs_out, int8_t *valid)
{
  unsigned sums[4] = { 0 };
  const kvz_pixel *pred_ptrs[4] = { preds[0], preds[1], preds[2], preds[3] };
  const kvz_pixel *orig_ptr = orig;
  costs_out[0] = 0; costs_out[1] = 0; costs_out[2] = 0; costs_out[3] = 0;

  const int width_mod_8 = width % 8;
  if (width_mod_8 != 0) {
    const kvz_pixel *pred_ptrs_tmp[4] = { preds[0], preds[1], preds[2], preds[3] };
    /* Process the first column using 4x4 blocks. */
    for (int y = 0; y < height; y += 4) {
      satd_4x4_subblock_quad_neon(pred_ptrs_tmp, stride, &orig[y * orig_stride], orig_stride, sums);
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
    /* Process the first row using 4x4 blocks. */
    for (int x = 0; x < width; x += 4) {
      satd_4x4_subblock_quad_neon(pred_ptrs_tmp, stride, &orig_ptr[x], orig_stride, sums);
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
  /* The rest can now be processed with 8x8 blocks. */
  for (int y = height % 8; y < height; y += 8) {
    orig_ptr = &orig[y * orig_stride + width_mod_8];
    pred_ptrs[0] = &preds[0][y * stride + width_mod_8];
    pred_ptrs[1] = &preds[1][y * stride + width_mod_8];
    pred_ptrs[2] = &preds[2][y * stride + width_mod_8];
    pred_ptrs[3] = &preds[3][y * stride + width_mod_8];
    for (int x = width_mod_8; x < width; x += 8) {
      satd_8x8_subblock_quad_neon(pred_ptrs, stride, orig_ptr, orig_stride, sums);
      orig_ptr += 8;
      pred_ptrs[0] += 8;
      pred_ptrs[1] += 8;
      pred_ptrs[2] += 8;
      pred_ptrs[3] += 8;
      costs_out[0] += sums[0];
      costs_out[1] += sums[1];
      costs_out[2] += sums[2];
      costs_out[3] += sums[3];
    }
  }
}

//////////////////////////////////////////////////////////////////////////
// SSD

static unsigned pixels_calc_ssd_neon(const kvz_pixel *const ref, const kvz_pixel *const rec,
                                     const int ref_stride, const int rec_stride,
                                     const int width)
{
  /* |a - b| fits in 8 bits and its square in 16, so square the absolute
   * difference with a widening multiply and accumulate pairwise into 32 bits. */
  uint32x4_t acc = vdupq_n_u32(0);
  if (width == 4) {
    for (int y = 0; y < 4; y += 2) {
      const uint8x8_t d = vabd_u8(load_row_pair4(ref + y * ref_stride, ref_stride),
                                  load_row_pair4(rec + y * rec_stride, rec_stride));
      acc = vpadalq_u16(acc, vmull_u8(d, d));
    }
    return vaddvq_u32(acc);
  }
  unsigned tail = 0;
  for (int y = 0; y < width; ++y) {
    const kvz_pixel *r = ref + y * ref_stride;
    const kvz_pixel *c = rec + y * rec_stride;
    int x = 0;
    for (; x + 16 <= width; x += 16) {
      const uint8x16_t d = vabdq_u8(vld1q_u8(r + x), vld1q_u8(c + x));
      acc = vpadalq_u16(acc, vmull_u8(vget_low_u8(d), vget_low_u8(d)));
      acc = vpadalq_u16(acc, vmull_high_u8(d, d));
    }
    if (x + 8 <= width) {
      const uint8x8_t d = vabd_u8(vld1_u8(r + x), vld1_u8(c + x));
      acc = vpadalq_u16(acc, vmull_u8(d, d));
      x += 8;
    }
    for (; x < width; ++x) {
      const int diff = (int)r[x] - (int)c[x];
      tail += (unsigned)(diff * diff);
    }
  }
  return vaddvq_u32(acc) + tail;
}

//////////////////////////////////////////////////////////////////////////
// Bi-prediction

/* Rounds two 16-bit inputs with the bi-prediction offset and shift, and
 * saturates the result to 8 bits. Each input is widened to 32 bits first, so
 * the sum cannot overflow. */
static INLINE uint8x8_t bipred_round8(int16x8_t a, int16x8_t b)
{
  const int32x4_t off = vdupq_n_s32(1 << 6);
  int32x4_t lo = vaddq_s32(vaddl_s16(vget_low_s16(a), vget_low_s16(b)), off);
  int32x4_t hi = vaddq_s32(vaddl_high_s16(a, b), off);
  lo = vshrq_n_s32(lo, 7);
  hi = vshrq_n_s32(hi, 7);
  return vqmovun_s16(vcombine_s16(vqmovn_s32(lo), vqmovn_s32(hi)));
}

/* Converts 8-bit samples to the 14-bit intermediate precision. */
static INLINE int16x8_t px_to_im8(const kvz_pixel *p)
{
  return vreinterpretq_s16_u16(vshlq_n_u16(vmovl_u8(vld1_u8(p)), 6));
}

static void bipred_average_px_px(kvz_pixel *dst,
  kvz_pixel *px_L0,
  kvz_pixel *px_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = px_to_im8(px_L0 + y * pu_w + x);
      int16x8_t b = px_to_im8(px_L1 + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_L0 = px_L0[y * pu_w + x] << 6;
      int16_t sample_L1 = px_L1[y * pu_w + x] << 6;
      int32_t rounded = (sample_L0 + sample_L1 + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_im_im(kvz_pixel *dst,
  kvz_pixel_im *im_L0,
  kvz_pixel_im *im_L1,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = vld1q_s16(im_L0 + y * pu_w + x);
      int16x8_t b = vld1q_s16(im_L1 + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_L0 = im_L0[y * pu_w + x];
      int16_t sample_L1 = im_L1[y * pu_w + x];
      int32_t rounded = (sample_L0 + sample_L1 + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_px_im(kvz_pixel *dst,
  kvz_pixel *px,
  kvz_pixel_im *im,
  unsigned pu_w,
  unsigned pu_h,
  unsigned dst_stride)
{
  for (unsigned y = 0; y < pu_h; ++y) {
    unsigned x = 0;
    for (; x + 8 <= pu_w; x += 8) {
      int16x8_t a = px_to_im8(px + y * pu_w + x);
      int16x8_t b = vld1q_s16(im + y * pu_w + x);
      vst1_u8(dst + y * dst_stride + x, bipred_round8(a, b));
    }
    for (; x < pu_w; ++x) {
      int16_t sample_px = px[y * pu_w + x] << 6;
      int16_t sample_im = im[y * pu_w + x];
      int32_t rounded = (sample_px + sample_im + 64) >> 7;
      int32_t clipped = rounded < 0 ? 0 : (rounded > 255 ? 255 : rounded);
      dst[y * dst_stride + x] = (kvz_pixel)clipped;
    }
  }
}

static void bipred_average_neon(lcu_t *const lcu,
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
  // After reconstruction, merge the predictors by taking an average of each pixel.
  if (predict_luma) {
    unsigned pb_offset = SUB_SCU(pu_y) * LCU_WIDTH + SUB_SCU(pu_x);

    if (!(im_flags_L0 & 1) && !(im_flags_L1 & 1)) {
      bipred_average_px_px(lcu->rec.y + pb_offset, px_L0->y, px_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else if ((im_flags_L0 & 1) && (im_flags_L1 & 1)) {
      bipred_average_im_im(lcu->rec.y + pb_offset, im_L0->y, im_L1->y, pu_w, pu_h, LCU_WIDTH);
    } else {
      kvz_pixel    *src_px = (im_flags_L0 & 1) ? px_L1->y : px_L0->y;
      kvz_pixel_im *src_im = (im_flags_L0 & 1) ? im_L0->y : im_L1->y;
      bipred_average_px_im(lcu->rec.y + pb_offset, src_px, src_im, pu_w, pu_h, LCU_WIDTH);
    }
  }
  if (predict_chroma) {
    unsigned pb_offset = ((SUB_SCU(pu_y) >> chroma_shift_h) * (LCU_WIDTH >> chroma_shift_w)) + (SUB_SCU(pu_x) >> chroma_shift_w);
    unsigned pb_w = pu_w >> chroma_shift_w;
    unsigned pb_h = pu_h >> chroma_shift_h;

    if (!(im_flags_L0 & 2) && !(im_flags_L1 & 2)) {
      bipred_average_px_px(lcu->rec.u + pb_offset, px_L0->u, px_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_px(lcu->rec.v + pb_offset, px_L0->v, px_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    } else if ((im_flags_L0 & 2) && (im_flags_L1 & 2)) {
      bipred_average_im_im(lcu->rec.u + pb_offset, im_L0->u, im_L1->u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_im_im(lcu->rec.v + pb_offset, im_L0->v, im_L1->v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    } else {
      kvz_pixel    *src_px_u = (im_flags_L0 & 2) ? px_L1->u : px_L0->u;
      kvz_pixel_im *src_im_u = (im_flags_L0 & 2) ? im_L0->u : im_L1->u;
      kvz_pixel    *src_px_v = (im_flags_L0 & 2) ? px_L1->v : px_L0->v;
      kvz_pixel_im *src_im_v = (im_flags_L0 & 2) ? im_L0->v : im_L1->v;
      bipred_average_px_im(lcu->rec.u + pb_offset, src_px_u, src_im_u, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
      bipred_average_px_im(lcu->rec.v + pb_offset, src_px_v, src_im_v, pb_w, pb_h, LCU_WIDTH >> chroma_shift_w);
    }
  }
}

//////////////////////////////////////////////////////////////////////////
// Variance

/*
 * Variance of an array of 8-bit samples. The sums are exact integers, and the
 * variance is n * sum(x^2) - (sum x)^2 over n^2, which equals the two-pass
 * definition of pixel_var_generic.
 */
static double pixel_var_neon(const kvz_pixel *arr, const uint32_t len)
{
  /* The encoder calls this on whole frames (VAQ), so the sum of squares needs
   * 64 bits. Each 16-pixel step adds at most 4 * 255^2 to a 32-bit lane, so
   * the lanes are flushed into 64-bit accumulators every 4096 steps. */
  uint64x2_t sum64 = vdupq_n_u64(0);
  uint64x2_t sq64 = vdupq_n_u64(0);
  uint32_t i = 0;

  while (i + 16 <= len) {
    uint32x4_t sum_acc = vdupq_n_u32(0);
    uint32x4_t sq_acc = vdupq_n_u32(0);
    for (int step = 0; step < 4096 && i + 16 <= len; ++step, i += 16) {
      uint8x16_t v = vld1q_u8(arr + i);
      sum_acc = vpadalq_u16(sum_acc, vpaddlq_u8(v));
      sq_acc = vpadalq_u16(sq_acc, vmull_u8(vget_low_u8(v), vget_low_u8(v)));
      sq_acc = vpadalq_u16(sq_acc, vmull_high_u8(v, v));
    }
    sum64 = vpadalq_u32(sum64, sum_acc);
    sq64 = vpadalq_u32(sq64, sq_acc);
  }

  uint64_t sum = vaddvq_u64(sum64);
  uint64_t sumsq = vaddvq_u64(sq64);
  for (; i < len; ++i) {
    sum += arr[i];
    sumsq += (uint32_t)arr[i] * arr[i];
  }

  /* Population variance (n * sum(x^2) - sum(x)^2) / n^2, with the numerator
   * computed exactly. Like the AVX2 version this is not bit-identical to the
   * generic two-pass double loop, only more accurate. */
  const unsigned __int128 num = (unsigned __int128)len * sumsq - (unsigned __int128)sum * sum;
  const double n = (double)len;
  return (double)num / (n * n);
}

} // namespace

#endif // COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

int kvz_strategy_register_picture_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH > 8
  if (bitdepth == KVZ_BIT_DEPTH) {
    success &= kvz_strategyselector_register(opaque, "reg_sad", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::reg_sad_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "sad_4x4", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_hbd_neon<64>));
    success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::sad_nxn_dual_hbd_neon<64>));
    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::pixels_calc_ssd_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "bipred_average", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::bipred_average_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::get_optimized_sad_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "ver_sad", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::ver_sad_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "hor_sad", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::hor_sad_hbd_neon));
    success &= kvz_strategyselector_register(opaque, "pixel_var", "neon", 40, KVZ_STRATEGY_FPTR(hbd_sad::pixel_var_hbd_neon));
  }
#endif // COMPILE_ARM_NEON && KVZ_BIT_DEPTH > 8
#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "reg_sad", "neon", 40, KVZ_STRATEGY_FPTR(reg_sad_neon));

    success &= kvz_strategyselector_register(opaque, "sad_4x4", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_neon<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_neon<64>));

    success &= kvz_strategyselector_register(opaque, "sad_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<4>));
    success &= kvz_strategyselector_register(opaque, "sad_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<8>));
    success &= kvz_strategyselector_register(opaque, "sad_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<16>));
    success &= kvz_strategyselector_register(opaque, "sad_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<32>));
    success &= kvz_strategyselector_register(opaque, "sad_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(sad_nxn_dual_neon<64>));

    success &= kvz_strategyselector_register(opaque, "satd_4x4_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_4x4_dual_neon));
    success &= kvz_strategyselector_register(opaque, "satd_8x8_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<8>));
    success &= kvz_strategyselector_register(opaque, "satd_16x16_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<16>));
    success &= kvz_strategyselector_register(opaque, "satd_32x32_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<32>));
    success &= kvz_strategyselector_register(opaque, "satd_64x64_dual", "neon", 40, KVZ_STRATEGY_FPTR(satd_nxn_dual_neon<64>));
    success &= kvz_strategyselector_register(opaque, "satd_any_size", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_neon));
    success &= kvz_strategyselector_register(opaque, "satd_any_size_quad", "neon", 40, KVZ_STRATEGY_FPTR(satd_any_size_quad_neon));

    success &= kvz_strategyselector_register(opaque, "pixels_calc_ssd", "neon", 40, KVZ_STRATEGY_FPTR(pixels_calc_ssd_neon));
    success &= kvz_strategyselector_register(opaque, "bipred_average", "neon", 40, KVZ_STRATEGY_FPTR(bipred_average_neon));

    success &= kvz_strategyselector_register(opaque, "get_optimized_sad", "neon", 40, KVZ_STRATEGY_FPTR(get_optimized_sad_neon));
    success &= kvz_strategyselector_register(opaque, "ver_sad", "neon", 40, KVZ_STRATEGY_FPTR(ver_sad_neon));
    success &= kvz_strategyselector_register(opaque, "hor_sad", "neon", 40, KVZ_STRATEGY_FPTR(hor_sad_neon_wrapper));

    success &= kvz_strategyselector_register(opaque, "pixel_var", "neon", 40, KVZ_STRATEGY_FPTR(pixel_var_neon));
  }
#endif
  return success;
}
