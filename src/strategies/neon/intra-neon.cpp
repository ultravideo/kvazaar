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

/**
 * \ingroup Optimization
 * \file
 * AArch64 NEON implementation of intra prediction (angular, planar, filtered DC).
 *
 * Bit-exact with the generic implementation for 8-bit video. The block width
 * is a template parameter, so every loop has a constant trip count and the
 * narrowing shifts that depend on the width can be immediates.
 */

#include "strategies/strategies-cpp.h"

extern "C" {
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "strategies/neon/intra-neon.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-intra.h"
#include "strategyselector.h"
}

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

namespace {

constexpr int8_t modedisp2sampledisp[9] = { 0, 2, 5, 9, 13, 17, 21, 26, 32 };
constexpr int16_t modedisp2invsampledisp[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / sampledisp

constexpr int log2_of(int w) { return w == 4 ? 2 : w == 8 ? 3 : w == 16 ? 4 : 5; }

/*
 * One row of an angular prediction with a fractional offset f (0 to 31):
 *   out[x] = ((32 - f) * r[x] + f * r[x + 1] + 16) >> 5
 * The 16-bit sum is at most (32 * 255) + 16, so widening u8 multiplies are
 * exact, and vrshrn by 5 is (sum + 16) >> 5. For W >= 8 only.
 */
template <int W>
inline void angular_filter_row(const kvz_pixel *r, int f, kvz_pixel *out)
{
  static_assert(W == 8 || W == 16 || W == 32, "angular_filter_row needs W >= 8");
  const uint8x8_t wa = vdup_n_u8((uint8_t)(32 - f));
  const uint8x8_t wb = vdup_n_u8((uint8_t)f);
  if constexpr (W == 8) {
    const uint16x8_t sum = vmlal_u8(vmull_u8(vld1_u8(r), wa), vld1_u8(r + 1), wb);
    vst1_u8(out, vrshrn_n_u16(sum, 5));
  } else {
    const uint8x16_t wa16 = vcombine_u8(wa, wa), wb16 = vcombine_u8(wb, wb);
    for (int x = 0; x < W; x += 16) {
      const uint8x16_t a = vld1q_u8(r + x);
      const uint8x16_t b = vld1q_u8(r + x + 1);
      const uint16x8_t lo = vmlal_u8(vmull_u8(vget_low_u8(a), wa), vget_low_u8(b), wb);
      const uint16x8_t hi = vmlal_high_u8(vmull_high_u8(a, wa16), b, wb16);
      vst1q_u8(out + x, vrshrn_high_n_u16(vrshrn_n_u16(lo, 5), hi, 5));
    }
  }
}

/* Two 4-wide angular rows at once: lanes 0-3 are row y (offset d0, fraction f0)
 * and lanes 4-7 are row y + 1. Only r[d .. d + 4] is read, like generic. A
 * fraction of 0 gives ((32 * r[x] + 16) >> 5) == r[x], the plain copy. */
inline uint8x8_t angular_rows_w4(const kvz_pixel *r, int d0, int f0, int d1, int f1)
{
  uint32_t a0, a1, b0, b1;
  memcpy(&a0, r + d0, 4);
  memcpy(&b0, r + d0 + 1, 4);
  memcpy(&a1, r + d1, 4);
  memcpy(&b1, r + d1 + 1, 4);
  const uint8x8_t a = vreinterpret_u8_u64(vcreate_u64((uint64_t)a0 | ((uint64_t)a1 << 32)));
  const uint8x8_t b = vreinterpret_u8_u64(vcreate_u64((uint64_t)b0 | ((uint64_t)b1 << 32)));
  const uint8x8_t wa = vreinterpret_u8_u32(vset_lane_u32((uint32_t)(32 - f1) * 0x01010101u,
                                                          vdup_n_u32((uint32_t)(32 - f0) * 0x01010101u), 1));
  const uint8x8_t wb = vreinterpret_u8_u32(vset_lane_u32((uint32_t)f1 * 0x01010101u,
                                                          vdup_n_u32((uint32_t)f0 * 0x01010101u), 1));
  return vrshrn_n_u16(vmlal_u8(vmull_u8(a, wa), b, wb), 5);
}

/* In-register transpose of an 8x8 byte block. */
inline void transpose8x8_u8(uint8x8_t r[8])
{
  const uint8x8x2_t t0 = vtrn_u8(r[0], r[1]);
  const uint8x8x2_t t1 = vtrn_u8(r[2], r[3]);
  const uint8x8x2_t t2 = vtrn_u8(r[4], r[5]);
  const uint8x8x2_t t3 = vtrn_u8(r[6], r[7]);
  const uint16x4x2_t u0 = vtrn_u16(vreinterpret_u16_u8(t0.val[0]), vreinterpret_u16_u8(t1.val[0]));
  const uint16x4x2_t u1 = vtrn_u16(vreinterpret_u16_u8(t0.val[1]), vreinterpret_u16_u8(t1.val[1]));
  const uint16x4x2_t u2 = vtrn_u16(vreinterpret_u16_u8(t2.val[0]), vreinterpret_u16_u8(t3.val[0]));
  const uint16x4x2_t u3 = vtrn_u16(vreinterpret_u16_u8(t2.val[1]), vreinterpret_u16_u8(t3.val[1]));
  const uint32x2x2_t v0 = vtrn_u32(vreinterpret_u32_u16(u0.val[0]), vreinterpret_u32_u16(u2.val[0]));
  const uint32x2x2_t v1 = vtrn_u32(vreinterpret_u32_u16(u1.val[0]), vreinterpret_u32_u16(u3.val[0]));
  const uint32x2x2_t v2 = vtrn_u32(vreinterpret_u32_u16(u0.val[1]), vreinterpret_u32_u16(u2.val[1]));
  const uint32x2x2_t v3 = vtrn_u32(vreinterpret_u32_u16(u1.val[1]), vreinterpret_u32_u16(u3.val[1]));
  r[0] = vreinterpret_u8_u32(v0.val[0]);
  r[4] = vreinterpret_u8_u32(v0.val[1]);
  r[1] = vreinterpret_u8_u32(v1.val[0]);
  r[5] = vreinterpret_u8_u32(v1.val[1]);
  r[2] = vreinterpret_u8_u32(v2.val[0]);
  r[6] = vreinterpret_u8_u32(v2.val[1]);
  r[3] = vreinterpret_u8_u32(v3.val[0]);
  r[7] = vreinterpret_u8_u32(v3.val[1]);
}

/* dst[y * W + x] = src[x * W + y] for a square block. */
template <int W>
inline void transpose_block(const kvz_pixel *src, kvz_pixel *dst)
{
  if constexpr (W == 4) {
    static const uint8_t idx[16] = { 0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15 };
    vst1q_u8(dst, vqtbl1q_u8(vld1q_u8(src), vld1q_u8(idx)));
  } else {
    for (int by = 0; by < W; by += 8) {
      for (int bx = 0; bx < W; bx += 8) {
        uint8x8_t r[8];
        for (int k = 0; k < 8; ++k) {
          r[k] = vld1_u8(src + (by + k) * W + bx);
        }
        transpose8x8_u8(r);
        for (int k = 0; k < 8; ++k) {
          vst1_u8(dst + (bx + k) * W + by, r[k]);
        }
      }
    }
  }
}

// Every row of a W x W block is the reference row (pure horizontal or vertical).
template <int W>
inline void fill_rows(const kvz_pixel *ref, kvz_pixel *buf)
{
  if constexpr (W == 4) {
    uint32_t r0;
    memcpy(&r0, ref, 4);
    vst1q_u8(buf, vreinterpretq_u8_u32(vdupq_n_u32(r0)));
  } else if constexpr (W == 8) {
    const uint8x8_t r0 = vld1_u8(ref);
    for (int y = 0; y < 8; ++y) vst1_u8(&buf[y * 8], r0);
  } else {
    const uint8x16_t r0 = vld1q_u8(ref);
    const uint8x16_t r1 = vld1q_u8(ref + (W == 32 ? 16 : 0));
    for (int y = 0; y < W; ++y) {
      vst1q_u8(&buf[y * W], r0);
      if constexpr (W == 32) vst1q_u8(&buf[y * W + 16], r1);
    }
  }
}

template <int W>
void angular_pred_w(const int_fast8_t intra_mode,
                    const kvz_pixel *const in_ref_above,
                    const kvz_pixel *const in_ref_left,
                    kvz_pixel *const dst)
{
  // Whether to swap references to always project on the left reference row.
  const bool vertical_mode = intra_mode >= 18;
  // Modes distance to horizontal or vertical mode.
  const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
  // Sample displacement per column in fractions of 32.
  const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

  // Temporary buffer for modes 11-25, indices from -W to W-1.
  kvz_pixel tmp_ref[2 * W];
  const kvz_pixel *ref_main;

  if (sample_disp < 0) {
    // Negative sample_disp means, we need to use both references.
    const kvz_pixel *ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

    for (int x = -1; x < W; ++x) {
      tmp_ref[x + W] = ref_main[x];
    }
    ref_main = &tmp_ref[W];

    // Extend the side reference to the negative indices of main reference.
    int_fast32_t col_sample_disp = 128; // rounding for the ">> 8"
    const int_fast16_t inv_abs_sample_disp = modedisp2invsampledisp[abs(mode_disp)];
    const int most_negative_index = (W * sample_disp) >> 5;
    for (int x = -2; x >= most_negative_index; --x) {
      col_sample_disp += inv_abs_sample_disp;
      const int side_index = (int)(col_sample_disp >> 8);
      tmp_ref[x + W] = ref_side[side_index - 1];
    }
  } else {
    ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
  }

  // Horizontal modes are computed in the vertical orientation and transposed.
  kvz_pixel tmp_out[W * W];
  kvz_pixel *buf = vertical_mode ? dst : tmp_out;

  if (sample_disp == 0) {
    fill_rows<W>(ref_main, buf);
  } else if constexpr (W == 4) {
    const int p0 = sample_disp, p1 = 2 * sample_disp, p2 = 3 * sample_disp, p3 = 4 * sample_disp;
    vst1_u8(&buf[0], angular_rows_w4(ref_main, p0 >> 5, p0 & 31, p1 >> 5, p1 & 31));
    vst1_u8(&buf[8], angular_rows_w4(ref_main, p2 >> 5, p2 & 31, p3 >> 5, p3 & 31));
  } else {
    int delta_pos = 0;
    for (int y = 0; y < W; ++y) {
      delta_pos += sample_disp;
      const kvz_pixel *r = ref_main + (delta_pos >> 5);
      const int frac = delta_pos & 31;
      // A fraction of 0 is a plain copy; it must not read r[W], which generic
      // does not read either (it is past the reference array for modes 2 and 34).
      if (frac) {
        angular_filter_row<W>(r, frac, &buf[y * W]);
      } else if constexpr (W == 8) {
        vst1_u8(&buf[y * W], vld1_u8(r));
      } else {
        for (int x = 0; x < W; x += 16) vst1q_u8(&buf[y * W + x], vld1q_u8(r + x));
      }
    }
  }

  if (!vertical_mode) {
    transpose_block<W>(tmp_out, dst);
  }
}

/*
 * Planar: out = ((W-1-x) * left[y] + (x+1) * top_right + (W-1-y) * top[x]
 *                + (y+1) * bottom_left + W) >> (log2(W) + 1).
 * All terms are non-negative and the sum is at most 4 * (W - 1/2) * 255 + W,
 * 32162 for W = 32, so 16-bit unsigned lanes are exact. The shift is an
 * immediate because W is a template parameter.
 */
template <int W>
void planar_w(const kvz_pixel *const ref_top, const kvz_pixel *const ref_left, kvz_pixel *const dst)
{
  constexpr int kShift = log2_of(W) + 1;
  const uint16x8_t vtr = vdupq_n_u16(ref_top[W + 1]);
  const uint16_t bottom_left = ref_left[W + 1];

  if constexpr (W == 4) {
    // Two rows per vector: lanes 0-3 are row y, lanes 4-7 are row y + 1.
    static const uint16_t a_tab[8] = { 3, 2, 1, 0, 3, 2, 1, 0 };   // W - 1 - x
    static const uint16_t b_tab[8] = { 1, 2, 3, 4, 1, 2, 3, 4 };   // x + 1
    const uint16x8_t a = vld1q_u16(a_tab), b = vld1q_u16(b_tab);
    uint32_t t4;
    memcpy(&t4, ref_top + 1, 4);
    const uint16x8_t top = vmovl_u8(vreinterpret_u8_u32(vdup_n_u32(t4)));
    const uint16x8_t base = vmlaq_u16(vdupq_n_u16(W), b, vtr);   // (x + 1) * top_right + W
    for (int y = 0; y < 4; y += 2) {
      const uint16x8_t left = vcombine_u16(vdup_n_u16(ref_left[y + 1]), vdup_n_u16(ref_left[y + 2]));
      const uint16x8_t wy = vcombine_u16(vdup_n_u16(W - 1 - y), vdup_n_u16(W - 2 - y));
      const uint16x8_t bl = vcombine_u16(vdup_n_u16((y + 1) * bottom_left), vdup_n_u16((y + 2) * bottom_left));
      uint16x8_t acc = vaddq_u16(base, bl);
      acc = vmlaq_u16(acc, a, left);
      acc = vmlaq_u16(acc, wy, top);
      vst1_u8(&dst[y * 4], vshrn_n_u16(acc, kShift));
    }
  } else {
    // Per-column vectors: a = W - 1 - x, top[x], and (x + 1) * top_right + W.
    uint16x8_t a[W / 8], top[W / 8], base[W / 8];
    for (int c = 0; c < W / 8; ++c) {
      static const uint16_t iota[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
      const uint16x8_t x = vaddq_u16(vld1q_u16(iota), vdupq_n_u16((uint16_t)(8 * c)));
      a[c] = vsubq_u16(vdupq_n_u16(W - 1), x);
      top[c] = vmovl_u8(vld1_u8(ref_top + 1 + 8 * c));
      base[c] = vmlaq_u16(vdupq_n_u16(W), vaddq_u16(x, vdupq_n_u16(1)), vtr);
    }
    for (int y = 0; y < W; ++y) {
      const uint16x8_t left = vdupq_n_u16(ref_left[y + 1]);
      const uint16x8_t wy = vdupq_n_u16((uint16_t)(W - 1 - y));
      const uint16x8_t bl = vdupq_n_u16((uint16_t)((y + 1) * bottom_left));
      for (int c = 0; c < W / 8; ++c) {
        uint16x8_t acc = vaddq_u16(base[c], bl);
        acc = vmlaq_u16(acc, a[c], left);
        acc = vmlaq_u16(acc, wy, top[c]);
        vst1_u8(&dst[y * W + 8 * c], vshrn_n_u16(acc, kShift));
      }
    }
  }
}

template <int W>
void filtered_dc_w(const kvz_pixel *const ref_top, const kvz_pixel *const ref_left, kvz_pixel *const out_block)
{
  constexpr int kShift = log2_of(W) + 1;
  const kvz_pixel *const top = ref_top + 1;
  const kvz_pixel *const left = ref_left + 1;

  // DC value: the rounded mean of the top and left references.
  uint32_t sum;
  if constexpr (W == 4) {
    uint32_t t, l;
    memcpy(&t, top, 4);
    memcpy(&l, left, 4);
    sum = vaddlv_u8(vreinterpret_u8_u64(vcreate_u64((uint64_t)t | ((uint64_t)l << 32))));
  } else {
    uint16x8_t acc = vdupq_n_u16(0);
    for (int i = 0; i < W; i += 8) {
      acc = vaddw_u8(acc, vld1_u8(top + i));
      acc = vaddw_u8(acc, vld1_u8(left + i));
    }
    sum = vaddlvq_u16(acc);
  }
  const kvz_pixel dc_val = (kvz_pixel)((sum + W) >> kShift);

  // Fill the block with the DC value. W * W is a multiple of 16.
  const uint8x16_t fill = vdupq_n_u8(dc_val);
  for (int i = 0; i < W * W; i += 16) {
    vst1q_u8(out_block + i, fill);
  }

  // Top row: (top[x] + 3 * dc + 2) >> 2, in 16-bit lanes.
  const uint16x8_t bias = vdupq_n_u16((uint16_t)(3 * dc_val + 2));
  if constexpr (W == 4) {
    uint32_t t;
    memcpy(&t, top, 4);
    const uint8x8_t r = vshrn_n_u16(vaddw_u8(bias, vreinterpret_u8_u32(vdup_n_u32(t))), 2);
    vst1_lane_u32(reinterpret_cast<uint32_t *>(out_block), vreinterpret_u32_u8(r), 0);
  } else {
    for (int x = 0; x < W; x += 8) {
      vst1_u8(out_block + x, vshrn_n_u16(vaddw_u8(bias, vld1_u8(top + x)), 2));
    }
  }
  // Left column, same filter.
  for (int y = 1; y < W; ++y) {
    out_block[y * W] = (kvz_pixel)((left[y] + 3 * dc_val + 2) >> 2);
  }
  // Top-left corner with ([1 2 1] / 4).
  out_block[0] = (kvz_pixel)((left[0] + 2 * dc_val + top[0] + 2) >> 2);
}

void angular_pred_neon(const int_fast8_t log2_width, const int_fast8_t intra_mode,
                       const kvz_pixel *const in_ref_above, const kvz_pixel *const in_ref_left,
                       kvz_pixel *const dst)
{
  switch (log2_width) {
    case 2: angular_pred_w<4>(intra_mode, in_ref_above, in_ref_left, dst); break;
    case 3: angular_pred_w<8>(intra_mode, in_ref_above, in_ref_left, dst); break;
    case 4: angular_pred_w<16>(intra_mode, in_ref_above, in_ref_left, dst); break;
    default: angular_pred_w<32>(intra_mode, in_ref_above, in_ref_left, dst); break;
  }
}

void intra_pred_planar_neon(const int_fast8_t log2_width, const kvz_pixel *const ref_top,
                            const kvz_pixel *const ref_left, kvz_pixel *const dst)
{
  switch (log2_width) {
    case 2: planar_w<4>(ref_top, ref_left, dst); break;
    case 3: planar_w<8>(ref_top, ref_left, dst); break;
    case 4: planar_w<16>(ref_top, ref_left, dst); break;
    default: planar_w<32>(ref_top, ref_left, dst); break;
  }
}

void intra_pred_filtered_dc_neon(const int_fast8_t log2_width, const kvz_pixel *const ref_top,
                                 const kvz_pixel *const ref_left, kvz_pixel *const out_block)
{
  switch (log2_width) {
    case 2: filtered_dc_w<4>(ref_top, ref_left, out_block); break;
    case 3: filtered_dc_w<8>(ref_top, ref_left, out_block); break;
    case 4: filtered_dc_w<16>(ref_top, ref_left, out_block); break;
    default: filtered_dc_w<32>(ref_top, ref_left, out_block); break;
  }
}

} // namespace

int kvz_strategy_register_intra_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "angular_pred", "neon", 40, KVZ_STRATEGY_FPTR(angular_pred_neon));
    success &= kvz_strategyselector_register(opaque, "intra_pred_planar", "neon", 40, KVZ_STRATEGY_FPTR(intra_pred_planar_neon));
    success &= kvz_strategyselector_register(opaque, "intra_pred_filtered_dc", "neon", 40, KVZ_STRATEGY_FPTR(intra_pred_filtered_dc_neon));
  }

  return success;
}

#else

int kvz_strategy_register_intra_neon(void* opaque, uint8_t bitdepth)
{
  (void)opaque;
  (void)bitdepth;
  return 1;
}

#endif
