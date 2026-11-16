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
 * AArch64 NEON implementation of sample adaptive offset (SAO) functions.
 *
 * Bit-exact with the generic implementation for 8, 10 and 12-bit video. For
 * 8-bit input the bit-depth shifts and offsets of the generic code are zero,
 * which lets the distortion of a pixel be written as off*off - 2*(orig - rec)*off,
 * where off is zero for pixels that receive no offset, so no per-lane masking
 * is needed. The 10 and 12-bit versions work on 16-bit lanes.
 */

#include "global.h" // IWYU pragma: keep
#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "sao.h"
#include "strategies/neon/sao-neon.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-sao.h"
#include "strategyselector.h"

#if COMPILE_ARM_NEON && KVZ_BIT_DEPTH == 8

// Maps edge_idx (0..4) to the SAO edge category, as in sao_shared_generics.h.
static const int8_t sao_eo_idx_to_category[5] = { 1, 2, 0, 3, 4 };

// Sign of (x - y) in {-1, 0, 1} for each lane.
static inline int8x8_t sign3_u8x8(uint8x8_t x, uint8x8_t y)
{
  int8x8_t gt = vreinterpret_s8_u8(vcgt_u8(x, y)); // -1 where x > y
  int8x8_t lt = vreinterpret_s8_u8(vclt_u8(x, y)); // -1 where x < y
  return vsub_s8(lt, gt);                          // (x > y) - (x < y)
}

// edge_idx = 2 + SIGN3(c - a) + SIGN3(c - b), in 0..4 for each lane.
static inline uint8x8_t eo_index_u8x8(uint8x8_t c, uint8x8_t a, uint8x8_t b)
{
  int8x8_t s = vadd_s8(sign3_u8x8(c, a), sign3_u8x8(c, b));
  return vreinterpret_u8_s8(vadd_s8(s, vdup_n_s8(2)));
}

// Accumulates off*off - 2*(orig - ref)*off over the lanes, which is the
// sum-of-squared-error change for an offset `off` (zero means no change).
static inline int32x4_t ddistortion_acc(int32x4_t acc, uint8x8_t orig, uint8x8_t ref, int8x8_t off)
{
  int16x8_t d = vsubq_s16(vreinterpretq_s16_u16(vmovl_u8(orig)),
                          vreinterpretq_s16_u16(vmovl_u8(ref)));
  int16x8_t o = vmovl_s8(off);
  int16x8_t term = vsubq_s16(vmulq_s16(o, o), vshlq_n_s16(vmulq_s16(d, o), 1));
  return vpadalq_s16(acc, term);
}

static int sao_edge_ddistortion_neon(const encoder_control_t* const encoder,
                                     const kvz_pixel *orig_data,
                                     const kvz_pixel *rec_data,
                                     int block_width,
                                     int block_height,
                                     int eo_class,
                                     const int offsets[NUM_SAO_EDGE_CATEGORIES])
{
  (void)encoder;
  const vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  const vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

  // Offset per edge_idx.
  int8_t tbl[16] = { 0 };
  for (int k = 0; k < 5; ++k) {
    tbl[k] = (int8_t)offsets[sao_eo_idx_to_category[k]];
  }
  const int8x16_t vt = vld1q_s8(tbl);

  int32x4_t acc = vdupq_n_s32(0);
  int32_t sum = 0;

  for (int y = 1; y < block_height - 1; ++y) {
    const kvz_pixel *c_row = rec_data + y * block_width;
    const kvz_pixel *a_row = rec_data + (y + a_ofs.y) * block_width + a_ofs.x;
    const kvz_pixel *b_row = rec_data + (y + b_ofs.y) * block_width + b_ofs.x;
    const kvz_pixel *o_row = orig_data + y * block_width;

    int x = 1;
    for (; x + 8 <= block_width - 1; x += 8) {
      uint8x8_t c = vld1_u8(c_row + x);
      uint8x8_t a = vld1_u8(a_row + x);
      uint8x8_t b = vld1_u8(b_row + x);
      uint8x8_t o = vld1_u8(o_row + x);
      int8x8_t off = vqtbl1_s8(vt, eo_index_u8x8(c, a, b));
      acc = ddistortion_acc(acc, o, c, off);
    }
    for (; x < block_width - 1; ++x) {
      int c = c_row[x];
      int a = a_row[x];
      int b = b_row[x];
      int o = o_row[x];
      int idx = 2 + SIGN3(c - a) + SIGN3(c - b);
      int off = tbl[idx];
      sum += off * off - 2 * (o - c) * off;
    }
  }

  return sum + vaddvq_s32(acc);
}

static void calc_sao_edge_dir_neon(const encoder_control_t* const encoder,
                                   const kvz_pixel *orig_data,
                                   const kvz_pixel *rec_data,
                                   int eo_class,
                                   int block_width,
                                   int block_height,
                                   int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  (void)encoder;
  const vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  const vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

  // Per edge_idx: sum of (orig - rec) and count, kept as signed lanes.
  int32x4_t sum_acc[5];
  int32x4_t cnt_acc[5];
  for (int k = 0; k < 5; ++k) {
    sum_acc[k] = vdupq_n_s32(0);
    cnt_acc[k] = vdupq_n_s32(0);
  }
  int32_t sum_s[5] = { 0, 0, 0, 0, 0 };
  int32_t cnt_s[5] = { 0, 0, 0, 0, 0 };

  for (int y = 1; y < block_height - 1; ++y) {
    const kvz_pixel *c_row = rec_data + y * block_width;
    const kvz_pixel *a_row = rec_data + (y + a_ofs.y) * block_width + a_ofs.x;
    const kvz_pixel *b_row = rec_data + (y + b_ofs.y) * block_width + b_ofs.x;
    const kvz_pixel *o_row = orig_data + y * block_width;

    int x = 1;
    for (; x + 8 <= block_width - 1; x += 8) {
      uint8x8_t c = vld1_u8(c_row + x);
      uint8x8_t a = vld1_u8(a_row + x);
      uint8x8_t b = vld1_u8(b_row + x);
      uint8x8_t o = vld1_u8(o_row + x);
      uint8x8_t idx = eo_index_u8x8(c, a, b);
      int16x8_t d = vsubq_s16(vreinterpretq_s16_u16(vmovl_u8(o)),
                              vreinterpretq_s16_u16(vmovl_u8(c)));
      for (int k = 0; k < 5; ++k) {
        // Lanes equal to k become -1 (all ones), others 0.
        int16x8_t m = vmovl_s8(vreinterpret_s8_u8(vceq_u8(idx, vdup_n_u8((uint8_t)k))));
        sum_acc[k] = vpadalq_s16(sum_acc[k], vandq_s16(d, m));
        cnt_acc[k] = vpadalq_s16(cnt_acc[k], m);
      }
    }
    for (; x < block_width - 1; ++x) {
      int c = c_row[x];
      int a = a_row[x];
      int b = b_row[x];
      int idx = 2 + SIGN3(c - a) + SIGN3(c - b);
      sum_s[idx] += o_row[x] - c;
      cnt_s[idx] += 1;
    }
  }

  for (int k = 0; k < 5; ++k) {
    const int cat = sao_eo_idx_to_category[k];
    cat_sum_cnt[0][cat] += vaddvq_s32(sum_acc[k]) + sum_s[k];
    cat_sum_cnt[1][cat] += -vaddvq_s32(cnt_acc[k]) + cnt_s[k];
  }
}

static void sao_reconstruct_color_neon(const encoder_control_t * const encoder,
                                       const kvz_pixel *rec_data,
                                       kvz_pixel *new_rec_data,
                                       const sao_info_t *sao,
                                       int stride,
                                       int new_stride,
                                       int block_width,
                                       int block_height,
                                       color_t color_i)
{
  (void)encoder;
  const int offset_v = color_i == COLOR_V ? 5 : 0;

  if (sao->type == SAO_TYPE_BAND) {
    // 8-bit: band index is the top 5 bits, i.e. sample >> 3.
    const int band_pos = color_i == COLOR_V ? 1 : 0;
    const int cur_bp = sao->band_position[band_pos];

    int8_t tbl[16] = { 0 };
    for (int k = 0; k < 4; ++k) {
      tbl[k] = (int8_t)sao->offsets[k + 1 + 5 * band_pos];
    }
    const int8x16_t vt = vld1q_s8(tbl);

    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *rec = rec_data + y * stride;
      kvz_pixel *out = new_rec_data + y * new_stride;

      int x = 0;
      for (; x + 8 <= block_width; x += 8) {
        uint8x8_t r = vld1_u8(rec + x);
        uint8x8_t t = vsub_u8(vshr_n_u8(r, 3), vdup_n_u8((uint8_t)cur_bp));
        // Lanes whose band is one of the four offset bands use table entry t,
        // all others use entry 4, which is zero.
        uint8x8_t valid = vclt_u8(t, vdup_n_u8(4));
        uint8x8_t idx = vbsl_u8(valid, t, vdup_n_u8(4));
        int8x8_t off = vqtbl1_s8(vt, idx);
        int16x8_t v = vaddq_s16(vreinterpretq_s16_u16(vmovl_u8(r)), vmovl_s8(off));
        vst1_u8(out + x, vqmovun_s16(v));
      }
      for (; x < block_width; ++x) {
        int val = rec[x];
        int band = (val >> 3) - cur_bp;
        int off = 0;
        if (band >= 0 && band <= 3) {
          off = sao->offsets[band + 1 + 5 * band_pos];
        }
        int res = val + off;
        out[x] = (kvz_pixel)(res < 0 ? 0 : (res > 255 ? 255 : res));
      }
    }
  }
  else {
    const vector2d_t a_ofs = g_sao_edge_offsets[sao->eo_class][0];
    const vector2d_t b_ofs = g_sao_edge_offsets[sao->eo_class][1];

    int8_t tbl[16] = { 0 };
    for (int k = 0; k < 5; ++k) {
      tbl[k] = (int8_t)sao->offsets[sao_eo_idx_to_category[k] + offset_v];
    }
    const int8x16_t vt = vld1q_s8(tbl);

    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *c_row = rec_data + y * stride;
      const kvz_pixel *a_row = rec_data + (y + a_ofs.y) * stride + a_ofs.x;
      const kvz_pixel *b_row = rec_data + (y + b_ofs.y) * stride + b_ofs.x;
      kvz_pixel *out = new_rec_data + y * new_stride;

      int x = 0;
      for (; x + 8 <= block_width; x += 8) {
        uint8x8_t c = vld1_u8(c_row + x);
        uint8x8_t a = vld1_u8(a_row + x);
        uint8x8_t b = vld1_u8(b_row + x);
        int8x8_t off = vqtbl1_s8(vt, eo_index_u8x8(c, a, b));
        int16x8_t v = vaddq_s16(vreinterpretq_s16_u16(vmovl_u8(c)), vmovl_s8(off));
        vst1_u8(out + x, vqmovun_s16(v));
      }
      for (; x < block_width; ++x) {
        int c = c_row[x];
        int a = a_row[x];
        int b = b_row[x];
        int idx = 2 + SIGN3(c - a) + SIGN3(c - b);
        int res = c + tbl[idx];
        out[x] = (kvz_pixel)(res < 0 ? 0 : (res > 255 ? 255 : res));
      }
    }
  }
}

static int sao_band_ddistortion_neon(const encoder_state_t * const state,
                                     const kvz_pixel *orig_data,
                                     const kvz_pixel *rec_data,
                                     int block_width,
                                     int block_height,
                                     int band_pos,
                                     const int sao_bands[4])
{
  (void)state;

  // 8-bit: band index is rec >> 3.
  int8_t tbl[16] = { 0 };
  for (int k = 0; k < 4; ++k) {
    tbl[k] = (int8_t)sao_bands[k];
  }
  const int8x16_t vt = vld1q_s8(tbl);

  int32x4_t acc = vdupq_n_s32(0);
  int32_t sum = 0;

  for (int y = 0; y < block_height; ++y) {
    const kvz_pixel *rec = rec_data + y * block_width;
    const kvz_pixel *orig = orig_data + y * block_width;

    int x = 0;
    for (; x + 8 <= block_width; x += 8) {
      uint8x8_t r = vld1_u8(rec + x);
      uint8x8_t o = vld1_u8(orig + x);
      uint8x8_t t = vsub_u8(vshr_n_u8(r, 3), vdup_n_u8((uint8_t)band_pos));
      uint8x8_t valid = vclt_u8(t, vdup_n_u8(4));
      uint8x8_t idx = vbsl_u8(valid, t, vdup_n_u8(4));
      int8x8_t off = vqtbl1_s8(vt, idx);
      acc = ddistortion_acc(acc, o, r, off);
    }
    for (; x < block_width; ++x) {
      int r = rec[x];
      int o = orig[x];
      int band = (r >> 3) - band_pos;
      int off = 0;
      if (band >= 0 && band <= 3) {
        off = sao_bands[band];
      }
      sum += off * off - 2 * (o - r) * off;
    }
  }

  return sum + vaddvq_s32(acc);
}

int kvz_strategy_register_sao_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "sao_edge_ddistortion", "neon", 40, &sao_edge_ddistortion_neon);
    success &= kvz_strategyselector_register(opaque, "calc_sao_edge_dir", "neon", 40, &calc_sao_edge_dir_neon);
    success &= kvz_strategyselector_register(opaque, "sao_reconstruct_color", "neon", 40, &sao_reconstruct_color_neon);
    success &= kvz_strategyselector_register(opaque, "sao_band_ddistortion", "neon", 40, &sao_band_ddistortion_neon);
  }

  return success;
}

#elif COMPILE_ARM_NEON && KVZ_BIT_DEPTH > 8

#include <assert.h>
#include <stdint.h>

#include "strategies/generic/sao_shared_generics.h"

/*
 * 10 and 12-bit SAO. kvz_pixel is uint16_t. The samples (at most 4095) and
 * their differences fit int16 lanes, so 8 samples are processed at once.
 *
 * The SAO categories are looked up from 16-bit tables of 8 entries with
 * vqtbl1q, indexed with edge_idx = 2 + SIGN3(c - a) + SIGN3(c - b) (0..4) for
 * edge offsets and with the band minus the band position (0..3, or 4 for
 * other bands) for band offsets. The entries past those are zero, and index 7
 * marks lanes outside the block.
 *
 * The distortion change of a sample with offset `off` is
 *   (diff - off)^2 - diff^2 == off * (off - 2 * diff),
 * which is zero for off == 0 like the masked generic code, and is summed in
 * 32-bit lanes. The sums wrap like the 32-bit generic sums.
 */

// Maps edge_idx to the SAO edge category, as in sao_shared_generics.h.
static const int sao_eo_idx_to_category[5] = { 1, 2, 0, 3, 4 };

// The distortion functions use the generic code for offsets larger than
// this, which never occur (SAO_ABS_OFFSET_MAX is 31). With |diff| <= 4095,
// off - 2 * diff fits int16.
#define SAO_HBD_MAX_SIMD_OFFSET 4096

static inline int16_t sat_int16(int v)
{
  return (int16_t)(v < INT16_MIN ? INT16_MIN : (v > INT16_MAX ? INT16_MAX : v));
}

static inline int clip_pixel_int(int v)
{
  return v < 0 ? 0 : (v > PIXEL_MAX ? PIXEL_MAX : v);
}

// tab[idx] for each 16-bit lane, idx in 0..7. tab holds 8 int16 entries.
static inline int16x8_t lookup_table16(uint8x16_t tab, uint16x8_t idx)
{
  // Byte indices 2 * idx and 2 * idx + 1.
  const uint16x8_t ctrl = vmlaq_n_u16(vdupq_n_u16(0x0100), idx, 0x0202);
  return vreinterpretq_s16_u8(vqtbl1q_u8(tab, vreinterpretq_u8_u16(ctrl)));
}

// edge_idx = 2 + SIGN3(c - a) + SIGN3(c - b), modulo 2^16 (the compares
// give all ones for true).
static inline uint16x8_t edge_idx_u16(uint16x8_t a, uint16x8_t b, uint16x8_t c)
{
  const uint16x8_t sa = vsubq_u16(vcgtq_u16(a, c), vcgtq_u16(c, a));
  const uint16x8_t sb = vsubq_u16(vcgtq_u16(b, c), vcgtq_u16(c, b));
  return vaddq_u16(vaddq_u16(sa, sb), vdupq_n_u16(2));
}

// Band index minus band_pos, or 4 for samples outside the four bands.
static inline uint16x8_t band_idx_u16(uint16x8_t rec, uint16x8_t band_pos)
{
  // Negative values are large as unsigned.
  return vminq_u16(vsubq_u16(vshrq_n_u16(rec, KVZ_BIT_DEPTH - 5), band_pos), vdupq_n_u16(4));
}

// Sets the index of the lanes outside the block (mask 0) to 7.
static inline uint16x8_t mask_idx(uint16x8_t idx, uint16x8_t valid)
{
  return vorrq_u16(idx, vbicq_u16(vdupq_n_u16(7), valid));
}

// Mask of the last `rem` lanes (1..7) of a vector.
static inline uint16x8_t last_lanes_mask(int rem)
{
  static const uint16_t iota[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  return vcgtq_u16(vld1q_u16(iota), vdupq_n_u16((uint16_t)(7 - rem)));
}

static inline int sign3_int(int v)
{
  return (v > 0) - (v < 0);
}

// acc += off * (off - 2 * diff) over the lanes.
static inline int32x4_t ddist_acc(int32x4_t acc, int16x8_t off, int16x8_t diff)
{
  const int16x8_t t = vsubq_s16(off, vaddq_s16(diff, diff));
  acc = vmlal_s16(acc, vget_low_s16(off), vget_low_s16(t));
  return vmlal_high_s16(acc, off, t);
}

// orig - c
static inline int16x8_t edge_diff(uint16x8_t orig, uint16x8_t c)
{
  return vreinterpretq_s16_u16(vsubq_u16(orig, c));
}

static int sao_edge_ddistortion_neon(const encoder_control_t* const encoder,
                                     const kvz_pixel *orig_data,
                                     const kvz_pixel *rec_data,
                                     int block_width,
                                     int block_height,
                                     int eo_class,
                                     const int offsets[NUM_SAO_EDGE_CATEGORIES])
{
  assert(encoder->bitdepth == KVZ_BIT_DEPTH);
  int16_t tab[8] = { 0 };
  for (int k = 0; k < 5; ++k) {
    const int off = offsets[sao_eo_idx_to_category[k]];
    if (off < -SAO_HBD_MAX_SIMD_OFFSET || off > SAO_HBD_MAX_SIMD_OFFSET) {
      return sao_edge_ddistortion_generic(encoder, orig_data, rec_data, block_width,
                                          block_height, eo_class, offsets);
    }
    tab[k] = (int16_t)off;
  }
  const uint8x16_t vtab = vreinterpretq_u8_s16(vld1q_s16(tab));
  const vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  const vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

  // Samples 1 .. block_width - 2 of each row.
  const int x_end = block_width - 1;
  const int inner = block_width - 2;
  const int rem = inner > 0 ? inner & 7 : 0;
  const uint16x8_t tail = last_lanes_mask(rem ? rem : 1);

  int32x4_t acc = vdupq_n_s32(0);
  int32_t sum = 0;

  for (int y = 1; y < block_height - 1; ++y) {
    const kvz_pixel *c_row = rec_data + y * block_width;
    const kvz_pixel *a_row = rec_data + (y + a_ofs.y) * block_width + a_ofs.x;
    const kvz_pixel *b_row = rec_data + (y + b_ofs.y) * block_width + b_ofs.x;
    const kvz_pixel *o_row = orig_data + y * block_width;

    int x = 1;
    for (; x + 8 <= x_end; x += 8) {
      const uint16x8_t c = vld1q_u16(c_row + x);
      const int16x8_t off = lookup_table16(vtab, edge_idx_u16(vld1q_u16(a_row + x), vld1q_u16(b_row + x), c));
      acc = ddist_acc(acc, off, edge_diff(vld1q_u16(o_row + x), c));
    }
    if (rem && inner >= 8) {
      // The last 8 samples of the row, of which the last rem are new.
      const int p = x_end - 8;
      const uint16x8_t c = vld1q_u16(c_row + p);
      const uint16x8_t idx = mask_idx(edge_idx_u16(vld1q_u16(a_row + p), vld1q_u16(b_row + p), c), tail);
      acc = ddist_acc(acc, lookup_table16(vtab, idx), edge_diff(vld1q_u16(o_row + p), c));
    } else {
      for (; x < x_end; ++x) {
        const int c = c_row[x];
        const int off = tab[2 + sign3_int(c - a_row[x]) + sign3_int(c - b_row[x])];
        const int diff = o_row[x] - c;
        sum += off * (off - 2 * diff);
      }
    }
  }

  return sum + vaddvq_s32(acc);
}

static void calc_sao_edge_dir_neon(const encoder_control_t* const encoder,
                                   const kvz_pixel *orig_data,
                                   const kvz_pixel *rec_data,
                                   int eo_class,
                                   int block_width,
                                   int block_height,
                                   int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  assert(encoder->bitdepth == KVZ_BIT_DEPTH);
  (void)encoder;
  const vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  const vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

  const int x_end = block_width - 1;
  const int inner = block_width - 2;
  const int rem = inner > 0 ? inner & 7 : 0;
  const uint16x8_t all = vdupq_n_u16(0xFFFF);
  const uint16x8_t tail = last_lanes_mask(rem ? rem : 1);

  // Per edge_idx: sum of the differences and the negated count.
  int32x4_t sum_acc[5];
  int32x4_t cnt_acc[5];
  for (int k = 0; k < 5; ++k) {
    sum_acc[k] = vdupq_n_s32(0);
    cnt_acc[k] = vdupq_n_s32(0);
  }
  int32_t sum_s[5] = { 0, 0, 0, 0, 0 };
  int32_t cnt_s[5] = { 0, 0, 0, 0, 0 };

  for (int y = 1; y < block_height - 1; ++y) {
    const kvz_pixel *c_row = rec_data + y * block_width;
    const kvz_pixel *a_row = rec_data + (y + a_ofs.y) * block_width + a_ofs.x;
    const kvz_pixel *b_row = rec_data + (y + b_ofs.y) * block_width + b_ofs.x;
    const kvz_pixel *o_row = orig_data + y * block_width;

    int x = 1;
    while (x + 8 <= x_end || (rem && inner >= 8 && x < x_end)) {
      uint16x8_t valid;
      int p;
      if (x + 8 <= x_end) {
        valid = all;
        p = x;
        x += 8;
      } else {
        // The last 8 samples of the row, of which the last rem are new.
        valid = tail;
        p = x_end - 8;
        x = x_end;
      }
      const uint16x8_t c = vld1q_u16(c_row + p);
      const uint16x8_t idx = mask_idx(edge_idx_u16(vld1q_u16(a_row + p), vld1q_u16(b_row + p), c), valid);
      const int16x8_t d = edge_diff(vld1q_u16(o_row + p), c);
      for (int k = 0; k < 5; ++k) {
        // Lanes equal to k become -1 (all ones), others 0.
        const int16x8_t m = vreinterpretq_s16_u16(vceqq_u16(idx, vdupq_n_u16((uint16_t)k)));
        sum_acc[k] = vpadalq_s16(sum_acc[k], vandq_s16(d, m));
        cnt_acc[k] = vpadalq_s16(cnt_acc[k], m);
      }
    }
    for (; x < x_end; ++x) {
      const int c = c_row[x];
      const int idx = 2 + sign3_int(c - a_row[x]) + sign3_int(c - b_row[x]);
      sum_s[idx] += o_row[x] - c;
      cnt_s[idx] += 1;
    }
  }

  for (int k = 0; k < 5; ++k) {
    const int cat = sao_eo_idx_to_category[k];
    cat_sum_cnt[0][cat] += vaddvq_s32(sum_acc[k]) + sum_s[k];
    cat_sum_cnt[1][cat] += -vaddvq_s32(cnt_acc[k]) + cnt_s[k];
  }
}

static int sao_band_ddistortion_neon(const encoder_state_t * const state,
                                     const kvz_pixel *orig_data,
                                     const kvz_pixel *rec_data,
                                     int block_width,
                                     int block_height,
                                     int band_pos,
                                     const int sao_bands[4])
{
  assert(state->encoder_control->bitdepth == KVZ_BIT_DEPTH);
  int16_t tab[8] = { 0 };
  for (int k = 0; k < 4; ++k) {
    if (sao_bands[k] < -SAO_HBD_MAX_SIMD_OFFSET || sao_bands[k] > SAO_HBD_MAX_SIMD_OFFSET) {
      return sao_band_ddistortion_generic(state, orig_data, rec_data, block_width,
                                          block_height, band_pos, sao_bands);
    }
    tab[k] = (int16_t)sao_bands[k];
  }
  const uint8x16_t vtab = vreinterpretq_u8_s16(vld1q_s16(tab));
  const uint16x8_t vbp = vdupq_n_u16((uint16_t)band_pos);
  const int rem = block_width & 7;
  const uint16x8_t tail = last_lanes_mask(rem ? rem : 1);

  int32x4_t acc = vdupq_n_s32(0);
  int32_t sum = 0;

  for (int y = 0; y < block_height; ++y) {
    const kvz_pixel *rec = rec_data + y * block_width;
    const kvz_pixel *orig = orig_data + y * block_width;
    int x = 0;
    for (; x + 8 <= block_width; x += 8) {
      const uint16x8_t r = vld1q_u16(rec + x);
      const int16x8_t off = lookup_table16(vtab, band_idx_u16(r, vbp));
      acc = ddist_acc(acc, off, vreinterpretq_s16_u16(vsubq_u16(vld1q_u16(orig + x), r)));
    }
    if (rem && block_width >= 8) {
      // The last 8 samples of the row, of which the last rem are new.
      const int p = block_width - 8;
      const uint16x8_t r = vld1q_u16(rec + p);
      const int16x8_t off = lookup_table16(vtab, mask_idx(band_idx_u16(r, vbp), tail));
      acc = ddist_acc(acc, off, vreinterpretq_s16_u16(vsubq_u16(vld1q_u16(orig + p), r)));
    } else {
      for (; x < block_width; ++x) {
        const int band = (rec[x] >> (KVZ_BIT_DEPTH - 5)) - band_pos;
        const int off = band >= 0 && band <= 3 ? tab[band] : 0;
        sum += off * (off - 2 * (orig[x] - rec[x]));
      }
    }
  }

  return sum + vaddvq_s32(acc);
}

// clip(rec + off) to [0, PIXEL_MAX]. The offsets are saturated to int16,
// which does not change the result for samples in 0 .. PIXEL_MAX.
static inline uint16x8_t add_offset_clip(uint16x8_t rec, int16x8_t off)
{
  const int16x8_t v = vqaddq_s16(vreinterpretq_s16_u16(rec), off);
  return vreinterpretq_u16_s16(vminq_s16(vmaxq_s16(v, vdupq_n_s16(0)), vdupq_n_s16(PIXEL_MAX)));
}

static void sao_reconstruct_color_neon(const encoder_control_t * const encoder,
                                       const kvz_pixel *rec_data,
                                       kvz_pixel *new_rec_data,
                                       const sao_info_t *sao,
                                       int stride,
                                       int new_stride,
                                       int block_width,
                                       int block_height,
                                       color_t color_i)
{
  assert(encoder->bitdepth == KVZ_BIT_DEPTH);
  (void)encoder;
  int16_t tab[8] = { 0 };

  // The rows are processed 8 and 4 samples at a time and the last 0..3
  // samples one at a time, so nothing outside the block is read or written
  // (except the neighbours of edge offsets, like in the generic code).
  if (sao->type == SAO_TYPE_BAND) {
    const int band_pos_idx = color_i == COLOR_V ? 1 : 0;
    const int band_pos = sao->band_position[band_pos_idx];
    for (int k = 0; k < 4; ++k) {
      tab[k] = sat_int16(sao->offsets[k + 1 + 5 * band_pos_idx]);
    }
    const uint8x16_t vtab = vreinterpretq_u8_s16(vld1q_s16(tab));
    const uint16x8_t vbp = vdupq_n_u16((uint16_t)band_pos);

    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *rec = rec_data + y * stride;
      kvz_pixel *out = new_rec_data + y * new_stride;
      int x = 0;
      for (; x + 8 <= block_width; x += 8) {
        const uint16x8_t r = vld1q_u16(rec + x);
        vst1q_u16(out + x, add_offset_clip(r, lookup_table16(vtab, band_idx_u16(r, vbp))));
      }
      if (x + 4 <= block_width) {
        const uint16x8_t r = vcombine_u16(vld1_u16(rec + x), vdup_n_u16(0));
        vst1_u16(out + x, vget_low_u16(add_offset_clip(r, lookup_table16(vtab, band_idx_u16(r, vbp)))));
        x += 4;
      }
      for (; x < block_width; ++x) {
        const int band = (rec[x] >> (KVZ_BIT_DEPTH - 5)) - band_pos;
        const int off = band >= 0 && band <= 3 ? tab[band] : 0;
        out[x] = (kvz_pixel)clip_pixel_int(rec[x] + off);
      }
    }
  } else {
    const int offset_v = color_i == COLOR_V ? 5 : 0;
    for (int k = 0; k < 5; ++k) {
      tab[k] = sat_int16(sao->offsets[sao_eo_idx_to_category[k] + offset_v]);
    }
    const uint8x16_t vtab = vreinterpretq_u8_s16(vld1q_s16(tab));
    const vector2d_t a_ofs = g_sao_edge_offsets[sao->eo_class][0];
    const vector2d_t b_ofs = g_sao_edge_offsets[sao->eo_class][1];

    for (int y = 0; y < block_height; ++y) {
      const kvz_pixel *c_row = rec_data + y * stride;
      const kvz_pixel *a_row = c_row + a_ofs.y * stride + a_ofs.x;
      const kvz_pixel *b_row = c_row + b_ofs.y * stride + b_ofs.x;
      kvz_pixel *out = new_rec_data + y * new_stride;
      int x = 0;
      for (; x + 8 <= block_width; x += 8) {
        const uint16x8_t c = vld1q_u16(c_row + x);
        const uint16x8_t idx = edge_idx_u16(vld1q_u16(a_row + x), vld1q_u16(b_row + x), c);
        vst1q_u16(out + x, add_offset_clip(c, lookup_table16(vtab, idx)));
      }
      if (x + 4 <= block_width) {
        const uint16x4_t zero = vdup_n_u16(0);
        const uint16x8_t c = vcombine_u16(vld1_u16(c_row + x), zero);
        const uint16x8_t idx = edge_idx_u16(vcombine_u16(vld1_u16(a_row + x), zero),
                                            vcombine_u16(vld1_u16(b_row + x), zero), c);
        vst1_u16(out + x, vget_low_u16(add_offset_clip(c, lookup_table16(vtab, idx))));
        x += 4;
      }
      for (; x < block_width; ++x) {
        const int c = c_row[x];
        const int off = tab[2 + sign3_int(c - a_row[x]) + sign3_int(c - b_row[x])];
        out[x] = (kvz_pixel)clip_pixel_int(c + off);
      }
    }
  }
}

int kvz_strategy_register_sao_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  if (bitdepth == KVZ_BIT_DEPTH) {
    success &= kvz_strategyselector_register(opaque, "sao_edge_ddistortion", "neon", 40, &sao_edge_ddistortion_neon);
    success &= kvz_strategyselector_register(opaque, "calc_sao_edge_dir", "neon", 40, &calc_sao_edge_dir_neon);
    success &= kvz_strategyselector_register(opaque, "sao_reconstruct_color", "neon", 40, &sao_reconstruct_color_neon);
    success &= kvz_strategyselector_register(opaque, "sao_band_ddistortion", "neon", 40, &sao_band_ddistortion_neon);
  }

  return success;
}

#else

int kvz_strategy_register_sao_neon(void* opaque, uint8_t bitdepth)
{
  (void)opaque;
  (void)bitdepth;
  return 1;
}

#endif
