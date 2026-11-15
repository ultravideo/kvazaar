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
 * Bit-exact with the generic implementation for 8-bit video. For 8-bit input
 * the bit-depth shifts and offsets of the generic code are zero, which lets the
 * distortion of a pixel be written as off*off - 2*(orig - rec)*off, where off is
 * zero for pixels that receive no offset, so no per-lane masking is needed.
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

#else

int kvz_strategy_register_sao_neon(void* opaque, uint8_t bitdepth)
{
  (void)opaque;
  (void)bitdepth;
  return 1;
}

#endif
