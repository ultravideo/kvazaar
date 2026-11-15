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

#include "strategies/neon/quant-neon.h"

#include "global.h" // IWYU pragma: keep

#if COMPILE_ARM_NEON

#include <stdlib.h>

#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "scalinglist.h"
#include "strategies/generic/quant-generic.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-quant.h"
#include "strategyselector.h"
#include "tables.h"
#include "transform.h"

#define QUANT_SHIFT 14
#define SCAN_SET_SIZE 16
#define LOG2_SCAN_SET_SIZE 4

/*
 * NEON versions of the element-wise loops of the quantization and
 * dequantization. The products are computed with the same widths as the
 * generic code: the quantizer multiplies 64-bit values, so it uses vmull_s32
 * to get 64-bit products. The sign-hiding search depends on the coefficients
 * in a data-dependent order, so it stays scalar and is the same code as the
 * generic version.
 */

/* Quantization of four coefficients: returns the unsigned levels in 32-bit
 * lanes, and the 64-bit (abs * quant + add) >> q_bits values are kept for the
 * sign-hiding pass. */
static void kvz_quant_neon(const encoder_state_t * const state, const coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t * const scan = kvz_g_sig_last_scan[scan_idx][log2_block_size - 1];

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6, encoder->cfg.chroma_format == KVZ_CSP_420);
  const uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
  const int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  const coeff_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
  const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  const int32_t q_bits = QUANT_SHIFT + qp_scaled / 6 + transform_shift;
  const int32_t add = ((state->frame->slicetype == KVZ_SLICE_I) ? 171 : 85) << (q_bits - 9);
  const int32_t q_bits8 = q_bits - 8;
  const int32_t n_total = width * height;

  uint32_t ac_sum = 0;
  int32_t n = 0;

  /* 32-bit lanes are exact: |coef| <= 32768 and the scaling factor is an int16,
   * so |abs * q| <= 2^30, and adding add (< 2^26 for 8-bit) still fits in int32.
   * vabsq_s16(-32768) wraps to -32768, which as u16 is the correct 32768. */
  const int32x4_t vadd = vdupq_n_s32(add);
  const int32x4_t vneg_qbits = vdupq_n_s32(-q_bits);
  int32x4_t vacc = vdupq_n_s32(0);
  for (; n + 8 <= n_total; n += 8) {
    const int16x8_t c = vld1q_s16(coef + n);
    const uint16x8_t a16 = vreinterpretq_u16_s16(vabsq_s16(c));
    const int16x8_t q = vld1q_s16(quant_coeff + n);
    int32x4_t ll = vmlaq_s32(vadd, vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(a16))), vmovl_s16(vget_low_s16(q)));
    int32x4_t lh = vmlaq_s32(vadd, vreinterpretq_s32_u32(vmovl_high_u16(a16)), vmovl_high_s16(q));
    ll = vshlq_s32(ll, vneg_qbits);
    lh = vshlq_s32(lh, vneg_qbits);
    vacc = vaddq_s32(vacc, vaddq_s32(ll, lh));

    // level * sign, with sign = -1 for negative coefficients (and +1 for zero).
    const int16x8_t neg = vshrq_n_s16(c, 15);
    const int32x4_t nl = vmovl_s16(vget_low_s16(neg));
    const int32x4_t nh = vmovl_high_s16(neg);
    ll = vsubq_s32(veorq_s32(ll, nl), nl);
    lh = vsubq_s32(veorq_s32(lh, nh), nh);
    vst1q_s16(q_coef + n, vqmovn_high_s32(vqmovn_s32(ll), lh));
  }
  ac_sum += (uint32_t)vaddvq_s32(vacc);
  for (; n < n_total; ++n) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t sign = (level < 0 ? -1 : 1);

    int32_t curr_quant_coeff = quant_coeff[n];
    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    ac_sum += level;

    level *= sign;
    q_coef[n] = (coeff_t)(CLIP(-32768, 32767, level));
  }

  if (!encoder->cfg.signhide_enable || ac_sum < 2) return;

  int32_t delta_u[LCU_WIDTH*LCU_WIDTH >> 2];

  n = 0;
  const int32x4_t vqbits = vdupq_n_s32(q_bits);
  const int32x4_t vneg_qbits8 = vdupq_n_s32(-q_bits8);
  for (; n + 8 <= n_total; n += 8) {
    const int16x8_t c = vld1q_s16(coef + n);
    const uint16x8_t a16 = vreinterpretq_u16_s16(vabsq_s16(c));
    const int16x8_t q = vld1q_s16(quant_coeff + n);
    const int32x4_t pl = vmulq_s32(vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(a16))), vmovl_s16(vget_low_s16(q)));
    const int32x4_t ph = vmulq_s32(vreinterpretq_s32_u32(vmovl_high_u16(a16)), vmovl_high_s16(q));
    const int32x4_t ll = vshlq_s32(vaddq_s32(pl, vadd), vneg_qbits);
    const int32x4_t lh = vshlq_s32(vaddq_s32(ph, vadd), vneg_qbits);
    // delta_u = (abs * q - (level << q_bits)) >> q_bits8; every term fits in int32.
    vst1q_s32(delta_u + n, vshlq_s32(vsubq_s32(pl, vshlq_s32(ll, vqbits)), vneg_qbits8));
    vst1q_s32(delta_u + n + 4, vshlq_s32(vsubq_s32(ph, vshlq_s32(lh, vqbits)), vneg_qbits8));
  }
  for (; n < n_total; ++n) {
    int32_t level = coef[n];
    int64_t abs_level = (int64_t)abs(level);
    int32_t curr_quant_coeff = quant_coeff[n];

    level = (abs_level * curr_quant_coeff + add) >> q_bits;
    delta_u[n] = (int32_t)((abs_level * curr_quant_coeff - ((int64_t)level << q_bits)) >> q_bits8);
  }

  // Sign data hiding: sequential, same as the generic code.
  {
    const uint32_t * const scan_ = scan;
    int32_t n2, last_cg = -1, abssum = 0, subset, subpos;
    for (subset = (width*height - 1) >> LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg = -1;
      subpos = subset << LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n2 = SCAN_SET_SIZE - 1; n2 >= 0; n2--)  {
        if (q_coef[scan_[n2 + subpos]])  {
          last_nz_pos_in_cg = n2;
          break;
        }
      }

      // First coeff pos
      for (n2 = 0; n2 < SCAN_SET_SIZE; n2++) {
        if (q_coef[scan_[n2 + subpos]]) {
          first_nz_pos_in_cg = n2;
          break;
        }
      }

      // Sum all kvz_quant coeffs between first and last
      for (n2 = first_nz_pos_in_cg; n2 <= last_nz_pos_in_cg; n2++) {
        abssum += q_coef[scan_[n2 + subpos]];
      }

      if (last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if (last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan_[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1);
        if (signbit != (abssum & 0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff, min_pos = -1, cur_cost = 0x7fffffff;
          int16_t final_change = 0, cur_change = 0;
          for (n2 = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n2 >= 0; n2--) {
            uint32_t blkPos = scan_[n2 + subpos];
            if (q_coef[blkPos] != 0) {
              if (delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change = 1;
              }
              else if (n2 == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost = 0x7fffffff;
              }
              else {
                cur_cost = delta_u[blkPos];
                cur_change = -1;
              }
            }
            else if (n2 < first_nz_pos_in_cg && ((coef[blkPos] >= 0) ? 0 : 1) != signbit) {
              cur_cost = 0x7fffffff;
            }
            else {
              cur_cost = -delta_u[blkPos];
              cur_change = 1;
            }

            if (cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos = blkPos;
            }
          } // CG loop

          if (q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if (coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;
        } // Hide
      }
      if (last_cg == 1) last_cg = 0;
    }
  }
}

static void kvz_dequant_neon(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width,
  int32_t height, int8_t type, int8_t block_type)
{
  const encoder_control_t * const encoder = state->encoder_control;
  int32_t shift, add;
  const int32_t n_total = width * height;
  int32_t transform_shift = 15 - encoder->bitdepth - (kvz_g_convert_to_bit[width] + 2);

  int32_t qp_scaled = kvz_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6, encoder->cfg.chroma_format == KVZ_CSP_420);

  shift = 20 - QUANT_SHIFT - transform_shift;

  int32_t n = 0;
  if (encoder->scaling_list.enable) {
    uint32_t log2_tr_size = kvz_g_convert_to_bit[width] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    const coeff_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size - 2][scalinglist_type][qp_scaled % 6];
    shift += 4;

    if (shift > qp_scaled / 6) {
      add = 1 << (shift - qp_scaled / 6 - 1);
      const int32x4_t sh = vdupq_n_s32(-(shift - qp_scaled / 6));
      for (; n + 4 <= n_total; n += 4) {
        const int32x4_t q32 = vmovl_s16(vld1_s16(q_coef + n));
        const int32x4_t d32 = vmovl_s16(vld1_s16(dequant_coef + n));
        const int32x4_t p = vmulq_s32(q32, d32);
        // (q * dq + add) >> (shift - qp/6), with add a rounding offset
        const int32x4_t r = vshlq_s32(vaddq_s32(p, vdupq_n_s32(add)), sh);
        vst1_s16(coef + n, vqmovn_s32(r));
      }
      for (; n < n_total; ++n) {
        int32_t coeff_q = ((q_coef[n] * dequant_coef[n]) + add) >> (shift - qp_scaled / 6);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
      }
    } else {
      const int32x4_t ls = vdupq_n_s32(qp_scaled / 6 - shift);
      for (; n + 4 <= n_total; n += 4) {
        const int32x4_t q32 = vmovl_s16(vld1_s16(q_coef + n));
        const int32x4_t d32 = vmovl_s16(vld1_s16(dequant_coef + n));
        // Clip to 16 bits before the left shift, as the generic code does.
        const int32x4_t p = vmovl_s16(vqmovn_s32(vmulq_s32(q32, d32)));
        const int32x4_t r = vshlq_s32(p, ls);
        vst1_s16(coef + n, vqmovn_s32(r));
      }
      for (; n < n_total; ++n) {
        int32_t coeff_q = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q << (qp_scaled / 6 - shift));
      }
    }
  } else {
    int32_t scale = kvz_g_inv_quant_scales[qp_scaled % 6] << (qp_scaled / 6);
    add = 1 << (shift - 1);
    const int32x4_t vscale = vdupq_n_s32(scale);
    const int32x4_t sh = vdupq_n_s32(-shift);
    const int32x4_t vadd = vdupq_n_s32(add);
    for (; n + 4 <= n_total; n += 4) {
      const int32x4_t q32 = vmovl_s16(vld1_s16(q_coef + n));
      const int32x4_t p = vaddq_s32(vmulq_s32(q32, vscale), vadd);
      vst1_s16(coef + n, vqmovn_s32(vshlq_s32(p, sh)));
    }
    for (; n < n_total; ++n) {
      int32_t coeff_q = (q_coef[n] * scale + add) >> shift;
      coef[n] = (coeff_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}

static uint32_t coeff_abs_sum_neon(const coeff_t *coeffs, size_t length)
{
  uint32x4_t acc = vdupq_n_u32(0);
  size_t i = 0;
  for (; i + 8 <= length; i += 8) {
    const int16x8_t v = vld1q_s16(coeffs + i);
    acc = vpadalq_u16(acc, vreinterpretq_u16_s16(vabsq_s16(v)));
  }
  uint32_t sum = vaddvq_u32(acc);
  for (; i < length; ++i) {
    sum += abs(coeffs[i]);
  }
  return sum;
}

static double fast_coeff_cost_neon(const coeff_t *coeff, int32_t width, uint64_t weights)
{
  // Weights are four 16-bit values packed in a 64-bit word. They are looked up
  // with a byte table, once for the low byte and once for the high byte.
  uint16_t w[4];
  w[0] = (weights >>  0) & 0xffff;
  w[1] = (weights >> 16) & 0xffff;
  w[2] = (weights >> 32) & 0xffff;
  w[3] = (weights >> 48) & 0xffff;

  const uint8_t lo_bytes[8] = { (uint8_t)(w[0] & 0xff), (uint8_t)(w[1] & 0xff), (uint8_t)(w[2] & 0xff), (uint8_t)(w[3] & 0xff), 0, 0, 0, 0 };
  const uint8_t hi_bytes[8] = { (uint8_t)(w[0] >> 8), (uint8_t)(w[1] >> 8), (uint8_t)(w[2] >> 8), (uint8_t)(w[3] >> 8), 0, 0, 0, 0 };
  const uint8x8_t tbl_lo = vld1_u8(lo_bytes);
  const uint8x8_t tbl_hi = vld1_u8(hi_bytes);
  const uint8x8_t three = vdup_n_u8(3);

  const int32_t n_total = width * width;
  uint32x4_t acc = vdupq_n_u32(0);
  int32_t i = 0;
  for (; i + 8 <= n_total; i += 8) {
    const int16x8_t v = vld1q_s16(coeff + i);
    const uint8x8_t idx = vmin_u8(vqmovun_s16(vqabsq_s16(v)), three);  // vqabs saturates -32768 to 32767, as the generic abs does for the > 3 clamp
    const uint8x8_t lo = vtbl1_u8(tbl_lo, idx);
    const uint8x8_t hi = vtbl1_u8(tbl_hi, idx);
    const uint16x8_t wv = vorrq_u16(vshlq_n_u16(vmovl_u8(hi), 8), vmovl_u8(lo));
    acc = vpadalq_u16(acc, wv);
  }
  uint32_t sum = vaddvq_u32(acc);
  for (; i < n_total; ++i) {
    int16_t curr = coeff[i];
    uint32_t curr_abs = abs(curr);
    if (curr_abs > 3) {
      curr_abs = 3;
    }
    sum += w[curr_abs];
  }
  return (double)sum / 256.0;
}

#endif // COMPILE_ARM_NEON

int kvz_strategy_register_quant_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_ARM_NEON
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "quant", "neon", 40, &kvz_quant_neon);
    success &= kvz_strategyselector_register(opaque, "dequant", "neon", 40, &kvz_dequant_neon);
    success &= kvz_strategyselector_register(opaque, "coeff_abs_sum", "neon", 40, &coeff_abs_sum_neon);
    success &= kvz_strategyselector_register(opaque, "fast_coeff_cost", "neon", 40, &fast_coeff_cost_neon);
  }
#endif
  return success;
}
