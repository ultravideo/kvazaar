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
 * AArch64 NEON forward and inverse transforms.
 *
 * Each 2-D transform is two 1-D passes over groups of four lines, one line
 * per vector lane. The arithmetic is the partial butterfly of dct-generic.c in
 * exact integer form, so the results are bit-identical:
 *
 * - Forward passes widen in the first butterfly stage (vaddl/vsubl) and
 *   multiply-accumulate in 32 bits. Each coefficient vector load serves four
 *   terms through lane-indexed multiplies. Output rounding wraps to 16 bits
 *   like the generic (short) cast.
 * - Inverse passes multiply the 16-bit inputs with 16-bit coefficients into
 *   32-bit accumulators (vmlal_laneq_s16), eight terms per coefficient load.
 *   Output rounding saturates like the generic CLIP.
 * - Inputs are transposed in 4x4 blocks on load (forward) or on store
 *   (inverse); for 4x4 transforms vld4/vst4 do the transpose.
 * - The strategies are registered for 8-bit only, so every pass shift is a
 *   compile-time constant (template parameter) and the rounding shift and
 *   narrowing of each output fuse into one vrshrn_n (wrapping, forward) or
 *   vqrshrn_n (saturating, inverse) instruction.
 */

#include "strategies/strategies-cpp.h"

extern "C" {
#include "strategies/neon/dct-neon.h"

#include "global.h" // IWYU pragma: keep

#if COMPILE_ARM_NEON
#include "kvazaar.h"
#include "strategies/generic/dct-generic.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-dct.h"
#include "strategyselector.h"
#include "tables.h"
#endif
}

#if COMPILE_ARM_NEON

namespace {

typedef int32x4_t v4;

// Pass shifts for 8-bit input (kvz_g_convert_to_bit[N] is log2(N) - 2).
template <int N> constexpr int log2_minus2() { return N == 4 ? 0 : N == 8 ? 1 : N == 16 ? 2 : 3; }
template <int N> constexpr int kFwdShift1 = log2_minus2<N>() + 1;
template <int N> constexpr int kFwdShift2 = log2_minus2<N>() + 8;
constexpr int kInvShift1 = 7;
constexpr int kInvShift2 = 12;

// 32-bit copies of the forward basis, so a row segment loads as a coefficient vector.
int32_t t4[4][4], t8[8][8], t16[16][16], t32[32][32];

void init_tables()
{
  for (int k = 0; k < 4; ++k) for (int i = 0; i < 4; ++i) t4[k][i] = kvz_g_dct_4[k][i];
  for (int k = 0; k < 8; ++k) for (int i = 0; i < 8; ++i) t8[k][i] = kvz_g_dct_8[k][i];
  for (int k = 0; k < 16; ++k) for (int i = 0; i < 16; ++i) t16[k][i] = kvz_g_dct_16[k][i];
  for (int k = 0; k < 32; ++k) for (int i = 0; i < 32; ++i) t32[k][i] = kvz_g_dct_32[k][i];
}

/* ---------------- Forward helpers ---------------- */

inline v4 dot2(v4 a, v4 b, const int32_t *c)
{
  const int32x2_t cv = vld1_s32(c);
  return vmlaq_lane_s32(vmulq_lane_s32(a, cv, 0), b, cv, 1);
}

inline v4 dot4_acc(v4 acc, const v4 *v, const int32_t *c)
{
  const int32x4_t cv = vld1q_s32(c);
  acc = vmlaq_laneq_s32(acc, v[0], cv, 0);
  acc = vmlaq_laneq_s32(acc, v[1], cv, 1);
  acc = vmlaq_laneq_s32(acc, v[2], cv, 2);
  return vmlaq_laneq_s32(acc, v[3], cv, 3);
}

inline v4 dot4(const v4 *v, const int32_t *c)
{
  const int32x4_t cv = vld1q_s32(c);
  v4 acc = vmulq_laneq_s32(v[0], cv, 0);
  acc = vmlaq_laneq_s32(acc, v[1], cv, 1);
  acc = vmlaq_laneq_s32(acc, v[2], cv, 2);
  return vmlaq_laneq_s32(acc, v[3], cv, 3);
}

template <int Terms>
inline v4 dotn(const v4 *v, const int32_t *c)
{
  v4 acc = dot4(v, c);
  KVZ_UNROLL
  for (int i = 4; i < Terms; i += 4) {
    acc = dot4_acc(acc, v + i, c + i);
  }
  return acc;
}

/* Forward 1-D kernels. x[n] holds sample n of four lines; out(k, val) rounds,
 * shifts and wraps output k of those lines to 16 bits and stores it. */
template <int N, int SH>
struct FwdOut {
  int16_t *d;
  inline void operator()(int k, v4 val) const { vst1_s16(d + k * N, vrshrn_n_s32(val, SH)); }
};

template <int SH>
KVZ_FORCE_INLINE void fwd4_1d(const int16x4_t *x, int16_t *d)
{
  const FwdOut<4, SH> out{d};
  const v4 e0 = vaddl_s16(x[0], x[3]), o0 = vsubl_s16(x[0], x[3]);
  const v4 e1 = vaddl_s16(x[1], x[2]), o1 = vsubl_s16(x[1], x[2]);
  out(0, dot2(e0, e1, t4[0]));
  out(2, dot2(e0, e1, t4[2]));
  out(1, dot2(o0, o1, t4[1]));
  out(3, dot2(o0, o1, t4[3]));
}

template <int SH>
KVZ_FORCE_INLINE void fwd8_1d(const int16x4_t *x, int16_t *d)
{
  const FwdOut<8, SH> out{d};
  v4 e[4], o[4];
  for (int k = 0; k < 4; ++k) {
    e[k] = vaddl_s16(x[k], x[7 - k]);
    o[k] = vsubl_s16(x[k], x[7 - k]);
  }
  const v4 ee0 = vaddq_s32(e[0], e[3]), eo0 = vsubq_s32(e[0], e[3]);
  const v4 ee1 = vaddq_s32(e[1], e[2]), eo1 = vsubq_s32(e[1], e[2]);
  out(0, dot2(ee0, ee1, t8[0]));
  out(4, dot2(ee0, ee1, t8[4]));
  out(2, dot2(eo0, eo1, t8[2]));
  out(6, dot2(eo0, eo1, t8[6]));
  for (int k = 1; k < 8; k += 2) out(k, dot4(o, t8[k]));
}

template <int SH>
KVZ_FORCE_INLINE void fwd16_1d(const int16x4_t *x, int16_t *d)
{
  const FwdOut<16, SH> out{d};
  v4 e[8], o[8], ee[4], eo[4];
  for (int k = 0; k < 8; ++k) {
    e[k] = vaddl_s16(x[k], x[15 - k]);
    o[k] = vsubl_s16(x[k], x[15 - k]);
  }
  for (int k = 1; k < 16; k += 2) out(k, dotn<8>(o, t16[k]));
  for (int k = 0; k < 4; ++k) {
    ee[k] = vaddq_s32(e[k], e[7 - k]);
    eo[k] = vsubq_s32(e[k], e[7 - k]);
  }
  for (int k = 2; k < 16; k += 4) out(k, dot4(eo, t16[k]));
  const v4 eee0 = vaddq_s32(ee[0], ee[3]), eeo0 = vsubq_s32(ee[0], ee[3]);
  const v4 eee1 = vaddq_s32(ee[1], ee[2]), eeo1 = vsubq_s32(ee[1], ee[2]);
  out(0, dot2(eee0, eee1, t16[0]));
  out(8, dot2(eee0, eee1, t16[8]));
  out(4, dot2(eeo0, eeo1, t16[4]));
  out(12, dot2(eeo0, eeo1, t16[12]));
}

template <int SH>
KVZ_FORCE_INLINE void fwd32_1d(const int16x4_t *x, int16_t *d)
{
  const FwdOut<32, SH> out{d};
  v4 e[16], o[16], ee[8], eo[8], eee[4], eeo[4];
  for (int k = 0; k < 16; ++k) {
    o[k] = vsubl_s16(x[k], x[31 - k]);
  }
  for (int k = 1; k < 32; k += 2) out(k, dotn<16>(o, t32[k]));
  for (int k = 0; k < 16; ++k) {
    e[k] = vaddl_s16(x[k], x[31 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    ee[k] = vaddq_s32(e[k], e[15 - k]);
    eo[k] = vsubq_s32(e[k], e[15 - k]);
  }
  for (int k = 2; k < 32; k += 4) out(k, dotn<8>(eo, t32[k]));
  for (int k = 0; k < 4; ++k) {
    eee[k] = vaddq_s32(ee[k], ee[7 - k]);
    eeo[k] = vsubq_s32(ee[k], ee[7 - k]);
  }
  for (int k = 4; k < 32; k += 8) out(k, dot4(eeo, t32[k]));
  const v4 eeee0 = vaddq_s32(eee[0], eee[3]), eeeo0 = vsubq_s32(eee[0], eee[3]);
  const v4 eeee1 = vaddq_s32(eee[1], eee[2]), eeeo1 = vsubq_s32(eee[1], eee[2]);
  out(0, dot2(eeee0, eeee1, t32[0]));
  out(16, dot2(eeee0, eeee1, t32[16]));
  out(8, dot2(eeeo0, eeeo1, t32[8]));
  out(24, dot2(eeeo0, eeeo1, t32[24]));
}

// 4-point forward DST (fast_forward_dst_4 in dct-generic.c).
template <int SH>
KVZ_FORCE_INLINE void fwd_dst4_1d(const int16x4_t *x, int16_t *d)
{
  const FwdOut<4, SH> out{d};
  const v4 c0 = vaddl_s16(x[0], x[3]);
  const v4 c1 = vaddl_s16(x[1], x[3]);
  const v4 c2 = vsubl_s16(x[0], x[1]);
  const v4 c3 = vmull_n_s16(x[2], 74);
  out(0, vaddq_s32(vmlaq_n_s32(vmulq_n_s32(c0, 29), c1, 55), c3));
  out(1, vmulq_n_s32(vsubq_s32(vaddl_s16(x[0], x[1]), vmovl_s16(x[3])), 74));
  out(2, vsubq_s32(vmlaq_n_s32(vmulq_n_s32(c2, 29), c0, 55), c3));
  out(3, vaddq_s32(vmlsq_n_s32(vmulq_n_s32(c2, 55), c1, 29), c3));
}

// Transposes a 4x4 block of int16 lanes: out[j] lane i = in[i] lane j.
inline void transpose4x4_s16(const int16x4_t in[4], int16x4_t out[4])
{
  const int16x4x2_t a = vtrn_s16(in[0], in[1]);
  const int16x4x2_t b = vtrn_s16(in[2], in[3]);
  const int32x2x2_t c0 = vtrn_s32(vreinterpret_s32_s16(a.val[0]), vreinterpret_s32_s16(b.val[0]));
  const int32x2x2_t c1 = vtrn_s32(vreinterpret_s32_s16(a.val[1]), vreinterpret_s32_s16(b.val[1]));
  out[0] = vreinterpret_s16_s32(c0.val[0]);
  out[1] = vreinterpret_s16_s32(c1.val[0]);
  out[2] = vreinterpret_s16_s32(c0.val[1]);
  out[3] = vreinterpret_s16_s32(c1.val[1]);
}

typedef void fwd1d_func(const int16x4_t *x, int16_t *d);

// One forward pass: line j of src (row j, stride N) becomes column j of dst.
template <int N, fwd1d_func Kernel>
KVZ_FORCE_INLINE void fwd_pass(const int16_t *src, int16_t *dst)
{
  int16x4_t x[N];
  for (int j0 = 0; j0 < N; j0 += 4) {
    if constexpr (N == 4) {
      const int16x4x4_t v = vld4_s16(src);
      x[0] = v.val[0]; x[1] = v.val[1]; x[2] = v.val[2]; x[3] = v.val[3];
    } else {
      for (int n0 = 0; n0 < N; n0 += 4) {
        const int16x4_t r[4] = {
          vld1_s16(src + (j0 + 0) * N + n0), vld1_s16(src + (j0 + 1) * N + n0),
          vld1_s16(src + (j0 + 2) * N + n0), vld1_s16(src + (j0 + 3) * N + n0),
        };
        transpose4x4_s16(r, x + n0);
      }
    }
    Kernel(x, dst + j0);
  }
}

/* ---------------- Inverse helpers ---------------- */

/* out[k] = sum over f = F0, F0 + STEP, ... < N of x[f] * M[f][k], for k < K.
 * Each coefficient load serves up to eight lane-indexed multiply-accumulates. */
template <int N, int K, int F0, int STEP>
KVZ_FORCE_INLINE void inv_sum(v4 *out, const int16x4_t *x, const int16_t *M)
{
  KVZ_UNROLL
  for (int k = 0; k < K; ++k) out[k] = vdupq_n_s32(0);
  KVZ_UNROLL
  for (int f = F0; f < N; f += STEP) {
    const int16_t *row = M + f * N;
    if constexpr (K >= 8) {
      KVZ_UNROLL
      for (int k0 = 0; k0 < K; k0 += 8) {
        const int16x8_t c = vld1q_s16(row + k0);
        v4 *acc = out + k0;
        acc[0] = vmlal_laneq_s16(acc[0], x[f], c, 0); acc[1] = vmlal_laneq_s16(acc[1], x[f], c, 1);
        acc[2] = vmlal_laneq_s16(acc[2], x[f], c, 2); acc[3] = vmlal_laneq_s16(acc[3], x[f], c, 3);
        acc[4] = vmlal_laneq_s16(acc[4], x[f], c, 4); acc[5] = vmlal_laneq_s16(acc[5], x[f], c, 5);
        acc[6] = vmlal_laneq_s16(acc[6], x[f], c, 6); acc[7] = vmlal_laneq_s16(acc[7], x[f], c, 7);
      }
    } else {
      const int16x4_t c = vld1_s16(row);
      out[0] = vmlal_lane_s16(out[0], x[f], c, 0);
      out[1] = vmlal_lane_s16(out[1], x[f], c, 1);
      if constexpr (K == 4) {
        out[2] = vmlal_lane_s16(out[2], x[f], c, 2);
        out[3] = vmlal_lane_s16(out[3], x[f], c, 3);
      }
    }
  }
}

/* Inverse 1-D kernels. x[k] holds coefficient k of four columns; y[n] is
 * output sample n of those columns. */
KVZ_FORCE_INLINE void inv4_1d(const int16x4_t *x, v4 *y)
{
  const int16x4_t c0 = vld1_s16(kvz_g_dct_4[0]), c1 = vld1_s16(kvz_g_dct_4[1]);
  const int16x4_t c2 = vld1_s16(kvz_g_dct_4[2]), c3 = vld1_s16(kvz_g_dct_4[3]);
  const v4 o0 = vmlal_lane_s16(vmull_lane_s16(x[1], c1, 0), x[3], c3, 0);
  const v4 o1 = vmlal_lane_s16(vmull_lane_s16(x[1], c1, 1), x[3], c3, 1);
  const v4 e0 = vmlal_lane_s16(vmull_lane_s16(x[0], c0, 0), x[2], c2, 0);
  const v4 e1 = vmlal_lane_s16(vmull_lane_s16(x[0], c0, 1), x[2], c2, 1);
  y[0] = vaddq_s32(e0, o0);
  y[1] = vaddq_s32(e1, o1);
  y[2] = vsubq_s32(e1, o1);
  y[3] = vsubq_s32(e0, o0);
}

KVZ_FORCE_INLINE void inv8_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_8[0][0];
  v4 o[4], eo[2], ee[2], e[4];
  inv_sum<8, 4, 1, 2>(o, x, M);
  inv_sum<8, 2, 2, 4>(eo, x, M);
  inv_sum<8, 2, 0, 4>(ee, x, M);
  e[0] = vaddq_s32(ee[0], eo[0]);
  e[3] = vsubq_s32(ee[0], eo[0]);
  e[1] = vaddq_s32(ee[1], eo[1]);
  e[2] = vsubq_s32(ee[1], eo[1]);
  for (int k = 0; k < 4; ++k) {
    y[k] = vaddq_s32(e[k], o[k]);
    y[k + 4] = vsubq_s32(e[3 - k], o[3 - k]);
  }
}

KVZ_FORCE_INLINE void inv16_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_16[0][0];
  v4 o[8], eo[4], eeo[2], eee[2], ee[4], e[8];
  inv_sum<16, 8, 1, 2>(o, x, M);
  inv_sum<16, 4, 2, 4>(eo, x, M);
  inv_sum<16, 2, 4, 8>(eeo, x, M);
  inv_sum<16, 2, 0, 8>(eee, x, M);
  for (int k = 0; k < 2; ++k) {
    ee[k] = vaddq_s32(eee[k], eeo[k]);
    ee[k + 2] = vsubq_s32(eee[1 - k], eeo[1 - k]);
  }
  for (int k = 0; k < 4; ++k) {
    e[k] = vaddq_s32(ee[k], eo[k]);
    e[k + 4] = vsubq_s32(ee[3 - k], eo[3 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    y[k] = vaddq_s32(e[k], o[k]);
    y[k + 8] = vsubq_s32(e[7 - k], o[7 - k]);
  }
}

KVZ_FORCE_INLINE void inv32_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_32[0][0];
  v4 o[16], eo[8], eeo[4], eeeo[2], eeee[2], eee[4], ee[8], e[16];
  inv_sum<32, 16, 1, 2>(o, x, M);
  inv_sum<32, 8, 2, 4>(eo, x, M);
  inv_sum<32, 4, 4, 8>(eeo, x, M);
  inv_sum<32, 2, 8, 16>(eeeo, x, M);
  inv_sum<32, 2, 0, 16>(eeee, x, M);
  eee[0] = vaddq_s32(eeee[0], eeeo[0]);
  eee[3] = vsubq_s32(eeee[0], eeeo[0]);
  eee[1] = vaddq_s32(eeee[1], eeeo[1]);
  eee[2] = vsubq_s32(eeee[1], eeeo[1]);
  for (int k = 0; k < 4; ++k) {
    ee[k] = vaddq_s32(eee[k], eeo[k]);
    ee[k + 4] = vsubq_s32(eee[3 - k], eeo[3 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    e[k] = vaddq_s32(ee[k], eo[k]);
    e[k + 8] = vsubq_s32(ee[7 - k], eo[7 - k]);
  }
  for (int k = 0; k < 16; ++k) {
    y[k] = vaddq_s32(e[k], o[k]);
    y[k + 16] = vsubq_s32(e[15 - k], o[15 - k]);
  }
}

// 4-point inverse DST (fast_inverse_dst_4 in dct-generic.c).
KVZ_FORCE_INLINE void inv_dst4_1d(const int16x4_t *x, v4 *y)
{
  const v4 c0 = vaddl_s16(x[0], x[2]);
  const v4 c1 = vaddl_s16(x[2], x[3]);
  const v4 c2 = vsubl_s16(x[0], x[3]);
  const v4 c3 = vmull_n_s16(x[1], 74);
  y[0] = vaddq_s32(vmlaq_n_s32(vmulq_n_s32(c0, 29), c1, 55), c3);
  y[1] = vaddq_s32(vmlsq_n_s32(vmulq_n_s32(c2, 55), c1, 29), c3);
  y[2] = vmulq_n_s32(vaddq_s32(vsubl_s16(x[0], x[2]), vmovl_s16(x[3])), 74);
  y[3] = vsubq_s32(vmlaq_n_s32(vmulq_n_s32(c0, 55), c2, 29), c3);
}

typedef void inv1d_func(const int16x4_t *x, v4 *y);

// One inverse pass: column j of src becomes line j of dst (row j, stride N).
template <int N, int SH, inv1d_func Kernel>
KVZ_FORCE_INLINE void inv_pass(const int16_t *src, int16_t *dst)
{
  int16x4_t x[N], cols[N];
  v4 y[N];
  for (int j0 = 0; j0 < N; j0 += 4) {
    for (int k = 0; k < N; ++k) {
      x[k] = vld1_s16(src + k * N + j0);
    }
    Kernel(x, y);
    for (int n = 0; n < N; ++n) {
      cols[n] = vqrshrn_n_s32(y[n], SH);
    }
    if constexpr (N == 4) {
      const int16x4x4_t v = { { cols[0], cols[1], cols[2], cols[3] } };
      vst4_s16(dst, v);
    } else {
      for (int n0 = 0; n0 < N; n0 += 4) {
        int16x4_t rows[4];
        transpose4x4_s16(cols + n0, rows);
        for (int l = 0; l < 4; ++l) {
          vst1_s16(dst + (j0 + l) * N + n0, rows[l]);
        }
      }
    }
  }
}

/* ---------------- 2-D transforms ---------------- */

// Passes and kernels are force-inlined into each transform (as in the C
// version) so that their vector arrays stay in registers.
template <int N> struct Dct;
template <> struct Dct<4> {
  template <int SH> KVZ_FORCE_INLINE static void fwd(const int16x4_t *x, int16_t *d) { fwd4_1d<SH>(x, d); }
  KVZ_FORCE_INLINE static void inv(const int16x4_t *x, v4 *y) { inv4_1d(x, y); }
};
template <> struct Dct<8> {
  template <int SH> KVZ_FORCE_INLINE static void fwd(const int16x4_t *x, int16_t *d) { fwd8_1d<SH>(x, d); }
  KVZ_FORCE_INLINE static void inv(const int16x4_t *x, v4 *y) { inv8_1d(x, y); }
};
template <> struct Dct<16> {
  template <int SH> KVZ_FORCE_INLINE static void fwd(const int16x4_t *x, int16_t *d) { fwd16_1d<SH>(x, d); }
  KVZ_FORCE_INLINE static void inv(const int16x4_t *x, v4 *y) { inv16_1d(x, y); }
};
template <> struct Dct<32> {
  template <int SH> KVZ_FORCE_INLINE static void fwd(const int16x4_t *x, int16_t *d) { fwd32_1d<SH>(x, d); }
  KVZ_FORCE_INLINE static void inv(const int16x4_t *x, v4 *y) { inv32_1d(x, y); }
};

template <int N>
void dct_nxn_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  (void)bitdepth; // registered for 8-bit only
  int16_t tmp[N * N];
  fwd_pass<N, Dct<N>::template fwd<kFwdShift1<N>>>(input, tmp);
  fwd_pass<N, Dct<N>::template fwd<kFwdShift2<N>>>(tmp, output);
}

template <int N>
void idct_nxn_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  (void)bitdepth;
  int16_t tmp[N * N];
  inv_pass<N, kInvShift1, Dct<N>::inv>(input, tmp);
  inv_pass<N, kInvShift2, Dct<N>::inv>(tmp, output);
}

void fast_forward_dst_4x4_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  (void)bitdepth;
  int16_t tmp[4 * 4];
  fwd_pass<4, fwd_dst4_1d<kFwdShift1<4>>>(input, tmp);
  fwd_pass<4, fwd_dst4_1d<kFwdShift2<4>>>(tmp, output);
}

void fast_inverse_dst_4x4_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  (void)bitdepth;
  int16_t tmp[4 * 4];
  inv_pass<4, kInvShift1, inv_dst4_1d>(input, tmp);
  inv_pass<4, kInvShift2, inv_dst4_1d>(tmp, output);
}

} // namespace

#endif //COMPILE_ARM_NEON

/* Defined on every platform so that non-ARM builds link. The call site only runs
 * when the NEON hardware flag is set, which is never the case off ARM. */
int kvz_strategy_register_dct_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

#if COMPILE_ARM_NEON
  if (bitdepth == 8) {
    init_tables();
    success &= kvz_strategyselector_register(opaque, "fast_forward_dst_4x4", "neon", 40, KVZ_STRATEGY_FPTR(fast_forward_dst_4x4_neon));
    success &= kvz_strategyselector_register(opaque, "dct_4x4", "neon", 40, KVZ_STRATEGY_FPTR(dct_nxn_neon<4>));
    success &= kvz_strategyselector_register(opaque, "dct_8x8", "neon", 40, KVZ_STRATEGY_FPTR(dct_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "dct_16x16", "neon", 40, KVZ_STRATEGY_FPTR(dct_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "dct_32x32", "neon", 40, KVZ_STRATEGY_FPTR(dct_nxn_neon<32>));

    success &= kvz_strategyselector_register(opaque, "fast_inverse_dst_4x4", "neon", 40, KVZ_STRATEGY_FPTR(fast_inverse_dst_4x4_neon));

    success &= kvz_strategyselector_register(opaque, "idct_4x4", "neon", 40, KVZ_STRATEGY_FPTR(idct_nxn_neon<4>));
    success &= kvz_strategyselector_register(opaque, "idct_8x8", "neon", 40, KVZ_STRATEGY_FPTR(idct_nxn_neon<8>));
    success &= kvz_strategyselector_register(opaque, "idct_16x16", "neon", 40, KVZ_STRATEGY_FPTR(idct_nxn_neon<16>));
    success &= kvz_strategyselector_register(opaque, "idct_32x32", "neon", 40, KVZ_STRATEGY_FPTR(idct_nxn_neon<32>));
  }
#else
  (void)opaque;
  (void)bitdepth;
#endif

  return success;
}
