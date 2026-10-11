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
 */

#include "strategies/neon/dct-neon.h"

#include "global.h" // IWYU pragma: keep

#if COMPILE_ARM_NEON

#include "kvazaar.h"
#include "strategies/generic/dct-generic.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-dct.h"
#include "strategyselector.h"
#include "tables.h"

typedef int32x4_t v4;

// The passes are inlined into each transform so that N is a constant and the
// 1-D kernel call is direct; otherwise the pass overhead dominates 4x4.
#if defined(_MSC_VER)
#define DCT_INLINE __forceinline
#else
#define DCT_INLINE inline __attribute__((always_inline))
#endif

// 32-bit copies of the forward basis, so a row segment loads as a coefficient vector.
static int32_t t4[4][4], t8[8][8], t16[16][16], t32[32][32];

static void init_tables(void)
{
  for (int k = 0; k < 4; ++k) for (int i = 0; i < 4; ++i) t4[k][i] = kvz_g_dct_4[k][i];
  for (int k = 0; k < 8; ++k) for (int i = 0; i < 8; ++i) t8[k][i] = kvz_g_dct_8[k][i];
  for (int k = 0; k < 16; ++k) for (int i = 0; i < 16; ++i) t16[k][i] = kvz_g_dct_16[k][i];
  for (int k = 0; k < 32; ++k) for (int i = 0; i < 32; ++i) t32[k][i] = kvz_g_dct_32[k][i];
}

/* ---------------- Forward helpers ---------------- */

static inline v4 dot2(v4 a, v4 b, const int32_t *c)
{
  const int32x2_t cv = vld1_s32(c);
  return vmlaq_lane_s32(vmulq_lane_s32(a, cv, 0), b, cv, 1);
}

static inline v4 dot4_acc(v4 acc, const v4 *v, const int32_t *c)
{
  const int32x4_t cv = vld1q_s32(c);
  acc = vmlaq_laneq_s32(acc, v[0], cv, 0);
  acc = vmlaq_laneq_s32(acc, v[1], cv, 1);
  acc = vmlaq_laneq_s32(acc, v[2], cv, 2);
  return vmlaq_laneq_s32(acc, v[3], cv, 3);
}

static inline v4 dot4(const v4 *v, const int32_t *c)
{
  const int32x4_t cv = vld1q_s32(c);
  v4 acc = vmulq_laneq_s32(v[0], cv, 0);
  acc = vmlaq_laneq_s32(acc, v[1], cv, 1);
  acc = vmlaq_laneq_s32(acc, v[2], cv, 2);
  return vmlaq_laneq_s32(acc, v[3], cv, 3);
}

static inline v4 dotn(const v4 *v, const int32_t *c, int n)
{
  v4 acc = dot4(v, c);
  for (int i = 4; i < n; i += 4) {
    acc = dot4_acc(acc, v + i, c + i);
  }
  return acc;
}

// Rounds, wraps to 16 bits and stores output row k of this group of lines.
#define FWD_OUT(N, k, val) vst1_s16(d + (k) * (N), vmovn_s32(vrshlq_s32((val), nsh)))

/* Forward 1-D kernels. x[n] holds sample n of four lines; output k of those
 * lines goes to d[k * N]. */
static DCT_INLINE void fwd4_1d(const int16x4_t *x, int16_t *d, int32x4_t nsh)
{
  const v4 e0 = vaddl_s16(x[0], x[3]), o0 = vsubl_s16(x[0], x[3]);
  const v4 e1 = vaddl_s16(x[1], x[2]), o1 = vsubl_s16(x[1], x[2]);
  FWD_OUT(4, 0, dot2(e0, e1, t4[0]));
  FWD_OUT(4, 2, dot2(e0, e1, t4[2]));
  FWD_OUT(4, 1, dot2(o0, o1, t4[1]));
  FWD_OUT(4, 3, dot2(o0, o1, t4[3]));
}

static DCT_INLINE void fwd8_1d(const int16x4_t *x, int16_t *d, int32x4_t nsh)
{
  v4 e[4], o[4];
  for (int k = 0; k < 4; ++k) {
    e[k] = vaddl_s16(x[k], x[7 - k]);
    o[k] = vsubl_s16(x[k], x[7 - k]);
  }
  const v4 ee0 = vaddq_s32(e[0], e[3]), eo0 = vsubq_s32(e[0], e[3]);
  const v4 ee1 = vaddq_s32(e[1], e[2]), eo1 = vsubq_s32(e[1], e[2]);
  FWD_OUT(8, 0, dot2(ee0, ee1, t8[0]));
  FWD_OUT(8, 4, dot2(ee0, ee1, t8[4]));
  FWD_OUT(8, 2, dot2(eo0, eo1, t8[2]));
  FWD_OUT(8, 6, dot2(eo0, eo1, t8[6]));
  for (int k = 1; k < 8; k += 2) FWD_OUT(8, k, dot4(o, t8[k]));
}

static DCT_INLINE void fwd16_1d(const int16x4_t *x, int16_t *d, int32x4_t nsh)
{
  v4 e[8], o[8], ee[4], eo[4];
  for (int k = 0; k < 8; ++k) {
    e[k] = vaddl_s16(x[k], x[15 - k]);
    o[k] = vsubl_s16(x[k], x[15 - k]);
  }
  for (int k = 1; k < 16; k += 2) FWD_OUT(16, k, dotn(o, t16[k], 8));
  for (int k = 0; k < 4; ++k) {
    ee[k] = vaddq_s32(e[k], e[7 - k]);
    eo[k] = vsubq_s32(e[k], e[7 - k]);
  }
  for (int k = 2; k < 16; k += 4) FWD_OUT(16, k, dot4(eo, t16[k]));
  const v4 eee0 = vaddq_s32(ee[0], ee[3]), eeo0 = vsubq_s32(ee[0], ee[3]);
  const v4 eee1 = vaddq_s32(ee[1], ee[2]), eeo1 = vsubq_s32(ee[1], ee[2]);
  FWD_OUT(16, 0, dot2(eee0, eee1, t16[0]));
  FWD_OUT(16, 8, dot2(eee0, eee1, t16[8]));
  FWD_OUT(16, 4, dot2(eeo0, eeo1, t16[4]));
  FWD_OUT(16, 12, dot2(eeo0, eeo1, t16[12]));
}

static DCT_INLINE void fwd32_1d(const int16x4_t *x, int16_t *d, int32x4_t nsh)
{
  v4 e[16], o[16], ee[8], eo[8], eee[4], eeo[4];
  for (int k = 0; k < 16; ++k) {
    o[k] = vsubl_s16(x[k], x[31 - k]);
  }
  for (int k = 1; k < 32; k += 2) FWD_OUT(32, k, dotn(o, t32[k], 16));
  for (int k = 0; k < 16; ++k) {
    e[k] = vaddl_s16(x[k], x[31 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    ee[k] = vaddq_s32(e[k], e[15 - k]);
    eo[k] = vsubq_s32(e[k], e[15 - k]);
  }
  for (int k = 2; k < 32; k += 4) FWD_OUT(32, k, dotn(eo, t32[k], 8));
  for (int k = 0; k < 4; ++k) {
    eee[k] = vaddq_s32(ee[k], ee[7 - k]);
    eeo[k] = vsubq_s32(ee[k], ee[7 - k]);
  }
  for (int k = 4; k < 32; k += 8) FWD_OUT(32, k, dot4(eeo, t32[k]));
  const v4 eeee0 = vaddq_s32(eee[0], eee[3]), eeeo0 = vsubq_s32(eee[0], eee[3]);
  const v4 eeee1 = vaddq_s32(eee[1], eee[2]), eeeo1 = vsubq_s32(eee[1], eee[2]);
  FWD_OUT(32, 0, dot2(eeee0, eeee1, t32[0]));
  FWD_OUT(32, 16, dot2(eeee0, eeee1, t32[16]));
  FWD_OUT(32, 8, dot2(eeeo0, eeeo1, t32[8]));
  FWD_OUT(32, 24, dot2(eeeo0, eeeo1, t32[24]));
}

// 4-point forward DST (fast_forward_dst_4 in dct-generic.c).
static DCT_INLINE void fwd_dst4_1d(const int16x4_t *x, int16_t *d, int32x4_t nsh)
{
  const v4 c0 = vaddl_s16(x[0], x[3]);
  const v4 c1 = vaddl_s16(x[1], x[3]);
  const v4 c2 = vsubl_s16(x[0], x[1]);
  const v4 c3 = vmull_n_s16(x[2], 74);
  FWD_OUT(4, 0, vaddq_s32(vmlaq_n_s32(vmulq_n_s32(c0, 29), c1, 55), c3));
  FWD_OUT(4, 1, vmulq_n_s32(vsubq_s32(vaddl_s16(x[0], x[1]), vmovl_s16(x[3])), 74));
  FWD_OUT(4, 2, vsubq_s32(vmlaq_n_s32(vmulq_n_s32(c2, 29), c0, 55), c3));
  FWD_OUT(4, 3, vaddq_s32(vmlsq_n_s32(vmulq_n_s32(c2, 55), c1, 29), c3));
}

typedef void (fwd1d_func)(const int16x4_t *x, int16_t *d, int32x4_t nsh);

// Transposes a 4x4 block of int16 lanes: out[j] lane i = in[i] lane j.
static inline void transpose4x4_s16(const int16x4_t in[4], int16x4_t out[4])
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

// One forward pass: line j of src (row j, stride N) becomes column j of dst.
static DCT_INLINE void fwd_pass(const int16_t *src, int16_t *dst, int N, fwd1d_func *fn, int32_t shift)
{
  const int32x4_t nsh = vdupq_n_s32(-shift);
  int16x4_t x[32];
  for (int j0 = 0; j0 < N; j0 += 4) {
    if (N == 4) {
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
    fn(x, dst + j0, nsh);
  }
}

/* ---------------- Inverse helpers ---------------- */

// acc[i] += x * c[i] for i < 8 (c holds row segment M[f][k0 .. k0 + 7]).
#define MAC8(acc, xv, c) do { \
  acc[0] = vmlal_laneq_s16(acc[0], xv, c, 0); acc[1] = vmlal_laneq_s16(acc[1], xv, c, 1); \
  acc[2] = vmlal_laneq_s16(acc[2], xv, c, 2); acc[3] = vmlal_laneq_s16(acc[3], xv, c, 3); \
  acc[4] = vmlal_laneq_s16(acc[4], xv, c, 4); acc[5] = vmlal_laneq_s16(acc[5], xv, c, 5); \
  acc[6] = vmlal_laneq_s16(acc[6], xv, c, 6); acc[7] = vmlal_laneq_s16(acc[7], xv, c, 7); } while (0)
#define MAC4(acc, xv, c) do { \
  acc[0] = vmlal_lane_s16(acc[0], xv, c, 0); acc[1] = vmlal_lane_s16(acc[1], xv, c, 1); \
  acc[2] = vmlal_lane_s16(acc[2], xv, c, 2); acc[3] = vmlal_lane_s16(acc[3], xv, c, 3); } while (0)
#define MAC2(acc, xv, c) do { \
  acc[0] = vmlal_lane_s16(acc[0], xv, c, 0); acc[1] = vmlal_lane_s16(acc[1], xv, c, 1); } while (0)

static inline void zero(v4 *a, int n)
{
  for (int i = 0; i < n; ++i) a[i] = vdupq_n_s32(0);
}

// out[k] = sum over f = f0, f0 + step, ... < N of x[f] * M[f][k], for k < K.
static DCT_INLINE void inv_sum(v4 *out, int K, const int16x4_t *x, const int16_t *M, int N, int f0, int step)
{
  zero(out, K);
  for (int f = f0; f < N; f += step) {
    const int16_t *row = M + f * N;
    if (K >= 8) {
      for (int k0 = 0; k0 < K; k0 += 8) {
        const int16x8_t c = vld1q_s16(row + k0);
        MAC8((out + k0), x[f], c);
      }
    } else if (K == 4) {
      const int16x4_t c = vld1_s16(row);
      MAC4(out, x[f], c);
    } else {
      const int16x4_t c = vld1_s16(row);
      MAC2(out, x[f], c);
    }
  }
}

/* Inverse 1-D kernels. x[k] holds coefficient k of four columns; y[n] is
 * output sample n of those columns. */
static DCT_INLINE void inv4_1d(const int16x4_t *x, v4 *y)
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

static DCT_INLINE void inv8_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_8[0][0];
  v4 o[4], eo[2], ee[2], e[4];
  inv_sum(o, 4, x, M, 8, 1, 2);
  inv_sum(eo, 2, x, M, 8, 2, 4);
  inv_sum(ee, 2, x, M, 8, 0, 4);
  e[0] = vaddq_s32(ee[0], eo[0]);
  e[3] = vsubq_s32(ee[0], eo[0]);
  e[1] = vaddq_s32(ee[1], eo[1]);
  e[2] = vsubq_s32(ee[1], eo[1]);
  for (int k = 0; k < 4; ++k) {
    y[k] = vaddq_s32(e[k], o[k]);
    y[k + 4] = vsubq_s32(e[3 - k], o[3 - k]);
  }
}

static DCT_INLINE void inv16_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_16[0][0];
  v4 o[8], eo[4], eeo[2], eee[2], ee[4], e[8];
  inv_sum(o, 8, x, M, 16, 1, 2);
  inv_sum(eo, 4, x, M, 16, 2, 4);
  inv_sum(eeo, 2, x, M, 16, 4, 8);
  inv_sum(eee, 2, x, M, 16, 0, 8);
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

static DCT_INLINE void inv32_1d(const int16x4_t *x, v4 *y)
{
  const int16_t *M = &kvz_g_dct_32[0][0];
  v4 o[16], eo[8], eeo[4], eeeo[2], eeee[2], eee[4], ee[8], e[16];
  inv_sum(o, 16, x, M, 32, 1, 2);
  inv_sum(eo, 8, x, M, 32, 2, 4);
  inv_sum(eeo, 4, x, M, 32, 4, 8);
  inv_sum(eeeo, 2, x, M, 32, 8, 16);
  inv_sum(eeee, 2, x, M, 32, 0, 16);
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
static DCT_INLINE void inv_dst4_1d(const int16x4_t *x, v4 *y)
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

typedef void (inv1d_func)(const int16x4_t *x, v4 *y);

// One inverse pass: column j of src becomes line j of dst (row j, stride N).
static DCT_INLINE void inv_pass(const int16_t *src, int16_t *dst, int N, inv1d_func *fn, int32_t shift)
{
  const int32x4_t nsh = vdupq_n_s32(-shift);
  int16x4_t x[32], cols[32];
  v4 y[32];
  for (int j0 = 0; j0 < N; j0 += 4) {
    for (int k = 0; k < N; ++k) {
      x[k] = vld1_s16(src + k * N + j0);
    }
    fn(x, y);
    for (int n = 0; n < N; ++n) {
      cols[n] = vqmovn_s32(vrshlq_s32(y[n], nsh));
    }
    if (N == 4) {
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

#define DCT_NXN_NEON(N) \
static void dct_ ## N ## x ## N ## _neon(int8_t bitdepth, const int16_t *input, int16_t *output) \
{ \
  int16_t tmp[N * N]; \
  fwd_pass(input, tmp, N, fwd ## N ## _1d, kvz_g_convert_to_bit[N] + 1 + (bitdepth - 8)); \
  fwd_pass(tmp, output, N, fwd ## N ## _1d, kvz_g_convert_to_bit[N] + 8); \
}

#define IDCT_NXN_NEON(N) \
static void idct_ ## N ## x ## N ## _neon(int8_t bitdepth, const int16_t *input, int16_t *output) \
{ \
  int16_t tmp[N * N]; \
  inv_pass(input, tmp, N, inv ## N ## _1d, 7); \
  inv_pass(tmp, output, N, inv ## N ## _1d, 12 - (bitdepth - 8)); \
}

DCT_NXN_NEON(4)
DCT_NXN_NEON(8)
DCT_NXN_NEON(16)
DCT_NXN_NEON(32)

IDCT_NXN_NEON(4)
IDCT_NXN_NEON(8)
IDCT_NXN_NEON(16)
IDCT_NXN_NEON(32)

static void fast_forward_dst_4x4_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int16_t tmp[4 * 4];
  fwd_pass(input, tmp, 4, fwd_dst4_1d, kvz_g_convert_to_bit[4] + 1 + (bitdepth - 8));
  fwd_pass(tmp, output, 4, fwd_dst4_1d, kvz_g_convert_to_bit[4] + 8);
}

static void fast_inverse_dst_4x4_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int16_t tmp[4 * 4];
  inv_pass(input, tmp, 4, inv_dst4_1d, 7);
  inv_pass(tmp, output, 4, inv_dst4_1d, 12 - (bitdepth - 8));
}

#endif //COMPILE_ARM_NEON

/* Defined on every platform so that non-ARM builds link. The call site only runs
 * when the NEON hardware flag is set, which is never the case off ARM. */
int kvz_strategy_register_dct_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

#if COMPILE_ARM_NEON
  if (bitdepth == 8) {
    init_tables();
    success &= kvz_strategyselector_register(opaque, "fast_forward_dst_4x4", "neon", 40, &fast_forward_dst_4x4_neon);
    success &= kvz_strategyselector_register(opaque, "dct_4x4", "neon", 40, &dct_4x4_neon);
    success &= kvz_strategyselector_register(opaque, "dct_8x8", "neon", 40, &dct_8x8_neon);
    success &= kvz_strategyselector_register(opaque, "dct_16x16", "neon", 40, &dct_16x16_neon);
    success &= kvz_strategyselector_register(opaque, "dct_32x32", "neon", 40, &dct_32x32_neon);

    success &= kvz_strategyselector_register(opaque, "fast_inverse_dst_4x4", "neon", 40, &fast_inverse_dst_4x4_neon);

    success &= kvz_strategyselector_register(opaque, "idct_4x4", "neon", 40, &idct_4x4_neon);
    success &= kvz_strategyselector_register(opaque, "idct_8x8", "neon", 40, &idct_8x8_neon);
    success &= kvz_strategyselector_register(opaque, "idct_16x16", "neon", 40, &idct_16x16_neon);
    success &= kvz_strategyselector_register(opaque, "idct_32x32", "neon", 40, &idct_32x32_neon);
  }
#else
  (void)opaque;
  (void)bitdepth;
#endif

  return success;
}
