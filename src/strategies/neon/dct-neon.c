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
 * The 1-D kernels are the partial butterflies of dct-generic.c with every
 * scalar replaced by a vector of four independent rows. All arithmetic is exact
 * int32 integer arithmetic in the same order as the generic code, so results are
 * bit-identical. Each 2-D transform is two passes. Forward passes transpose the
 * input in 4x4 blocks and store the output rows directly. Inverse passes load the
 * input rows directly and transpose the output. Forward rounding wraps to 16 bits
 * like the (short) casts of the generic code, and inverse rounding saturates like
 * its CLIP.
 *
 * The implementation is bit-exact but is not registered by default. Under QEMU
 * user-mode emulation it runs slower than the generic code (see the notes in the
 * project documentation), so it is compiled only when KVZ_NEON_DCT is defined.
 * Enable it only after measuring a speedup on real AArch64 hardware.
 */

#include "strategies/neon/dct-neon.h"

#include "global.h" // IWYU pragma: keep

#if COMPILE_ARM_NEON

#include <string.h>

#include "kvazaar.h"
#include "strategies/generic/dct-generic.h"
#include "strategies/neon/neon-common.h"
#include "strategies/strategies-dct.h"
#include "strategyselector.h"
#include "tables.h"

#ifdef KVZ_NEON_DCT

typedef int32x4_t v4;

#define ADD(a, b)      vaddq_s32((a), (b))
#define SUB(a, b)      vsubq_s32((a), (b))
#define MULC(a, c)     vmulq_n_s32((a), (int32_t)(c))
#define MACC(acc, a, c) vmlaq_n_s32((acc), (a), (int32_t)(c))

// Transposes a 4x4 block of int16 lanes: out[j] lane i = in[i] lane j.
static inline void transpose4x4_s16(const int16x4_t in[4], int16x4_t out[4])
{
  int16x4x2_t a = vtrn_s16(in[0], in[1]);
  int16x4x2_t b = vtrn_s16(in[2], in[3]);
  int32x2x2_t c0 = vtrn_s32(vreinterpret_s32_s16(a.val[0]), vreinterpret_s32_s16(b.val[0]));
  int32x2x2_t c1 = vtrn_s32(vreinterpret_s32_s16(a.val[1]), vreinterpret_s32_s16(b.val[1]));
  out[0] = vreinterpret_s16_s32(c0.val[0]);
  out[1] = vreinterpret_s16_s32(c1.val[0]);
  out[2] = vreinterpret_s16_s32(c0.val[1]);
  out[3] = vreinterpret_s16_s32(c1.val[1]);
}

/* ---------------- Forward 1-D kernels: y[k] = sum_n M[k][n] * x[n] ---------------- */

static void fwd4_1d(const v4 *x, v4 *y)
{
  const v4 e0 = ADD(x[0], x[3]), o0 = SUB(x[0], x[3]);
  const v4 e1 = ADD(x[1], x[2]), o1 = SUB(x[1], x[2]);
  y[0] = MACC(MULC(e0, kvz_g_dct_4[0][0]), e1, kvz_g_dct_4[0][1]);
  y[2] = MACC(MULC(e0, kvz_g_dct_4[2][0]), e1, kvz_g_dct_4[2][1]);
  y[1] = MACC(MULC(o0, kvz_g_dct_4[1][0]), o1, kvz_g_dct_4[1][1]);
  y[3] = MACC(MULC(o0, kvz_g_dct_4[3][0]), o1, kvz_g_dct_4[3][1]);
}

static void fwd8_1d(const v4 *x, v4 *y)
{
  v4 e[4], o[4];
  for (int k = 0; k < 4; ++k) {
    e[k] = ADD(x[k], x[7 - k]);
    o[k] = SUB(x[k], x[7 - k]);
  }
  const v4 ee0 = ADD(e[0], e[3]), eo0 = SUB(e[0], e[3]);
  const v4 ee1 = ADD(e[1], e[2]), eo1 = SUB(e[1], e[2]);
  y[0] = MACC(MULC(ee0, kvz_g_dct_8[0][0]), ee1, kvz_g_dct_8[0][1]);
  y[4] = MACC(MULC(ee0, kvz_g_dct_8[4][0]), ee1, kvz_g_dct_8[4][1]);
  y[2] = MACC(MULC(eo0, kvz_g_dct_8[2][0]), eo1, kvz_g_dct_8[2][1]);
  y[6] = MACC(MULC(eo0, kvz_g_dct_8[6][0]), eo1, kvz_g_dct_8[6][1]);
  for (int k = 1; k < 8; k += 2) {
    v4 t = MULC(o[0], kvz_g_dct_8[k][0]);
    t = MACC(t, o[1], kvz_g_dct_8[k][1]);
    t = MACC(t, o[2], kvz_g_dct_8[k][2]);
    t = MACC(t, o[3], kvz_g_dct_8[k][3]);
    y[k] = t;
  }
}

static void fwd16_1d(const v4 *x, v4 *y)
{
  v4 e[8], o[8], ee[4], eo[4];
  for (int k = 0; k < 8; ++k) {
    e[k] = ADD(x[k], x[15 - k]);
    o[k] = SUB(x[k], x[15 - k]);
  }
  for (int k = 0; k < 4; ++k) {
    ee[k] = ADD(e[k], e[7 - k]);
    eo[k] = SUB(e[k], e[7 - k]);
  }
  const v4 eee0 = ADD(ee[0], ee[3]), eeo0 = SUB(ee[0], ee[3]);
  const v4 eee1 = ADD(ee[1], ee[2]), eeo1 = SUB(ee[1], ee[2]);
  y[0] = MACC(MULC(eee0, kvz_g_dct_16[0][0]), eee1, kvz_g_dct_16[0][1]);
  y[8] = MACC(MULC(eee0, kvz_g_dct_16[8][0]), eee1, kvz_g_dct_16[8][1]);
  y[4] = MACC(MULC(eeo0, kvz_g_dct_16[4][0]), eeo1, kvz_g_dct_16[4][1]);
  y[12] = MACC(MULC(eeo0, kvz_g_dct_16[12][0]), eeo1, kvz_g_dct_16[12][1]);
  for (int k = 2; k < 16; k += 4) {
    v4 t = MULC(eo[0], kvz_g_dct_16[k][0]);
    t = MACC(t, eo[1], kvz_g_dct_16[k][1]);
    t = MACC(t, eo[2], kvz_g_dct_16[k][2]);
    t = MACC(t, eo[3], kvz_g_dct_16[k][3]);
    y[k] = t;
  }
  for (int k = 1; k < 16; k += 2) {
    v4 t = MULC(o[0], kvz_g_dct_16[k][0]);
    for (int i = 1; i < 8; ++i) {
      t = MACC(t, o[i], kvz_g_dct_16[k][i]);
    }
    y[k] = t;
  }
}

static void fwd32_1d(const v4 *x, v4 *y)
{
  v4 e[16], o[16], ee[8], eo[8], eee[4], eeo[4];
  for (int k = 0; k < 16; ++k) {
    e[k] = ADD(x[k], x[31 - k]);
    o[k] = SUB(x[k], x[31 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    ee[k] = ADD(e[k], e[15 - k]);
    eo[k] = SUB(e[k], e[15 - k]);
  }
  for (int k = 0; k < 4; ++k) {
    eee[k] = ADD(ee[k], ee[7 - k]);
    eeo[k] = SUB(ee[k], ee[7 - k]);
  }
  const v4 eeee0 = ADD(eee[0], eee[3]), eeeo0 = SUB(eee[0], eee[3]);
  const v4 eeee1 = ADD(eee[1], eee[2]), eeeo1 = SUB(eee[1], eee[2]);

  y[0] = MACC(MULC(eeee0, kvz_g_dct_32[0][0]), eeee1, kvz_g_dct_32[0][1]);
  y[16] = MACC(MULC(eeee0, kvz_g_dct_32[16][0]), eeee1, kvz_g_dct_32[16][1]);
  y[8] = MACC(MULC(eeeo0, kvz_g_dct_32[8][0]), eeeo1, kvz_g_dct_32[8][1]);
  y[24] = MACC(MULC(eeeo0, kvz_g_dct_32[24][0]), eeeo1, kvz_g_dct_32[24][1]);
  for (int k = 4; k < 32; k += 8) {
    v4 t = MULC(eeo[0], kvz_g_dct_32[k][0]);
    for (int i = 1; i < 4; ++i) {
      t = MACC(t, eeo[i], kvz_g_dct_32[k][i]);
    }
    y[k] = t;
  }
  for (int k = 2; k < 32; k += 4) {
    v4 t = MULC(eo[0], kvz_g_dct_32[k][0]);
    for (int i = 1; i < 8; ++i) {
      t = MACC(t, eo[i], kvz_g_dct_32[k][i]);
    }
    y[k] = t;
  }
  for (int k = 1; k < 32; k += 2) {
    v4 t = MULC(o[0], kvz_g_dct_32[k][0]);
    for (int i = 1; i < 16; ++i) {
      t = MACC(t, o[i], kvz_g_dct_32[k][i]);
    }
    y[k] = t;
  }
}

/* ---------------- Inverse 1-D kernels: y[n] = sum_k M[k][n] * x[k] ---------------- */

static void inv4_1d(const v4 *x, v4 *y)
{
  const v4 o0 = MACC(MULC(x[1], kvz_g_dct_4[1][0]), x[3], kvz_g_dct_4[3][0]);
  const v4 o1 = MACC(MULC(x[1], kvz_g_dct_4[1][1]), x[3], kvz_g_dct_4[3][1]);
  const v4 e0 = MACC(MULC(x[0], kvz_g_dct_4[0][0]), x[2], kvz_g_dct_4[2][0]);
  const v4 e1 = MACC(MULC(x[0], kvz_g_dct_4[0][1]), x[2], kvz_g_dct_4[2][1]);
  y[0] = ADD(e0, o0);
  y[1] = ADD(e1, o1);
  y[2] = SUB(e1, o1);
  y[3] = SUB(e0, o0);
}

static void inv8_1d(const v4 *x, v4 *y)
{
  v4 o[4], e[4];
  for (int k = 0; k < 4; ++k) {
    v4 t = MULC(x[1], kvz_g_dct_8[1][k]);
    t = MACC(t, x[3], kvz_g_dct_8[3][k]);
    t = MACC(t, x[5], kvz_g_dct_8[5][k]);
    t = MACC(t, x[7], kvz_g_dct_8[7][k]);
    o[k] = t;
  }
  const v4 eo0 = MACC(MULC(x[2], kvz_g_dct_8[2][0]), x[6], kvz_g_dct_8[6][0]);
  const v4 eo1 = MACC(MULC(x[2], kvz_g_dct_8[2][1]), x[6], kvz_g_dct_8[6][1]);
  const v4 ee0 = MACC(MULC(x[0], kvz_g_dct_8[0][0]), x[4], kvz_g_dct_8[4][0]);
  const v4 ee1 = MACC(MULC(x[0], kvz_g_dct_8[0][1]), x[4], kvz_g_dct_8[4][1]);
  e[0] = ADD(ee0, eo0);
  e[3] = SUB(ee0, eo0);
  e[1] = ADD(ee1, eo1);
  e[2] = SUB(ee1, eo1);
  for (int k = 0; k < 4; ++k) {
    y[k] = ADD(e[k], o[k]);
    y[k + 4] = SUB(e[3 - k], o[3 - k]);
  }
}

static void inv16_1d(const v4 *x, v4 *y)
{
  v4 o[8], eo[4], e[8], ee[4];
  for (int k = 0; k < 8; ++k) {
    v4 t = MULC(x[1], kvz_g_dct_16[1][k]);
    for (int f = 3; f < 16; f += 2) {
      t = MACC(t, x[f], kvz_g_dct_16[f][k]);
    }
    o[k] = t;
  }
  for (int k = 0; k < 4; ++k) {
    v4 t = MULC(x[2], kvz_g_dct_16[2][k]);
    for (int f = 6; f < 16; f += 4) {
      t = MACC(t, x[f], kvz_g_dct_16[f][k]);
    }
    eo[k] = t;
  }
  const v4 eeo0 = MACC(MULC(x[4], kvz_g_dct_16[4][0]), x[12], kvz_g_dct_16[12][0]);
  const v4 eeo1 = MACC(MULC(x[4], kvz_g_dct_16[4][1]), x[12], kvz_g_dct_16[12][1]);
  const v4 eee0 = MACC(MULC(x[0], kvz_g_dct_16[0][0]), x[8], kvz_g_dct_16[8][0]);
  const v4 eee1 = MACC(MULC(x[0], kvz_g_dct_16[0][1]), x[8], kvz_g_dct_16[8][1]);
  const v4 eee[2] = { eee0, eee1 };
  const v4 eeo[2] = { eeo0, eeo1 };
  for (int k = 0; k < 2; ++k) {
    ee[k] = ADD(eee[k], eeo[k]);
    ee[k + 2] = SUB(eee[1 - k], eeo[1 - k]);
  }
  for (int k = 0; k < 4; ++k) {
    e[k] = ADD(ee[k], eo[k]);
    e[k + 4] = SUB(ee[3 - k], eo[3 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    y[k] = ADD(e[k], o[k]);
    y[k + 8] = SUB(e[7 - k], o[7 - k]);
  }
}

static void inv32_1d(const v4 *x, v4 *y)
{
  v4 o[16], eo[8], eeo[4], eeeo[2], eeee[2], eee[4], ee[8], e[16];
  for (int k = 0; k < 16; ++k) {
    v4 t = MULC(x[1], kvz_g_dct_32[1][k]);
    for (int f = 3; f < 32; f += 2) {
      t = MACC(t, x[f], kvz_g_dct_32[f][k]);
    }
    o[k] = t;
  }
  for (int k = 0; k < 8; ++k) {
    v4 t = MULC(x[2], kvz_g_dct_32[2][k]);
    for (int f = 6; f < 32; f += 4) {
      t = MACC(t, x[f], kvz_g_dct_32[f][k]);
    }
    eo[k] = t;
  }
  for (int k = 0; k < 4; ++k) {
    v4 t = MULC(x[4], kvz_g_dct_32[4][k]);
    for (int f = 12; f < 32; f += 8) {
      t = MACC(t, x[f], kvz_g_dct_32[f][k]);
    }
    eeo[k] = t;
  }
  eeeo[0] = MACC(MULC(x[8], kvz_g_dct_32[8][0]), x[24], kvz_g_dct_32[24][0]);
  eeeo[1] = MACC(MULC(x[8], kvz_g_dct_32[8][1]), x[24], kvz_g_dct_32[24][1]);
  eeee[0] = MACC(MULC(x[0], kvz_g_dct_32[0][0]), x[16], kvz_g_dct_32[16][0]);
  eeee[1] = MACC(MULC(x[0], kvz_g_dct_32[0][1]), x[16], kvz_g_dct_32[16][1]);

  eee[0] = ADD(eeee[0], eeeo[0]);
  eee[3] = SUB(eeee[0], eeeo[0]);
  eee[1] = ADD(eeee[1], eeeo[1]);
  eee[2] = SUB(eeee[1], eeeo[1]);
  for (int k = 0; k < 4; ++k) {
    ee[k] = ADD(eee[k], eeo[k]);
    ee[k + 4] = SUB(eee[3 - k], eeo[3 - k]);
  }
  for (int k = 0; k < 8; ++k) {
    e[k] = ADD(ee[k], eo[k]);
    e[k + 8] = SUB(ee[7 - k], eo[7 - k]);
  }
  for (int k = 0; k < 16; ++k) {
    y[k] = ADD(e[k], o[k]);
    y[k + 16] = SUB(e[15 - k], o[15 - k]);
  }
}

/* ---------------- 4-point DST (fast_forward_dst_4 / fast_inverse_dst_4) ---------------- */

static void fwd_dst4_1d(const v4 *x, v4 *y)
{
  const v4 c0 = ADD(x[0], x[3]);
  const v4 c1 = ADD(x[1], x[3]);
  const v4 c2 = SUB(x[0], x[1]);
  const v4 c3 = MULC(x[2], 74);
  y[0] = ADD(ADD(MULC(c0, 29), MULC(c1, 55)), c3);
  y[1] = MULC(SUB(ADD(x[0], x[1]), x[3]), 74);
  y[2] = SUB(ADD(MULC(c2, 29), MULC(c0, 55)), c3);
  y[3] = ADD(SUB(MULC(c2, 55), MULC(c1, 29)), c3);
}

static void inv_dst4_1d(const v4 *x, v4 *y)
{
  const v4 c0 = ADD(x[0], x[2]);
  const v4 c1 = ADD(x[2], x[3]);
  const v4 c2 = SUB(x[0], x[3]);
  const v4 c3 = MULC(x[1], 74);
  y[0] = ADD(ADD(MULC(c0, 29), MULC(c1, 55)), c3);
  y[1] = ADD(SUB(MULC(c2, 55), MULC(c1, 29)), c3);
  y[2] = MULC(ADD(SUB(x[0], x[2]), x[3]), 74);
  y[3] = SUB(ADD(MULC(c0, 55), MULC(c2, 29)), c3);
}

/* ---------------- 2-D passes ---------------- */

typedef void (fwd1d_func)(const v4 *x, v4 *y);
typedef void (inv1d_func)(const v4 *x, v4 *y);

// Forward pass over an n x n block of int16 with a wrapping 16-bit store.
static void fwd_pass(const int16_t *src, int16_t *dst, int n, fwd1d_func *fn, int32_t shift)
{
  const int32x4_t nsh = vdupq_n_s32(-shift);
  v4 x[32], y[32];
  for (int j0 = 0; j0 < n; j0 += 4) {
    for (int c0 = 0; c0 < n; c0 += 4) {
      int16x4_t rows[4], cols[4];
      for (int l = 0; l < 4; ++l) {
        rows[l] = vld1_s16(src + (j0 + l) * n + c0);
      }
      transpose4x4_s16(rows, cols);
      for (int i = 0; i < 4; ++i) {
        x[c0 + i] = vmovl_s16(cols[i]);
      }
    }
    fn(x, y);
    for (int k = 0; k < n; ++k) {
      vst1_s16(dst + k * n + j0, vmovn_s32(vrshlq_s32(y[k], nsh)));
    }
  }
}

// Inverse pass over an n x n block of int16 with a saturating 16-bit store.
static void inv_pass(const int16_t *src, int16_t *dst, int n, inv1d_func *fn, int32_t shift)
{
  const int32x4_t nsh = vdupq_n_s32(-shift);
  v4 x[32], y[32];
  for (int j0 = 0; j0 < n; j0 += 4) {
    for (int k = 0; k < n; ++k) {
      x[k] = vmovl_s16(vld1_s16(src + k * n + j0));
    }
    fn(x, y);
    for (int c0 = 0; c0 < n; c0 += 4) {
      int16x4_t cols[4], rows[4];
      for (int i = 0; i < 4; ++i) {
        cols[i] = vqmovn_s32(vrshlq_s32(y[c0 + i], nsh));
      }
      transpose4x4_s16(cols, rows);
      for (int l = 0; l < 4; ++l) {
        vst1_s16(dst + (j0 + l) * n + c0, rows[l]);
      }
    }
  }
}

#define DCT_NXN_NEON(N) \
static void dct_ ## N ## x ## N ## _neon(int8_t bitdepth, const int16_t *input, int16_t *output) \
{ \
  int16_t tmp[N * N]; \
  const int32_t sh1 = kvz_g_convert_to_bit[N] + 1 + (bitdepth - 8); \
  const int32_t sh2 = kvz_g_convert_to_bit[N] + 8; \
  fwd_pass(input, tmp, N, fwd ## N ## _1d, sh1); \
  fwd_pass(tmp, output, N, fwd ## N ## _1d, sh2); \
}

#define IDCT_NXN_NEON(N) \
static void idct_ ## N ## x ## N ## _neon(int8_t bitdepth, const int16_t *input, int16_t *output) \
{ \
  int16_t tmp[N * N]; \
  const int32_t sh1 = 7; \
  const int32_t sh2 = 12 - (bitdepth - 8); \
  inv_pass(input, tmp, N, inv ## N ## _1d, sh1); \
  inv_pass(tmp, output, N, inv ## N ## _1d, sh2); \
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
  const int32_t sh1 = kvz_g_convert_to_bit[4] + 1 + (bitdepth - 8);
  const int32_t sh2 = kvz_g_convert_to_bit[4] + 8;
  fwd_pass(input, tmp, 4, fwd_dst4_1d, sh1);
  fwd_pass(tmp, output, 4, fwd_dst4_1d, sh2);
}

static void fast_inverse_dst_4x4_neon(int8_t bitdepth, const int16_t *input, int16_t *output)
{
  int16_t tmp[4 * 4];
  const int32_t sh1 = 7;
  const int32_t sh2 = 12 - (bitdepth - 8);
  inv_pass(input, tmp, 4, inv_dst4_1d, sh1);
  inv_pass(tmp, output, 4, inv_dst4_1d, sh2);
}

#endif // KVZ_NEON_DCT

#endif //COMPILE_ARM_NEON

/* Defined on every platform so that non-ARM builds link. The call site only runs
 * when the NEON hardware flag is set, which is never the case off ARM. */
int kvz_strategy_register_dct_neon(void* opaque, uint8_t bitdepth)
{
  bool success = true;

#if COMPILE_ARM_NEON && defined(KVZ_NEON_DCT)
  if (bitdepth == 8) {
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
