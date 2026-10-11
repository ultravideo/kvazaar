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
 * Differential tests for picture and quant strategies that have no dedicated
 * test file: every non-generic implementation is compared with the generic one
 * on random and extreme inputs.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "strategies/strategies-picture.h"
#include "strategies/strategies-quant.h"
#include "strategies/strategies-dct.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define BUF_STRIDE 128
#define BUF_SIZE (BUF_STRIDE * 80)
#define BUF_ORIGIN (4 * BUF_STRIDE + 8)

static kvz_pixel buf_a[BUF_SIZE];
static kvz_pixel buf_b[BUF_SIZE];
static kvz_pixel pred_mem[4][64 * 64 + 1024];
static coeff_t coeffs[32 * 32];

// Large enough for a 1080p luma plane, to cover whole-frame variance (VAQ).
#define VAR_LEN (1920 * 1080)
static kvz_pixel var_buf[VAR_LEN];

static uint32_t rng_state = 0x2545F491u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

// mode 0: random, mode 1: random 0/255, mode 2: a=255 and b=0 (maximum differences).
static void fill(int mode)
{
  for (int i = 0; i < BUF_SIZE; ++i) {
    if (mode == 0) {
      buf_a[i] = (kvz_pixel)(rnd() & 0xFF);
      buf_b[i] = (kvz_pixel)(rnd() & 0xFF);
    } else if (mode == 1) {
      buf_a[i] = (rnd() & 1) ? 255 : 0;
      buf_b[i] = (rnd() & 1) ? 255 : 0;
    } else {
      buf_a[i] = 255;
      buf_b[i] = 0;
    }
  }
  for (int k = 0; k < 4; ++k) {
    for (int i = 0; i < 64 * 64 + 1024; ++i) {
      pred_mem[k][i] = mode == 2 ? 0 : (kvz_pixel)(rnd() & 0xFF);
    }
  }
}

static void *find_impl(const char *type, const char *name)
{
  for (int i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, name) == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

static const strategy_t *cur;   // strategy under test
static void *ref_fn;            // generic implementation of the same type

#define CHECK_EQ_U(ref, opt, ...) do { \
  if ((ref) != (opt)) { \
    printf("\n    %s:%s mismatch ref=%u opt=%u ", cur->type, cur->strategy_name, \
           (unsigned)(ref), (unsigned)(opt)); \
    printf(__VA_ARGS__); printf("\n"); \
    FAILm("mismatch"); \
  } } while (0)

TEST test_pixel_var(void)
{
  pixel_var_func *ref = ref_fn, *opt = cur->fptr;
  static const uint32_t lens[] = { 1, 15, 16, 17, 1024, 4096, 4096 + 7, 352 * 288, VAR_LEN };
  for (int mode = 0; mode < 3; ++mode) {
    for (uint32_t i = 0; i < VAR_LEN; ++i) {
      var_buf[i] = mode == 0 ? (kvz_pixel)(rnd() & 0xFF) :
                   mode == 1 ? ((rnd() & 1) ? 255 : 0) : 255;
    }
    for (unsigned li = 0; li < sizeof(lens) / sizeof(lens[0]); ++li) {
      const double r = ref(var_buf, lens[li]);
      const double o = opt(var_buf, lens[li]);
      // Implementations need not be bit-identical (AVX2 works in float), but must
      // agree closely: an overflow shows up as a large error.
      const double tol = 1e-4 * fabs(r) + 1e-3;
      if (fabs(r - o) > tol) {
        printf("\n    pixel_var:%s len=%u mode=%d ref=%.9g opt=%.9g\n",
               cur->strategy_name, lens[li], mode, r, o);
        FAILm("pixel_var mismatch");
      }
    }
  }
  PASS();
}

TEST test_ver_sad(void)
{
  ver_sad_func *ref = ref_fn, *opt = cur->fptr;
  static const int sizes[] = { 4, 8, 12, 16, 24, 32, 48, 64 };
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    for (int wi = 0; wi < 8; ++wi) for (int h = 1; h <= 64; h += 7) {
      const int w = sizes[wi];
      CHECK_EQ_U(ref(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, w, h, BUF_STRIDE),
                 opt(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, w, h, BUF_STRIDE),
                 "w=%d h=%d mode=%d", w, h, mode);
    }
  }
  PASS();
}

TEST test_hor_sad(void)
{
  hor_sad_func *ref = ref_fn, *opt = cur->fptr;
  static const int sizes[] = { 4, 8, 12, 16, 24, 32, 48, 64 };
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    for (int wi = 0; wi < 8; ++wi) for (int h = 1; h <= 64; h += 9) {
      const int w = sizes[wi];
      // image.c calls hor_sad only when left or right is non-zero; the no-border
      // case goes through reg_sad.
      for (int side = 1; side < 3; ++side) for (int border = 1; border < w; border += 3) {
        const uint32_t left = side == 1 ? (uint32_t)border : 0;
        const uint32_t right = side == 2 ? (uint32_t)border : 0;
        CHECK_EQ_U(ref(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, w, h, BUF_STRIDE, BUF_STRIDE, left, right),
                   opt(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, w, h, BUF_STRIDE, BUF_STRIDE, left, right),
                   "w=%d h=%d left=%u right=%u mode=%d", w, h, left, right, mode);
      }
    }
  }
  PASS();
}

TEST test_get_optimized_sad(void)
{
  get_optimized_sad_func *opt = cur->fptr;
  reg_sad_func *ref = find_impl("reg_sad", "generic");
  static const int widths[] = { 0, 4, 8, 12, 16, 24, 32, 48, 64 };
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    for (int wi = 0; wi < 9; ++wi) {
      optimized_sad_func_ptr_t f = opt(widths[wi]);
      if (!f) continue; // No specialized function: the caller falls back to reg_sad.
      for (int h = 1; h <= 64; h += 3) {
        CHECK_EQ_U(widths[wi] ? ref(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, widths[wi], h, BUF_STRIDE, BUF_STRIDE) : 0u,
                   f(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, h, BUF_STRIDE, BUF_STRIDE),
                   "w=%d h=%d mode=%d", widths[wi], h, mode);
      }
    }
  }
  PASS();
}

TEST test_pixels_calc_ssd(void)
{
  pixels_calc_ssd_func *ref = ref_fn, *opt = cur->fptr;
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    for (int w = 4; w <= 64; w *= 2) {
      CHECK_EQ_U(ref(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, BUF_STRIDE, BUF_STRIDE, w),
                 opt(buf_a + BUF_ORIGIN, buf_b + BUF_ORIGIN, BUF_STRIDE, BUF_STRIDE, w),
                 "w=%d mode=%d", w, mode);
    }
  }
  PASS();
}

TEST test_nxn_dual(void)
{
  cost_pixel_nxn_multi_func *ref = ref_fn, *opt = cur->fptr;
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    unsigned rc[2] = { 0, 0 }, oc[2] = { 1, 1 };
    ref((pred_buffer)pred_mem[0], buf_a, 2, rc);
    opt((pred_buffer)pred_mem[0], buf_a, 2, oc);
    CHECK_EQ_U(rc[0], oc[0], "cost 0 mode=%d", mode);
    CHECK_EQ_U(rc[1], oc[1], "cost 1 mode=%d", mode);
  }
  PASS();
}

TEST test_satd_any_size_quad(void)
{
  cost_pixel_any_size_multi_func *ref = ref_fn, *opt = cur->fptr;
  static const int sizes[] = { 4, 8, 12, 16, 24, 32, 64 };
  for (int mode = 0; mode < 3; ++mode) {
    fill(mode);
    const kvz_pixel *preds[4] = { pred_mem[0], pred_mem[1], pred_mem[2], pred_mem[3] };
    for (int wi = 0; wi < 7; ++wi) for (int hi = 0; hi < 7; ++hi) {
      const int w = sizes[wi], h = sizes[hi];
      int8_t valid[4] = { 1, 1, 1, 1 };
      unsigned rc[4] = { 0 }, oc[4] = { 0 };
      ref(w, h, preds, 64, buf_a + BUF_ORIGIN, BUF_STRIDE, 4, rc, valid);
      opt(w, h, preds, 64, buf_a + BUF_ORIGIN, BUF_STRIDE, 4, oc, valid);
      for (int k = 0; k < 4; ++k) {
        CHECK_EQ_U(rc[k], oc[k], "w=%d h=%d pred=%d mode=%d", w, h, k, mode);
      }
    }
  }
  PASS();
}

TEST test_fast_coeff_cost(void)
{
  fast_coeff_cost_func *ref = ref_fn, *opt = cur->fptr;
  for (int mode = 0; mode < 3; ++mode) {
    for (int w = 4; w <= 32; w *= 2) {
      for (int i = 0; i < w * w; ++i) {
        coeffs[i] = mode == 0 ? (coeff_t)((int)(rnd() % 9) - 4) :
                    mode == 1 ? (coeff_t)((rnd() & 1) ? 32767 : -32768) : 0;
      }
      const uint64_t weights = ((uint64_t)rnd() << 32) ^ rnd();
      const double r = ref(coeffs, w, weights), o = opt(coeffs, w, weights);
      if (r != o) {
        printf("\n    fast_coeff_cost:%s w=%d mode=%d ref=%.6f opt=%.6f\n",
               cur->strategy_name, w, mode, r, o);
        FAILm("fast_coeff_cost mismatch");
      }
    }
  }
  PASS();
}


static int16_t tr_in[32 * 32], tr_ref[32 * 32], tr_opt[32 * 32];

// Forward transforms take 8-bit residuals (|x| <= 255). Inverse transforms take
// dequantized coefficients, which dequant clips to the full int16 range, so
// they are tested there, including all-extreme blocks that saturate.
TEST test_transform(void)
{
  dct_func *ref = ref_fn, *opt = cur->fptr;
  const bool inverse = strstr(cur->type, "idct") || strstr(cur->type, "inverse");
  int n = 4;
  if (strstr(cur->type, "8x8")) n = 8;
  if (strstr(cur->type, "16x16")) n = 16;
  if (strstr(cur->type, "32x32")) n = 32;
  for (int mode = 0; mode < 5; ++mode) {
    for (int trial = 0; trial < 16; ++trial) {
      for (int i = 0; i < n * n; ++i) {
        int v;
        if (!inverse) {
          v = mode == 0 ? (int)(rnd() % 511) - 255 : mode == 1 ? ((rnd() & 1) ? 255 : -255) :
              mode == 2 ? 255 : mode == 3 ? -255 : ((i ^ (i >> 3)) & 1 ? 255 : -255);
        } else {
          v = mode == 0 ? (int)(rnd() & 0xFFFF) - 32768 : mode == 1 ? ((rnd() & 1) ? 32767 : -32768) :
              mode == 2 ? 32767 : mode == 3 ? -32768 : (int)(rnd() % 2001) - 1000;
        }
        tr_in[i] = (int16_t)v;
      }
      memset(tr_ref, 0x55, sizeof(tr_ref));
      memset(tr_opt, 0x55, sizeof(tr_opt));
      ref(8, tr_in, tr_ref);
      opt(8, tr_in, tr_opt);
      for (int i = 0; i < n * n; ++i) {
        if (tr_ref[i] != tr_opt[i]) {
          printf("\n    %s:%s mode=%d trial=%d idx=%d ref=%d opt=%d\n", cur->type,
                 cur->strategy_name, mode, trial, i, tr_ref[i], tr_opt[i]);
          FAILm("transform mismatch");
        }
      }
    }
  }
  PASS();
}

// The dual functions take pred_buffer, kvz_pixel (*)[32 * 32], and are only
// used by intra search for blocks up to 32x32. The 64x64 variants cannot be
// called through that type without indexing past its rows, so they are not
// tested here.
static const char *dual_types[] = {
  "sad_4x4_dual", "sad_8x8_dual", "sad_16x16_dual", "sad_32x32_dual",
  "satd_4x4_dual", "satd_8x8_dual", "satd_16x16_dual", "satd_32x32_dual",
};

SUITE(strategy_diff_tests)
{
  // hor_sad_generic calls the selected kvz_reg_sad, so point it at generic.
  kvz_reg_sad = find_impl("reg_sad", "generic");

  for (volatile int i = 0; i < strategies.count; ++i) {
    cur = &strategies.strategies[i];
    if (strcmp(cur->strategy_name, "generic") == 0) continue;
    ref_fn = find_impl(cur->type, "generic");
    if (!ref_fn && strcmp(cur->type, "get_optimized_sad") != 0) continue;

    if (strcmp(cur->type, "pixel_var") == 0) {
      RUN_TEST(test_pixel_var);
    } else if (strcmp(cur->type, "ver_sad") == 0) {
      RUN_TEST(test_ver_sad);
    } else if (strcmp(cur->type, "hor_sad") == 0) {
      RUN_TEST(test_hor_sad);
    } else if (strcmp(cur->type, "get_optimized_sad") == 0) {
      RUN_TEST(test_get_optimized_sad);
    } else if (strcmp(cur->type, "pixels_calc_ssd") == 0) {
      RUN_TEST(test_pixels_calc_ssd);
    } else if (strcmp(cur->type, "satd_any_size_quad") == 0) {
      RUN_TEST(test_satd_any_size_quad);
    } else if (strcmp(cur->type, "fast_coeff_cost") == 0) {
      RUN_TEST(test_fast_coeff_cost);
    } else if (!strncmp(cur->type, "dct_", 4) || !strncmp(cur->type, "idct_", 5) ||
               !strcmp(cur->type, "fast_forward_dst_4x4") || !strcmp(cur->type, "fast_inverse_dst_4x4")) {
      RUN_TEST(test_transform);
    } else {
      for (unsigned k = 0; k < sizeof(dual_types) / sizeof(dual_types[0]); ++k) {
        if (strcmp(cur->type, dual_types[k]) == 0) {
          RUN_TEST(test_nxn_dual);
          break;
        }
      }
    }
  }
}
