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
 * Differential tests of the SIMD SATD and temporal filter strategies: every
 * non-generic implementation of these types is compared with the generic one
 * at KVZ_BIT_DEPTH, with random, extreme (0 / PIXEL_MAX, maximum differences,
 * Hadamard basis patterns that maximize the transform coefficients) and
 * edge-case inputs, for all block sizes, strides and alignments.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "strategies/strategies-picture.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define SS_STRIDE 96
#define SS_ROWS 84
#define SS_SIZE (SS_STRIDE * SS_ROWS)
// Origin of the blocks inside the buffers, with room for negative motion
// vectors and the extra column and row of the bilinear interpolation.
#define SS_ORIGIN (8 * SS_STRIDE + 8)
#define SS_NUM_MODES 8

static kvz_pixel ss_a[SS_SIZE];
static kvz_pixel ss_b[SS_SIZE];
static kvz_pixel ss_p[4][SS_SIZE];
// The dual functions take pred_buffer, kvz_pixel (*)[32 * 32]: the second
// prediction starts 1024 pixels after the first. A 64x64 block of the first
// prediction overlaps the second, so the buffer has room for that.
static kvz_pixel ss_dual[64 * 64 + 1024];
static kvz_pixel ss_dual_orig[64 * 64];

static uint32_t ss_rng = 0x1234567u;
static unsigned ss_u, ss_v;

static uint32_t ss_rand(void)
{
  ss_rng = ss_rng * 1664525u + 1013904223u;
  return ss_rng >> 8;
}

static kvz_pixel ss_rand_px(void)
{
  return (kvz_pixel)(ss_rand() % (PIXEL_MAX + 1));
}

static int ss_parity(unsigned v)
{
  int p = 0;
  for (; v; v >>= 1) p ^= v & 1;
  return p;
}

/* A Hadamard basis pattern of the 8x8 (or 4x4) blocks: 1 or 0. */
static int ss_basis(int x, int y)
{
  return ss_parity(((unsigned)x & ss_u) ^ ((unsigned)y & ss_v)) ? 0 : 1;
}

static kvz_pixel ss_clip(int v)
{
  return (kvz_pixel)(v < 0 ? 0 : (v > PIXEL_MAX ? PIXEL_MAX : v));
}

/*
 * Fills a (the original) and b (the prediction) according to mode:
 * 0 random, 1 random 0 / PIXEL_MAX, 2 a = PIXEL_MAX and b = 0, 3 a = 0 and
 * b = PIXEL_MAX, 4 a Hadamard basis pattern with maximum differences,
 * 5 small differences, 6 a basis pattern against a constant, 7 random against 0.
 * x and y are relative to (x0, y0).
 */
static void ss_fill_pair(kvz_pixel *a, kvz_pixel *b, int stride, int rows, int x0, int y0, int mode)
{
  ss_u = ss_rand() & 7;
  ss_v = ss_rand() & 7;
  for (int y = 0; y < rows; ++y) {
    for (int x = 0; x < stride; ++x) {
      const int i = y * stride + x;
      const int bx = x - x0, by = y - y0;
      switch (mode) {
        case 0: a[i] = ss_rand_px(); b[i] = ss_rand_px(); break;
        case 1: a[i] = (ss_rand() & 1) ? PIXEL_MAX : 0; b[i] = (ss_rand() & 1) ? PIXEL_MAX : 0; break;
        case 2: a[i] = PIXEL_MAX; b[i] = 0; break;
        case 3: a[i] = 0; b[i] = PIXEL_MAX; break;
        case 4: a[i] = ss_basis(bx, by) ? PIXEL_MAX : 0; b[i] = (kvz_pixel)(PIXEL_MAX - a[i]); break;
        case 5: a[i] = ss_rand_px(); b[i] = ss_clip(a[i] + (int)(ss_rand() % 7) - 3); break;
        case 6: a[i] = ss_basis(bx, by) ? PIXEL_MAX : 0; b[i] = (kvz_pixel)(PIXEL_MAX / 2); break;
        default: a[i] = ss_rand_px(); b[i] = 0; break;
      }
    }
  }
}

/* Fills ss_a as the original and the four predictions of ss_p. */
static void ss_fill_all(int mode)
{
  ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
  for (int k = 0; k < 4; ++k) {
    for (int i = 0; i < SS_SIZE; ++i) {
      switch (k) {
        case 0: ss_p[k][i] = ss_b[i]; break;
        case 1: ss_p[k][i] = (kvz_pixel)(PIXEL_MAX - ss_a[i]); break;
        case 2: ss_p[k][i] = mode == 2 ? 0 : (mode == 3 ? PIXEL_MAX : ss_rand_px()); break;
        default: ss_p[k][i] = (ss_rand() & 1) ? PIXEL_MAX : 0; break;
      }
    }
  }
}

static void *ss_find(const char *type, const char *name)
{
  for (unsigned i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, name) == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

static const strategy_t *ss_cur;  // strategy under test
static void *ss_ref;              // generic implementation of the same type

#define SS_CHECK(ref, opt, ...) do { \
  const unsigned long long ss_r_ = (unsigned long long)(ref), ss_o_ = (unsigned long long)(opt); \
  if (ss_r_ != ss_o_) { \
    printf("\n    %s:%s mismatch ref=%llu opt=%llu ", ss_cur->type, ss_cur->strategy_name, ss_r_, ss_o_); \
    printf(__VA_ARGS__); printf("\n"); \
    FAILm("mismatch"); \
  } } while (0)

static int ss_nxn_size(const char *type)
{
  int n = 0;
  sscanf(type, "satd_%d", &n);
  return n;
}

TEST ss_test_satd_nxn(void)
{
  cost_pixel_nxn_func *ref = ss_ref, *opt = ss_cur->fptr;
  const int n = ss_nxn_size(ss_cur->type);
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    for (int trial = 0; trial < 8; ++trial) {
      // Contiguous n x n blocks, at an odd offset for unaligned loads.
      ss_fill_pair(ss_a, ss_b, n, n + 1, 0, 0, mode);
      for (int off = 0; off < 2; ++off) {
        SS_CHECK(ref(ss_a + off, ss_b + off), opt(ss_a + off, ss_b + off), "n=%d mode=%d off=%d", n, mode, off);
        SS_CHECK(ref(ss_b + off, ss_a + off), opt(ss_b + off, ss_a + off), "n=%d mode=%d off=%d swapped", n, mode, off);
      }
    }
  }
  PASS();
}

TEST ss_test_satd_nxn_dual(void)
{
  cost_pixel_nxn_multi_func *ref = ss_ref, *opt = ss_cur->fptr;
  int n = 0;
  sscanf(ss_cur->type, "satd_%d", &n);
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    for (int trial = 0; trial < 8; ++trial) {
      ss_fill_pair(ss_dual_orig, ss_dual, n, n, 0, 0, mode);
      // The second prediction: the complement or random.
      for (int i = 0; i < n * n; ++i) {
        ss_dual[1024 + i] = (trial & 1) ? (kvz_pixel)(PIXEL_MAX - ss_dual_orig[i]) : ss_rand_px();
      }
      unsigned rc[2] = { 0, 0 }, oc[2] = { 1, 1 };
      ref((pred_buffer)ss_dual, ss_dual_orig, 2, rc);
      opt((pred_buffer)ss_dual, ss_dual_orig, 2, oc);
      SS_CHECK(rc[0], oc[0], "n=%d mode=%d cost 0", n, mode);
      SS_CHECK(rc[1], oc[1], "n=%d mode=%d cost 1", n, mode);
    }
  }
  PASS();
}

TEST ss_test_satd_any_size(void)
{
  cost_pixel_any_size_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
    for (int w = 4; w <= 64; w += 4) {
      for (int h = 4; h <= 64; h += 4) {
        for (int v = 0; v < 2; ++v) {
          // Different strides and an odd offset of the second block.
          const int s1 = SS_STRIDE, s2 = v ? SS_STRIDE - 3 : SS_STRIDE;
          const kvz_pixel *b1 = ss_a + SS_ORIGIN, *b2 = ss_b + SS_ORIGIN + v;
          SS_CHECK(ref(w, h, b1, s1, b2, s2), opt(w, h, b1, s1, b2, s2), "w=%d h=%d v=%d mode=%d", w, h, v, mode);
        }
      }
    }
  }
  PASS();
}

TEST ss_test_satd_any_size_quad(void)
{
  cost_pixel_any_size_multi_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_all(mode);
    for (int w = 4; w <= 64; w += 4) {
      for (int h = 4; h <= 64; h += 4) {
        for (int v = 0; v < 2; ++v) {
          const kvz_pixel *preds[4];
          for (int k = 0; k < 4; ++k) preds[k] = ss_p[k] + SS_ORIGIN + v * k;
          const int stride = v ? 67 : 64;
          int8_t valid[4] = { 1, 1, 1, 1 };
          unsigned rc[4] = { 0 }, oc[4] = { 1, 1, 1, 1 };
          ref(w, h, preds, stride, ss_a + SS_ORIGIN, SS_STRIDE, 4, rc, valid);
          opt(w, h, preds, stride, ss_a + SS_ORIGIN, SS_STRIDE, 4, oc, valid);
          for (int k = 0; k < 4; ++k) {
            SS_CHECK(rc[k], oc[k], "w=%d h=%d v=%d pred=%d mode=%d", w, h, v, k, mode);
          }
        }
      }
    }
  }
  PASS();
}

TEST ss_test_bilinear_ssd(void)
{
  bilinear_ssd_func *ref = ss_ref, *opt = ss_cur->fptr;
  static const int sizes[] = { 8, 16 };
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
    for (int si = 0; si < 2; ++si) {
      for (int f = 0; f < 16; ++f) {
        const int fx = f & 3, fy = f >> 2;
        for (int off = 0; off < 2; ++off) {
          const kvz_pixel *o = ss_a + SS_ORIGIN, *r = ss_b + SS_ORIGIN + off;
          SS_CHECK(ref(o, SS_STRIDE, r, SS_STRIDE, sizes[si], fx, fy),
                   opt(o, SS_STRIDE, r, SS_STRIDE, sizes[si], fx, fy),
                   "size=%d fx=%d fy=%d off=%d mode=%d", sizes[si], fx, fy, off, mode);
        }
      }
    }
  }
  PASS();
}

TEST ss_test_bilinear_ssd_multi(void)
{
  bilinear_ssd_multi_func *ref = ss_ref, *opt = ss_cur->fptr;
  static const int sizes[] = { 8, 16, 4 };
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
    for (int si = 0; si < 3; ++si) {
      for (int trial = 0; trial < 8; ++trial) {
        int mv_x[25], mv_y[25];
        const int num = 1 + (int)(ss_rand() % 25);
        for (int i = 0; i < num; ++i) {
          // Quarter-pel vectors within +-7 pixels, every other one integer.
          mv_x[i] = (int)(ss_rand() % 57) - 28;
          mv_y[i] = (int)(ss_rand() % 57) - 28;
          if (i & 1) { mv_x[i] &= ~3; mv_y[i] &= ~3; }
        }
        unsigned ro[25], oo[25];
        const kvz_pixel *o = ss_a + SS_ORIGIN, *r = ss_b + SS_ORIGIN + 8;
        ref(o, SS_STRIDE, r, SS_STRIDE, sizes[si], mv_x, mv_y, num, ro);
        opt(o, SS_STRIDE, r, SS_STRIDE, sizes[si], mv_x, mv_y, num, oo);
        for (int i = 0; i < num; ++i) {
          SS_CHECK(ro[i], oo[i], "size=%d mv=(%d,%d) mode=%d", sizes[si], mv_x[i], mv_y[i], mode);
        }
      }
    }
  }
  PASS();
}

static kvz_pixel ss_dst_ref[SS_SIZE];
static kvz_pixel ss_dst_opt[SS_SIZE];

TEST ss_test_tf_bilinear_pred(void)
{
  tf_bilinear_pred_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
    for (int w = 1; w <= 17; ++w) {
      for (int h = 1; h <= 10; ++h) {
        // Fraction denominators of 4 or 8 in each direction (shift 4 to 6),
        // and 2 (shift 1 to 3) for the extremes of the documented range.
        const int sx = 1 + (int)(ss_rand() % 3), sy = 1 + (int)(ss_rand() % 3);
        const int dx = 1 << sx, dy = 1 << sy;
        const int fx = (h & 1) ? dx - 1 : (int)(ss_rand() % dx), fy = (w & 1) ? dy - 1 : (int)(ss_rand() % dy);
        const int w00 = (dx - fx) * (dy - fy), w10 = fx * (dy - fy), w01 = (dx - fx) * fy, w11 = fx * fy;
        const int off = (int)(ss_rand() % 3);
        memset(ss_dst_ref, 0x5A, sizeof(ss_dst_ref));
        memset(ss_dst_opt, 0x5A, sizeof(ss_dst_opt));
        ref(ss_dst_ref + 1, SS_STRIDE, ss_b + SS_ORIGIN + off, SS_STRIDE, w, h, w00, w10, w01, w11, sx + sy);
        opt(ss_dst_opt + 1, SS_STRIDE, ss_b + SS_ORIGIN + off, SS_STRIDE, w, h, w00, w10, w01, w11, sx + sy);
        if (memcmp(ss_dst_ref, ss_dst_opt, sizeof(ss_dst_ref))) {
          printf("\n    %s:%s w=%d h=%d shift=%d mode=%d\n", ss_cur->type, ss_cur->strategy_name, w, h, sx + sy, mode);
          FAILm("mismatch");
        }
      }
    }
  }
  // Weights of 1 << shift for the largest shift (6), a full-scale copy.
  ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, 2);
  for (int k = 0; k < 4; ++k) {
    const int w[4] = { k == 0 ? 64 : 0, k == 1 ? 64 : 0, k == 2 ? 64 : 0, k == 3 ? 64 : 0 };
    memset(ss_dst_ref, 0, sizeof(ss_dst_ref));
    memset(ss_dst_opt, 0, sizeof(ss_dst_opt));
    ref(ss_dst_ref, SS_STRIDE, ss_a + SS_ORIGIN, SS_STRIDE, 8, 8, w[0], w[1], w[2], w[3], 6);
    opt(ss_dst_opt, SS_STRIDE, ss_a + SS_ORIGIN, SS_STRIDE, 8, 8, w[0], w[1], w[2], w[3], 6);
    ASSERT_EQ(0, memcmp(ss_dst_ref, ss_dst_opt, sizeof(ss_dst_ref)));
  }
  PASS();
}

TEST ss_test_tf_block_stats(void)
{
  tf_block_stats_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    for (int trial = 0; trial < 4; ++trial) {
      ss_fill_pair(ss_a, ss_b, SS_STRIDE, SS_ROWS, 8, 8, mode);
      if (trial == 3) {
        // Alternating maximum differences: the largest gradients.
        for (int i = 0; i < SS_SIZE; ++i) {
          const int x = i % SS_STRIDE, y = i / SS_STRIDE;
          ss_a[i] = ((x + y) & 1) ? PIXEL_MAX : 0;
          ss_b[i] = (kvz_pixel)(PIXEL_MAX - ss_a[i]);
        }
      }
      for (int w = 1; w <= 16; ++w) {
        for (int h = 1; h <= 16; ++h) {
          if ((w > 8 || h > 8) && w != h) continue;
          int64_t ro[5], oo[5];
          ref(ss_a + SS_ORIGIN, SS_STRIDE, ss_b + SS_ORIGIN + 1, SS_STRIDE - 2, w, h, ro);
          opt(ss_a + SS_ORIGIN, SS_STRIDE, ss_b + SS_ORIGIN + 1, SS_STRIDE - 2, w, h, oo);
          for (int k = 0; k < 5; ++k) {
            SS_CHECK(ro[k], oo[k], "w=%d h=%d out[%d] mode=%d trial=%d", w, h, k, mode, trial);
          }
        }
      }
    }
  }
  PASS();
}

static float ss_luts[9][2 * PIXEL_MAX + 1];

TEST ss_test_tf_filter_block(void)
{
  tf_filter_block_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_NUM_MODES; ++mode) {
    ss_fill_all(mode);
    for (int trial = 0; trial < 40; ++trial) {
      const int num_refs = 1 + trial % 9;  // up to 9 for the fallback
      const kvz_pixel *preds[9];
      const float *luts[9];
      float weights[9];
      for (int r = 0; r < num_refs; ++r) {
        // Like the encoder: a Gaussian of the difference, and weights that
        // may be 0.
        const float sigma = (float)(1 + ss_rand() % (2 * PIXEL_MAX));
        const float scale = 0.1f + (float)(ss_rand() % 100) / 100.0f;
        for (int d = -PIXEL_MAX; d <= PIXEL_MAX; ++d) {
          ss_luts[r][d + PIXEL_MAX] = scale * expf(-(float)d * (float)d / (2.0f * sigma * sigma));
        }
        luts[r] = &ss_luts[r][PIXEL_MAX];
        weights[r] = (ss_rand() % 5) == 0 ? 0.0f : (float)(ss_rand() % 2000) / 997.0f;
        preds[r] = ss_p[r & 3] + SS_ORIGIN + (r >> 2) * 3;
      }
      for (int w = 1; w <= 17; ++w) {
        const int h = 1 + (int)(ss_rand() % 9);
        memset(ss_dst_ref, 0x5A, sizeof(ss_dst_ref));
        memset(ss_dst_opt, 0x5A, sizeof(ss_dst_opt));
        ref(ss_dst_ref + 1, SS_STRIDE, ss_a + SS_ORIGIN, SS_STRIDE, preds, SS_STRIDE,
            num_refs, weights, luts, w, h, PIXEL_MAX);
        opt(ss_dst_opt + 1, SS_STRIDE, ss_a + SS_ORIGIN, SS_STRIDE, preds, SS_STRIDE,
            num_refs, weights, luts, w, h, PIXEL_MAX);
        if (memcmp(ss_dst_ref, ss_dst_opt, sizeof(ss_dst_ref))) {
          printf("\n    %s:%s w=%d h=%d refs=%d mode=%d\n", ss_cur->type, ss_cur->strategy_name, w, h, num_refs, mode);
          FAILm("mismatch");
        }
      }
    }
  }
  PASS();
}

SUITE(simd_satd_tests)
{
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    ss_cur = &strategies.strategies[i];
    if (strcmp(ss_cur->strategy_name, "generic") == 0) continue;
    const char *type = ss_cur->type;
    ss_ref = ss_find(type, "generic");
    if (!ss_ref) continue;
    ss_rng = 0x1234567u;

    if (!strncmp(type, "satd_", 5) && strstr(type, "_dual")) {
      RUN_TEST(ss_test_satd_nxn_dual);
    } else if (!strcmp(type, "satd_any_size")) {
      RUN_TEST(ss_test_satd_any_size);
    } else if (!strcmp(type, "satd_any_size_quad")) {
      RUN_TEST(ss_test_satd_any_size_quad);
    } else if (!strncmp(type, "satd_", 5) && ss_nxn_size(type) >= 4) {
      RUN_TEST(ss_test_satd_nxn);
    } else if (!strcmp(type, "bilinear_ssd")) {
      RUN_TEST(ss_test_bilinear_ssd);
    } else if (!strcmp(type, "bilinear_ssd_multi")) {
      RUN_TEST(ss_test_bilinear_ssd_multi);
    } else if (!strcmp(type, "tf_bilinear_pred")) {
      RUN_TEST(ss_test_tf_bilinear_pred);
    } else if (!strcmp(type, "tf_block_stats")) {
      RUN_TEST(ss_test_tf_block_stats);
    } else if (!strcmp(type, "tf_filter_block")) {
      RUN_TEST(ss_test_tf_filter_block);
    }
  }
}
