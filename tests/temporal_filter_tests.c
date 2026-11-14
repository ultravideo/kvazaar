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
 * \file
 * Tests that the SIMD versions of the temporal filter functions give the
 * same results as the generic ones.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/strategies-picture.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TF_TEST_STRIDE 48
#define TF_TEST_ROUNDS 200

static kvz_pixel tf_org[TF_TEST_STRIDE * TF_TEST_STRIDE];
static kvz_pixel tf_ref[TF_TEST_STRIDE * TF_TEST_STRIDE];
static kvz_pixel tf_preds[8][TF_TEST_STRIDE * TF_TEST_STRIDE];
static float tf_luts[8][2 * PIXEL_MAX + 1];
static unsigned tf_seed = 1;

static void *tf_tested;

static int tf_rand(void)
{
  tf_seed = tf_seed * 1103515245u + 12345u;
  return (int)((tf_seed >> 16) & 0x7fff);
}

static void tf_fill(kvz_pixel *buf, int smooth)
{
  for (int i = 0; i < TF_TEST_STRIDE * TF_TEST_STRIDE; i++) {
    // Mostly smooth content with noise, or full range noise. CLIP evaluates
    // its argument more than once, so draw the random number first.
    if (smooth) {
      const int val = 128 + (i % TF_TEST_STRIDE) * 2 - (i / TF_TEST_STRIDE) + tf_rand() % 21 - 10;
      buf[i] = (kvz_pixel)CLIP(0, PIXEL_MAX, val);
    } else {
      buf[i] = (kvz_pixel)(tf_rand() % (PIXEL_MAX + 1));
    }
  }
}


TEST tf_bilinear_ssd_multi_test(void)
{
  bilinear_ssd_multi_func *f = tf_tested;
  for (int round = 0; round < TF_TEST_ROUNDS; round++) {
    tf_fill(tf_org, round & 1);
    tf_fill(tf_ref, round & 2);
    const int size = (round % 3) == 0 ? 16 : 8;
    int mv_x[24], mv_y[24];
    const int num = 1 + tf_rand() % 24;
    for (int i = 0; i < num; i++) {
      mv_x[i] = tf_rand() % 33 - 16;
      mv_y[i] = tf_rand() % 33 - 16;
    }
    unsigned expected[24], got[24];
    const kvz_pixel *ref = &tf_ref[5 * TF_TEST_STRIDE + 5];
    kvz_bilinear_ssd_multi_generic(tf_org, TF_TEST_STRIDE, ref, TF_TEST_STRIDE, size, mv_x, mv_y, num, expected);
    f(tf_org, TF_TEST_STRIDE, ref, TF_TEST_STRIDE, size, mv_x, mv_y, num, got);
    for (int i = 0; i < num; i++) ASSERT_EQ(expected[i], got[i]);
  }
  PASS();
}


TEST tf_bilinear_pred_test(void)
{
  tf_bilinear_pred_func *f = tf_tested;
  static kvz_pixel expected[TF_TEST_STRIDE * TF_TEST_STRIDE], got[TF_TEST_STRIDE * TF_TEST_STRIDE];
  for (int round = 0; round < TF_TEST_ROUNDS; round++) {
    tf_fill(tf_ref, round & 1);
    static const int widths[] = { 8, 4, 2, 6, 16 };
    const int width = widths[round % 5];
    const int height = 1 + tf_rand() % 9;
    // Chroma uses 1/8 fractions in each direction (4:2:0), or 1/8 and 1/4
    // (4:2:2); luma uses 1/4.
    const int sx = 2 + tf_rand() % 2, sy = 2 + tf_rand() % 2;
    const int fx = tf_rand() % (1 << sx), fy = tf_rand() % (1 << sy);
    const int dx = 1 << sx, dy = 1 << sy;
    const int w00 = (dx - fx) * (dy - fy), w10 = fx * (dy - fy), w01 = (dx - fx) * fy, w11 = fx * fy;
    memset(expected, 0, sizeof(expected));
    memset(got, 0, sizeof(got));
    kvz_tf_bilinear_pred_generic(expected, TF_TEST_STRIDE, tf_ref, TF_TEST_STRIDE, width, height,
                                 w00, w10, w01, w11, sx + sy);
    f(got, TF_TEST_STRIDE, tf_ref, TF_TEST_STRIDE, width, height, w00, w10, w01, w11, sx + sy);
    ASSERT_EQ(0, memcmp(expected, got, sizeof(expected)));
  }
  PASS();
}


TEST tf_block_stats_test(void)
{
  tf_block_stats_func *f = tf_tested;
  for (int round = 0; round < TF_TEST_ROUNDS; round++) {
    tf_fill(tf_org, round & 1);
    tf_fill(tf_ref, round & 2);
    const int width = (round % 4) ? 8 : 1 + tf_rand() % 8;
    const int height = (round % 4) ? 8 : 1 + tf_rand() % 8;
    int64_t expected[5], got[5];
    kvz_tf_block_stats_generic(tf_org, TF_TEST_STRIDE, tf_ref, TF_TEST_STRIDE, width, height, expected);
    f(tf_org, TF_TEST_STRIDE, tf_ref, TF_TEST_STRIDE, width, height, got);
    for (int i = 0; i < 5; i++) ASSERT_EQ(expected[i], got[i]);
  }
  PASS();
}


TEST tf_filter_block_test(void)
{
  tf_filter_block_func *f = tf_tested;
  static kvz_pixel expected[TF_TEST_STRIDE * TF_TEST_STRIDE], got[TF_TEST_STRIDE * TF_TEST_STRIDE];
  for (int round = 0; round < TF_TEST_ROUNDS; round++) {
    tf_fill(tf_org, round & 1);
    const int num_refs = 1 + tf_rand() % 8;
    const kvz_pixel *preds[8];
    const float *luts[8];
    float weights[8];
    for (int r = 0; r < num_refs; r++) {
      tf_fill(tf_preds[r], (round >> 1) & 1);
      // Similar predictions in most rounds so that the weights matter.
      if (round & 4) {
        for (int i = 0; i < TF_TEST_STRIDE * TF_TEST_STRIDE; i++) {
          const int val = tf_org[i] + tf_rand() % 9 - 4;
          tf_preds[r][i] = (kvz_pixel)CLIP(0, PIXEL_MAX, val);
        }
      }
      const float sigma = (float)(20 + tf_rand() % 400);
      for (int d = -PIXEL_MAX; d <= PIXEL_MAX; d++) {
        tf_luts[r][d + PIXEL_MAX] = 0.4f * expf(-(float)(d * d) / (2.0f * sigma));
      }
      preds[r] = tf_preds[r];
      luts[r] = &tf_luts[r][PIXEL_MAX];
      weights[r] = (float)(tf_rand() % 1000) / 500.0f;
    }
    static const int widths[] = { 8, 4, 16, 3, 5 };
    const int width = widths[round % 5];
    const int height = 1 + tf_rand() % 8;
    memset(expected, 0, sizeof(expected));
    memset(got, 0, sizeof(got));
    kvz_tf_filter_block_generic(expected, TF_TEST_STRIDE, tf_org, TF_TEST_STRIDE, preds, TF_TEST_STRIDE,
                                num_refs, weights, luts, width, height, PIXEL_MAX);
    f(got, TF_TEST_STRIDE, tf_org, TF_TEST_STRIDE, preds, TF_TEST_STRIDE,
      num_refs, weights, luts, width, height, PIXEL_MAX);
    ASSERT_EQ(0, memcmp(expected, got, sizeof(expected)));
  }
  PASS();
}


SUITE(temporal_filter_tests)
{
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const char *type = strategies.strategies[i].type;
    tf_tested = strategies.strategies[i].fptr;
    tf_seed = 1;
    if (strcmp(type, "bilinear_ssd_multi") == 0) {
      RUN_TEST(tf_bilinear_ssd_multi_test);
    } else if (strcmp(type, "tf_bilinear_pred") == 0) {
      RUN_TEST(tf_bilinear_pred_test);
    } else if (strcmp(type, "tf_block_stats") == 0) {
      RUN_TEST(tf_block_stats_test);
    } else if (strcmp(type, "tf_filter_block") == 0) {
      RUN_TEST(tf_filter_block_test);
    }
  }
}
