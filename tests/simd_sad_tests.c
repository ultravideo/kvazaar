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
 * Differential tests of the SAD, SSD, variance and bi-prediction strategies:
 * every non-generic implementation of reg_sad, sad_NxN(_dual),
 * get_optimized_sad, ver_sad, hor_sad, pixels_calc_ssd, pixel_var and
 * bipred_average is compared with the generic one at KVZ_BIT_DEPTH, with
 * random and extreme pixels, every block width up to 64, odd heights,
 * misaligned blocks and a zero reference stride.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "cu.h"
#include "image.h"
#include "strategies/strategies-picture.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define SS_STRIDE 160
#define SS_ROWS 72
#define SS_SIZE (SS_STRIDE * SS_ROWS)
#define SS_ORIGIN (4 * SS_STRIDE + 16)

static ALIGNED(32) kvz_pixel ss_a[SS_SIZE];
static ALIGNED(32) kvz_pixel ss_b[SS_SIZE];
// sad_64x64_dual reads preds[1] = preds[0] + 1024 up to 1024 + 64 * 64.
static ALIGNED(32) kvz_pixel ss_preds[6 * 1024];

#define SS_VAR_LEN (1920 * 1080 + 3)
static ALIGNED(32) kvz_pixel ss_var[SS_VAR_LEN];

static uint32_t ss_rng = 0x9E3779B9u;
static uint32_t ss_rnd(void)
{
  ss_rng = ss_rng * 1664525u + 1013904223u;
  return ss_rng >> 8;
}

#define SS_MODES 6

/*
 * mode 0: random, 1: random 0 or PIXEL_MAX, 2: a = PIXEL_MAX and b = 0,
 * 3: a = 0 and b = PIXEL_MAX, 4: all PIXEL_MAX, 5: a and b close to each other
 */
static kvz_pixel ss_pixel(int mode, int which)
{
  switch (mode) {
    case 0: return (kvz_pixel)(ss_rnd() & PIXEL_MAX);
    case 1: return (ss_rnd() & 1) ? PIXEL_MAX : 0;
    case 2: return which ? 0 : PIXEL_MAX;
    case 3: return which ? PIXEL_MAX : 0;
    case 4: return PIXEL_MAX;
    default: return (kvz_pixel)(PIXEL_MAX / 2 + (int)(ss_rnd() % 9) - 4);
  }
}

static void ss_fill(int mode)
{
  for (int i = 0; i < SS_SIZE; ++i) {
    ss_a[i] = ss_pixel(mode, 0);
    ss_b[i] = ss_pixel(mode, 1);
  }
  for (int i = 0; i < 6 * 1024; ++i) {
    ss_preds[i] = ss_pixel(mode, 1);
  }
}

static void *ss_find(const char *type, const char *name)
{
  for (int i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, name) == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

static const strategy_t *ss_cur;
static void *ss_ref;

#define SS_CHECK(ref, opt, ...) do { \
  const unsigned r_ = (unsigned)(ref), o_ = (unsigned)(opt); \
  if (r_ != o_) { \
    printf("\n    %s:%s mismatch ref=%u opt=%u ", ss_cur->type, ss_cur->strategy_name, r_, o_); \
    printf(__VA_ARGS__); printf("\n"); \
    FAILm("mismatch"); \
  } } while (0)

static const int ss_heights[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 23, 24, 31, 32, 33, 47, 48, 63, 64 };
#define SS_NUM_HEIGHTS ((int)(sizeof(ss_heights) / sizeof(ss_heights[0])))

TEST test_reg_sad(void)
{
  reg_sad_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    for (int w = 0; w <= 64; ++w) {
      for (int hi = 0; hi < SS_NUM_HEIGHTS; ++hi) {
        const int h = ss_heights[hi];
        const int off = (w + h) % 5;
        const kvz_pixel *a = ss_a + SS_ORIGIN + off;
        const kvz_pixel *b = ss_b + SS_ORIGIN + 3 - (off & 3);
        SS_CHECK(ref(a, b, w, h, SS_STRIDE, SS_STRIDE), opt(a, b, w, h, SS_STRIDE, SS_STRIDE),
                 "w=%d h=%d mode=%d", w, h, mode);
        SS_CHECK(ref(a, b, w, h, SS_STRIDE, 64), opt(a, b, w, h, SS_STRIDE, 64),
                 "w=%d h=%d mode=%d stride2=64", w, h, mode);
        SS_CHECK(ref(a, b, w, h, 64, 0), opt(a, b, w, h, 64, 0),
                 "w=%d h=%d mode=%d stride2=0", w, h, mode);
      }
    }
  }
  PASS();
}

TEST test_get_optimized_sad(void)
{
  get_optimized_sad_func *opt = ss_cur->fptr;
  reg_sad_func *ref = ss_find("reg_sad", "generic");
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    for (int w = 0; w <= 64; ++w) {
      optimized_sad_func_ptr_t f = opt(w);
      if (!f) continue; // The caller falls back to reg_sad.
      for (int hi = 0; hi < SS_NUM_HEIGHTS; ++hi) {
        const int h = ss_heights[hi];
        const kvz_pixel *a = ss_a + SS_ORIGIN + (h & 3);
        const kvz_pixel *b = ss_b + SS_ORIGIN + 1;
        SS_CHECK(w ? ref(a, b, w, h, SS_STRIDE, SS_STRIDE) : 0u,
                 f(a, b, h, SS_STRIDE, SS_STRIDE), "w=%d h=%d mode=%d", w, h, mode);
        SS_CHECK(w ? ref(a, b, w, h, 64, 96) : 0u,
                 f(a, b, h, 64, 96), "w=%d h=%d mode=%d strides 64/96", w, h, mode);
      }
    }
  }
  PASS();
}

TEST test_sad_nxn(void)
{
  cost_pixel_nxn_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    // The blocks are aligned like the encoder buffers: the 8-bit SIMD
    // versions use aligned loads.
    for (int k = 0; k < 4; ++k) {
      const kvz_pixel *a = ss_a + 64 * k, *b = ss_b + 128 * k;
      SS_CHECK(ref(a, b), opt(a, b), "mode=%d k=%d", mode, k);
    }
  }
  PASS();
}

TEST test_sad_nxn_dual(void)
{
  cost_pixel_nxn_multi_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    unsigned rc[2] = { 0, 0 }, oc[2] = { 1, 1 };
    ref((pred_buffer)ss_preds, ss_a, 2, rc);
    opt((pred_buffer)ss_preds, ss_a, 2, oc);
    SS_CHECK(rc[0], oc[0], "cost 0 mode=%d", mode);
    SS_CHECK(rc[1], oc[1], "cost 1 mode=%d", mode);
  }
  PASS();
}

TEST test_ver_sad(void)
{
  ver_sad_func *ref = ss_ref, *opt = ss_cur->fptr;
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    for (int w = 1; w <= 64; ++w) {
      for (int hi = 1; hi < SS_NUM_HEIGHTS; ++hi) {
        const int h = ss_heights[hi];
        const kvz_pixel *a = ss_a + SS_ORIGIN + (w & 3);
        const kvz_pixel *b = ss_b + SS_ORIGIN + 2;
        SS_CHECK(ref(a, b, w, h, SS_STRIDE), opt(a, b, w, h, SS_STRIDE), "w=%d h=%d mode=%d", w, h, mode);
      }
    }
  }
  PASS();
}

TEST test_hor_sad(void)
{
  hor_sad_func *ref = ss_ref, *opt = ss_cur->fptr;
  static const int widths[] = { 4, 8, 12, 16, 24, 32, 48, 64 };
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    for (int wi = 0; wi < 8; ++wi) {
      const int w = widths[wi];
      for (int hi = 1; hi < SS_NUM_HEIGHTS; hi += 2) {
        const int h = ss_heights[hi];
        // image.c calls hor_sad only when left or right is non-zero.
        for (int border = 1; border < w; ++border) {
          for (int side = 0; side < 2; ++side) {
            const uint32_t left = side == 0 ? (uint32_t)border : 0;
            const uint32_t right = side == 1 ? (uint32_t)border : 0;
            const kvz_pixel *a = ss_a + SS_ORIGIN + (border & 3);
            const kvz_pixel *b = ss_b + SS_ORIGIN + 1;
            SS_CHECK(ref(a, b, w, h, SS_STRIDE, 96, left, right),
                     opt(a, b, w, h, SS_STRIDE, 96, left, right),
                     "w=%d h=%d left=%u right=%u mode=%d", w, h, left, right, mode);
          }
        }
      }
    }
  }
  PASS();
}

TEST test_pixels_calc_ssd(void)
{
  pixels_calc_ssd_func *ref = ss_ref, *opt = ss_cur->fptr;
  // The encoder uses the power-of-two widths, which are all that the 8-bit
  // SIMD versions support.
#if KVZ_BIT_DEPTH > 8
  static const int widths[] = { 4, 8, 16, 32, 64, 1, 2, 3, 6, 12, 24, 48 };
#else
  static const int widths[] = { 4, 8, 16, 32, 64 };
#endif
  for (int mode = 0; mode < SS_MODES; ++mode) {
    ss_fill(mode);
    for (int wi = 0; wi < (int)(sizeof(widths) / sizeof(widths[0])); ++wi) {
      const int w = widths[wi];
      for (int off = 0; off < 3; ++off) {
        const kvz_pixel *a = ss_a + SS_ORIGIN + off;
        const kvz_pixel *b = ss_b + SS_ORIGIN + 2 * off;
        SS_CHECK(ref(a, b, SS_STRIDE, SS_STRIDE, w), opt(a, b, SS_STRIDE, SS_STRIDE, w),
                 "w=%d off=%d mode=%d", w, off, mode);
        SS_CHECK(ref(a, b, SS_STRIDE, w, w), opt(a, b, SS_STRIDE, w, w),
                 "w=%d off=%d mode=%d rec_stride=w", w, off, mode);
      }
    }
  }
  PASS();
}

TEST test_pixel_var(void)
{
  pixel_var_func *ref = ss_ref, *opt = ss_cur->fptr;
  static const uint32_t lens[] = { 1, 2, 7, 15, 16, 17, 31, 33, 64, 1024, 4096, 4096 + 7,
                                   16 * 1024 + 5, 65536, 65537, 352 * 288, SS_VAR_LEN };
  for (int mode = 0; mode < SS_MODES; ++mode) {
    for (uint32_t i = 0; i < SS_VAR_LEN; ++i) {
      ss_var[i] = ss_pixel(mode, i & 1);
    }
    for (unsigned li = 0; li < sizeof(lens) / sizeof(lens[0]); ++li) {
      for (int off = 0; off < 2; ++off) {
        const uint32_t len = lens[li] - (lens[li] > 1 ? off : 0);
        const double r = ref(ss_var + off, len);
        const double o = opt(ss_var + off, len);
        // Same tolerance as strategy_diff_tests: the implementations need not
        // round the same way, but an overflow shows up as a large error.
        const double tol = 1e-4 * fabs(r) + 1e-3;
        if (!(fabs(r - o) <= tol)) {
          printf("\n    pixel_var:%s len=%u mode=%d ref=%.9g opt=%.9g\n",
                 ss_cur->strategy_name, len, mode, r, o);
          FAILm("pixel_var mismatch");
        }
      }
    }
  }
  PASS();
}

// Bi-prediction: the source buffers are filled with pixels and with
// intermediate samples of various ranges, including the full int16 range.
static ALIGNED(32) kvz_pixel bp_px[2][3][LCU_LUMA_SIZE];
static ALIGNED(32) kvz_pixel_im bp_im[2][3][LCU_LUMA_SIZE];
static ALIGNED(32) kvz_pixel bp_rec_ref[3][LCU_LUMA_SIZE];
static ALIGNED(32) kvz_pixel bp_rec_opt[3][LCU_LUMA_SIZE];
static lcu_t bp_lcu_ref;
static lcu_t bp_lcu_opt;

static kvz_pixel_im bp_im_sample(int mode)
{
  switch (mode) {
    case 0: return (kvz_pixel_im)((int)(ss_rnd() % (1 << 14)) - 1024);     // typical range
    case 1: return (kvz_pixel_im)(int16_t)(ss_rnd() & 0xFFFF);              // any int16
    case 2: return (ss_rnd() & 1) ? 32767 : -32768;
    case 3: return 32767;
    case 4: return -32768;
    default: return (kvz_pixel_im)((ss_rnd() & 1) ? (16383 + (int)(ss_rnd() % 64)) : -(int)(ss_rnd() % 64));
  }
}

static void bp_fill(int mode)
{
  for (int l = 0; l < 2; ++l) {
    for (int c = 0; c < 3; ++c) {
      for (int i = 0; i < LCU_LUMA_SIZE; ++i) {
        bp_px[l][c][i] = ss_pixel(mode == 2 ? 2 + l : mode == 3 ? 4 : mode == 4 ? 1 : 0, l);
        bp_im[l][c][i] = bp_im_sample(mode);
      }
    }
  }
}

TEST test_bipred_average(void)
{
  inter_recon_bipred_func *ref = ss_ref, *opt = ss_cur->fptr;
  // All PU sizes that allow bi-prediction (not 8x4 and 4x8): square, SMP
  // and AMP.
  static const int sizes[][2] = {
    { 8, 8 }, { 16, 16 }, { 32, 32 }, { 64, 64 },
    { 16, 8 }, { 8, 16 }, { 32, 16 }, { 16, 32 }, { 64, 32 }, { 32, 64 },
    { 16, 4 }, { 16, 12 }, { 4, 16 }, { 12, 16 }, { 32, 8 }, { 32, 24 }, { 8, 32 }, { 24, 32 },
    { 64, 16 }, { 64, 48 }, { 16, 64 }, { 48, 64 },
  };
  const int num_sizes = (int)(sizeof(sizes) / sizeof(sizes[0]));
  static const uint8_t shifts[3][2] = { { 1, 1 }, { 1, 0 }, { 0, 0 } }; // 4:2:0, 4:2:2, 4:4:4

  yuv_t px[2];
  yuv_im_t im[2];
  for (int l = 0; l < 2; ++l) {
    px[l].size = LCU_LUMA_SIZE;
    px[l].y = bp_px[l][0]; px[l].u = bp_px[l][1]; px[l].v = bp_px[l][2];
    im[l].size = LCU_LUMA_SIZE;
    im[l].y = bp_im[l][0]; im[l].u = bp_im[l][1]; im[l].v = bp_im[l][2];
  }
  bp_lcu_ref.rec.y = bp_rec_ref[0]; bp_lcu_ref.rec.u = bp_rec_ref[1]; bp_lcu_ref.rec.v = bp_rec_ref[2];
  bp_lcu_opt.rec.y = bp_rec_opt[0]; bp_lcu_opt.rec.u = bp_rec_opt[1]; bp_lcu_opt.rec.v = bp_rec_opt[2];

  for (int mode = 0; mode < 6; ++mode) {
    bp_fill(mode);
    for (int si = 0; si < num_sizes; ++si) for (int pos = 0; pos < 2; ++pos) {
      const int w = sizes[si][0], h = sizes[si][1];
      // The top-left or the bottom-right corner of an LCU.
      const unsigned pu_x = 64 * (1 + si % 2) + (pos ? 64 - w : 0);
      const unsigned pu_y = 64 * (si % 3) + (pos ? 64 - h : 0);
      for (int csp = 0; csp < 3; ++csp) {
        for (unsigned flags = 0; flags < 16; ++flags) {
          const unsigned f0 = flags & 3, f1 = flags >> 2;
          const bool luma = (flags + si) % 3 != 0;
          const bool chroma = (flags + pos) % 3 != 1;
          memset(bp_rec_ref, 0x5A, sizeof(bp_rec_ref));
          memset(bp_rec_opt, 0x5A, sizeof(bp_rec_opt));
          ref(&bp_lcu_ref, &px[0], &px[1], &im[0], &im[1], pu_x, pu_y, w, h, f0, f1,
              luma, chroma, shifts[csp][0], shifts[csp][1]);
          opt(&bp_lcu_opt, &px[0], &px[1], &im[0], &im[1], pu_x, pu_y, w, h, f0, f1,
              luma, chroma, shifts[csp][0], shifts[csp][1]);
          for (int c = 0; c < 3; ++c) {
            for (int i = 0; i < LCU_LUMA_SIZE; ++i) {
              if (bp_rec_ref[c][i] != bp_rec_opt[c][i]) {
                printf("\n    bipred_average:%s w=%d h=%d x=%u y=%u csp=%d flags=%u/%u luma=%d chroma=%d "
                       "mode=%d plane=%d idx=%d ref=%u opt=%u\n",
                       ss_cur->strategy_name, w, h, pu_x, pu_y, csp, f0, f1, luma, chroma,
                       mode, c, i, (unsigned)bp_rec_ref[c][i], (unsigned)bp_rec_opt[c][i]);
                FAILm("bipred_average mismatch");
              }
            }
          }
        }
      }
    }
  }
  PASS();
}

static const char *ss_nxn_types[] = { "sad_4x4", "sad_8x8", "sad_16x16", "sad_32x32", "sad_64x64" };
static const char *ss_dual_types[] = {
  "sad_4x4_dual", "sad_8x8_dual", "sad_16x16_dual", "sad_32x32_dual", "sad_64x64_dual"
};

SUITE(simd_sad_tests)
{
  // hor_sad_generic calls the selected kvz_reg_sad, so point it at generic.
  reg_sad_func *saved_reg_sad = kvz_reg_sad;
  kvz_reg_sad = ss_find("reg_sad", "generic");

  for (volatile int i = 0; i < strategies.count; ++i) {
    ss_cur = &strategies.strategies[i];
    if (strcmp(ss_cur->strategy_name, "generic") == 0) continue;
    ss_ref = ss_find(ss_cur->type, "generic");
    if (!ss_ref) continue;

    const char *type = ss_cur->type;
    if (strcmp(type, "reg_sad") == 0) {
      RUN_TEST(test_reg_sad);
    } else if (strcmp(type, "get_optimized_sad") == 0) {
      RUN_TEST(test_get_optimized_sad);
    } else if (strcmp(type, "ver_sad") == 0) {
      RUN_TEST(test_ver_sad);
    } else if (strcmp(type, "hor_sad") == 0) {
      RUN_TEST(test_hor_sad);
    } else if (strcmp(type, "pixels_calc_ssd") == 0) {
      RUN_TEST(test_pixels_calc_ssd);
    } else if (strcmp(type, "pixel_var") == 0) {
      RUN_TEST(test_pixel_var);
    } else if (strcmp(type, "bipred_average") == 0) {
      RUN_TEST(test_bipred_average);
    } else {
      for (unsigned k = 0; k < 5; ++k) {
        if (strcmp(type, ss_nxn_types[k]) == 0) {
          RUN_TEST(test_sad_nxn);
        } else if (strcmp(type, ss_dual_types[k]) == 0) {
          RUN_TEST(test_sad_nxn_dual);
        }
      }
    }
  }

  kvz_reg_sad = saved_reg_sad;
}
