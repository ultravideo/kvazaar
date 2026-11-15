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
 * Differential tests for the interpolation strategies: every non-generic
 * implementation of the block filters (fractional motion estimation) and of
 * the sample functions (motion compensation) is compared with the generic one
 * at the bit depth of the build, with random, extreme and worst-case inputs,
 * all block sizes, fractional positions and offsets.
 *
 * The source blocks are allocated with exactly the samples the generic
 * functions read (plus a little slack for the 8-bit AVX2 code, which reads
 * one sample more), so a sanitizer build also catches out-of-range reads.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "encoder.h"
#include "strategies/strategies-ipol.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if KVZ_BIT_DEPTH == 8
#define SRC_SLACK 16
#else
#define SRC_SLACK 0
#endif

static const strategy_t *cur;  // strategy under test
static void *ref_fn;           // generic implementation of the same type

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

static encoder_control_t enc;

static uint32_t rng_state = 0x9E3779B9u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

enum { FILL_RANDOM, FILL_EXTREME, FILL_HPEL_WORST, FILL_MAX, FILL_ZERO, FILL_SMOOTH, NUM_FILLS };

// Fills rows x cols samples. FILL_HPEL_WORST puts PIXEL_MAX under the positive
// taps of the half-pel filter (-1, 4, -11, 40, 40, -11, 4, -1) and 0 under the
// negative ones in both directions, which makes the generic diagonal sum wrap
// in its int16 store.
static void fill_block(kvz_pixel *buf, int stride, int rows, int cols, int mode)
{
  static const int pat[8] = { 0, 1, 0, 1, 1, 0, 1, 0 };
  const int phase_x = rnd() & 7, phase_y = rnd() & 7;
  for (int y = 0; y < rows; ++y) {
    for (int x = 0; x < cols; ++x) {
      kvz_pixel v;
      switch (mode) {
        case FILL_RANDOM: v = (kvz_pixel)(rnd() & PIXEL_MAX); break;
        case FILL_EXTREME: v = (rnd() & 1) ? PIXEL_MAX : 0; break;
        case FILL_HPEL_WORST: {
          const int px = pat[(x + phase_x) & 7], py = pat[(y + phase_y) & 7];
          v = (px == py) ? PIXEL_MAX : 0;
          break;
        }
        case FILL_MAX: v = PIXEL_MAX; break;
        case FILL_ZERO: v = 0; break;
        default: {
          // A gradient with noise. CLIP evaluates its argument several times.
          const int g = (x * 37 + y * 11) % (PIXEL_MAX + 1) + (int)(rnd() & 15) - 8;
          v = (kvz_pixel)CLIP(0, PIXEL_MAX, g);
          break;
        }
      }
      buf[y * stride + x] = v;
    }
  }
}

// A source block with pad samples on each side of a width x height block.
typedef struct {
  kvz_pixel *mem;
  kvz_pixel *origin;
  int stride;
} src_block_t;

static void alloc_src(src_block_t *b, int width, int height, int pad_l, int pad_r, int pad_t, int pad_b,
                      int extra_stride, int mode)
{
  const int cols = pad_l + width + pad_r;
  const int rows = pad_t + height + pad_b;
  b->stride = cols + extra_stride;
  const size_t size = (size_t)(rows - 1) * b->stride + cols + SRC_SLACK;
  b->mem = malloc(size * sizeof(kvz_pixel));
  // The gap between rows (extra_stride) is never used: fill it too, so that
  // a read there would at least be deterministic.
  for (size_t i = 0; i < size; ++i) b->mem[i] = (kvz_pixel)(rnd() & PIXEL_MAX);
  fill_block(b->mem, b->stride, rows, cols, mode);
  b->origin = b->mem + pad_t * b->stride + pad_l;
}

#define PRINT_MISMATCH(what, idx, ref, opt) \
  printf("\n    %s:%s %s mismatch at %d: ref=%d opt=%d ", cur->type, cur->strategy_name, \
         (what), (int)(idx), (int)(ref), (int)(opt))

/* ---- Sample functions (motion compensation) ---- */

#define DST_SIZE (LCU_WIDTH * LCU_WIDTH + 64)

static kvz_pixel dst_ref[DST_SIZE], dst_opt[DST_SIZE];
static int16_t dst_ref_hi[DST_SIZE], dst_opt_hi[DST_SIZE];

static const int luma_sizes[] = { 4, 8, 12, 16, 24, 32, 48, 64 };
static const int chroma_sizes[] = { 2, 4, 6, 8, 12, 16, 24, 32, 48, 64 };
#define NUM_LUMA_SIZES ((int)(sizeof(luma_sizes) / sizeof(luma_sizes[0])))
#define NUM_CHROMA_SIZES ((int)(sizeof(chroma_sizes) / sizeof(chroma_sizes[0])))

// Chroma configurations: (shift_w, shift_h). 1,1 is 4:2:0; 0,0 is 4:4:4; 1,0 is 4:2:2.
static const int chroma_shifts[][2] = { { 1, 1 }, { 0, 0 }, { 1, 0 } };

// Runs one sample function call for both implementations and compares the
// whole destination buffers (also catching writes outside the block).
static int run_sample(int luma, int hi, int width, int height, const int16_t mv[2], int mode, int variant)
{
  const int pad = luma ? KVZ_LUMA_FILTER_OFFSET : KVZ_CHROMA_FILTER_OFFSET;
  const int ext = luma ? KVZ_EXT_PADDING_LUMA : KVZ_EXT_PADDING_CHROMA;
  // The encoder provides one extra row for luma and three for chroma
  // (pad_b_simd in inter.c), which the 8-bit AVX2 functions read.
  const int pad_b_simd = luma ? 1 : 3;
  src_block_t b;
  alloc_src(&b, width, height, pad, ext - pad, pad, ext - pad + pad_b_simd, (variant % 3) * 5, mode);

  const int dst_stride = (variant & 1) ? width : LCU_WIDTH;
  const int8_t hor_flag = (int8_t)(mv[0] & (luma ? 3 : 7));
  const int8_t ver_flag = (int8_t)(mv[1] & (luma ? 3 : 7));
  int ok = 1;

  if (!hi) {
    memset(dst_ref, 0x5A, sizeof(dst_ref));
    memset(dst_opt, 0x5A, sizeof(dst_opt));
    kvz_sample_quarterpel_luma_func *ref = ref_fn, *opt = cur->fptr;
    ref(&enc, b.origin, (int16_t)b.stride, width, height, dst_ref, (int16_t)dst_stride, hor_flag, ver_flag, mv);
    opt(&enc, b.origin, (int16_t)b.stride, width, height, dst_opt, (int16_t)dst_stride, hor_flag, ver_flag, mv);
    for (int i = 0; i < DST_SIZE && ok; ++i) {
      if (dst_ref[i] != dst_opt[i]) {
        PRINT_MISMATCH("sample", i, dst_ref[i], dst_opt[i]);
        ok = 0;
      }
    }
  } else {
    memset(dst_ref_hi, 0x5A, sizeof(dst_ref_hi));
    memset(dst_opt_hi, 0x5A, sizeof(dst_opt_hi));
    kvz_sample_quarterpel_luma_hi_func *ref = ref_fn, *opt = cur->fptr;
    ref(&enc, b.origin, (int16_t)b.stride, width, height, dst_ref_hi, (int16_t)dst_stride, hor_flag, ver_flag, mv);
    opt(&enc, b.origin, (int16_t)b.stride, width, height, dst_opt_hi, (int16_t)dst_stride, hor_flag, ver_flag, mv);
    for (int i = 0; i < DST_SIZE && ok; ++i) {
      if (dst_ref_hi[i] != dst_opt_hi[i]) {
        PRINT_MISMATCH("sample", i, dst_ref_hi[i], dst_opt_hi[i]);
        ok = 0;
      }
    }
  }
  if (!ok) {
    printf("(%dx%d mv %d,%d fill %d stride %d dst_stride %d shift %d,%d)\n", width, height, mv[0], mv[1],
           mode, b.stride, dst_stride, enc.cfg.chroma_shift_w, enc.cfg.chroma_shift_h);
  }
  free(b.mem);
  return ok;
}

static int test_luma_sample(int hi)
{
  int variant = 0;
  memset(&enc, 0, sizeof(enc));
  for (int mode = 0; mode < NUM_FILLS; ++mode) {
    for (int wi = 0; wi < NUM_LUMA_SIZES; ++wi) {
      for (int hi_ = 0; hi_ < NUM_LUMA_SIZES; ++hi_) {
        for (int f = 0; f < 16; ++f) {
          // Negative and positive vectors with the same fraction.
          const int16_t mv[2] = { (int16_t)((f & 3) - ((variant & 2) ? 8 : -4)),
                                  (int16_t)((f >> 2) - ((variant & 4) ? 12 : -8)) };
          if (!run_sample(1, hi, luma_sizes[wi], luma_sizes[hi_], mv, mode, variant++)) return 0;
        }
      }
    }
  }
  return 1;
}

static int test_chroma_sample(int hi)
{
  int variant = 0;
  memset(&enc, 0, sizeof(enc));
  for (int s = 0; s < 3; ++s) {
    enc.cfg.chroma_shift_w = (uint8_t)chroma_shifts[s][0];
    enc.cfg.chroma_shift_h = (uint8_t)chroma_shifts[s][1];
    for (int mode = 0; mode < NUM_FILLS; ++mode) {
      for (int wi = 0; wi < NUM_CHROMA_SIZES; ++wi) {
        for (int hi_ = 0; hi_ < NUM_CHROMA_SIZES; ++hi_) {
          // The generic functions assume chroma blocks of at most LCU_WIDTH_C.
          if (chroma_sizes[wi] > (LCU_WIDTH >> chroma_shifts[s][0])) continue;
          if (chroma_sizes[hi_] > (LCU_WIDTH >> chroma_shifts[s][1])) continue;
          for (int f = 0; f < 64; ++f) {
            // Half of the fractions for the random inputs, a rotating subset otherwise.
            const int every = (mode == FILL_RANDOM || mode == FILL_EXTREME) ? 2 : 8;
            if (((f + wi * 3 + hi_ * 5) % every) != 0) continue;
            const int16_t mv[2] = { (int16_t)((f & 7) - ((variant & 2) ? 16 : -8)),
                                    (int16_t)((f >> 3) - ((variant & 4) ? 24 : -16)) };
            if (!run_sample(0, hi, chroma_sizes[wi], chroma_sizes[hi_], mv, mode, variant++)) return 0;
          }
        }
      }
    }
  }
  return 1;
}

TEST test_sample_quarterpel_luma(void)
{
  ASSERT(test_luma_sample(0));
  PASS();
}

TEST test_sample_quarterpel_luma_hi(void)
{
  ASSERT(test_luma_sample(1));
  PASS();
}

TEST test_sample_octpel_chroma(void)
{
  ASSERT(test_chroma_sample(0));
  PASS();
}

TEST test_sample_octpel_chroma_hi(void)
{
  ASSERT(test_chroma_sample(1));
  PASS();
}

/* ---- Block filters (fractional motion estimation) ---- */

typedef struct {
  ALIGNED(64) kvz_pixel filtered[4][LCU_LUMA_SIZE];
  ALIGNED(64) int16_t intermediate[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
  int16_t first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1];
} fme_state_t;

static fme_state_t fme_ref, fme_opt;

static const char *block_types[4] = {
  "filter_hpel_blocks_hor_ver_luma",
  "filter_hpel_blocks_diag_luma",
  "filter_qpel_blocks_hor_ver_luma",
  "filter_qpel_blocks_diag_luma",
};

static ipol_blocks_func *find_ref_block(int step)
{
  return (ipol_blocks_func *)find_impl(block_types[step], "generic");
}

static const int fme_sizes[] = { 8, 16, 24, 32, 48, 64 };
#define NUM_FME_SIZES ((int)(sizeof(fme_sizes) / sizeof(fme_sizes[0])))

// Compares the outputs and the intermediates the generic step writes. The
// intermediates are inputs of the later steps, so they must match too.
static int compare_fme(int step, int width, int height, int fme_level)
{
  for (int k = 0; k < 4; ++k) {
    for (int i = 0; i < LCU_LUMA_SIZE; ++i) {
      if (fme_ref.filtered[k][i] != fme_opt.filtered[k][i]) {
        PRINT_MISMATCH("filtered", i, fme_ref.filtered[k][i], fme_opt.filtered[k][i]);
        printf("(block %d, x %d, y %d)\n", k, i % LCU_WIDTH, i / LCU_WIDTH);
        return 0;
      }
    }
  }
  // Intermediates written by the horizontal steps: rows and the first column.
  int ims[2] = { -1, -1 };
  int cols[2] = { -1, -1 };
  int first_y = 0;
  if (step == 0) {
    ims[0] = 0; ims[1] = 1; cols[0] = 0; cols[1] = 2;
    first_y = fme_level > 1 ? 0 : 1;
  } else if (step == 2) {
    ims[0] = 3; ims[1] = 4; cols[0] = 1; cols[1] = 3;
  }
  const int rows = height + KVZ_EXT_PADDING_LUMA + 1;
  for (int j = 0; j < 2; ++j) {
    if (ims[j] < 0) continue;
    // Row 0 of the half-pel intermediate is only written for fme_level > 1.
    const int y0 = (step == 0 && j == 1) ? first_y : 0;
    for (int y = y0; y < rows; ++y) {
      for (int x = 0; x < width; ++x) {
        const int16_t r = fme_ref.intermediate[ims[j]][y * LCU_WIDTH + x];
        const int16_t o = fme_opt.intermediate[ims[j]][y * LCU_WIDTH + x];
        if (r != o) {
          PRINT_MISMATCH("intermediate", y * LCU_WIDTH + x, r, o);
          printf("(buffer %d)\n", ims[j]);
          return 0;
        }
      }
      const int16_t r = fme_ref.first_cols[cols[j]][y];
      const int16_t o = fme_opt.first_cols[cols[j]][y];
      if (r != o) {
        PRINT_MISMATCH("first column", y, r, o);
        printf("(buffer %d)\n", cols[j]);
        return 0;
      }
    }
  }
  return 1;
}

/*
 * Runs the four steps of the fractional search like search_frac does, the
 * generic functions on fme_ref and the implementations under test on fme_opt.
 * Before each step the generic state is copied to fme_opt, so every step is
 * tested on the same input, and the step under test uses intermediates
 * computed by the generic functions (the strategies are selected separately).
 */
static int run_fme(int width, int height, int fme_level, int8_t off_x, int8_t off_y, int mode, int extra_stride)
{
  src_block_t b;
  // One extra row and column on the top-left for the search; see search_frac.
  alloc_src(&b, width + 1, height + 1, KVZ_LUMA_FILTER_OFFSET, KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET,
            KVZ_LUMA_FILTER_OFFSET, KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET, extra_stride, mode);

  int ok = 1;
  memset(&fme_ref, 0x33, sizeof(fme_ref));
  for (int step = 0; step < fme_level && ok; ++step) {
    ipol_blocks_func *ref = find_ref_block(step);
    // Hpel steps use offsets 0; the qpel steps the best half-pel position.
    const int8_t ox = step < 2 ? 0 : off_x;
    const int8_t oy = step < 2 ? 0 : off_y;
    memcpy(&fme_opt, &fme_ref, sizeof(fme_ref));
    ref(&enc, b.origin, (int16_t)b.stride, width, height, fme_ref.filtered, fme_ref.intermediate,
        (int8_t)fme_level, fme_ref.first_cols, ox, oy);
    if (strcmp(cur->type, block_types[step]) == 0) {
      ipol_blocks_func *opt = cur->fptr;
      opt(&enc, b.origin, (int16_t)b.stride, width, height, fme_opt.filtered, fme_opt.intermediate,
          (int8_t)fme_level, fme_opt.first_cols, ox, oy);
      ok = compare_fme(step, width, height, fme_level);
      if (!ok) {
        printf("    (%dx%d fme_level %d off %d,%d fill %d stride %d)\n", width, height, fme_level,
               ox, oy, mode, b.stride);
      }
    }
  }
  free(b.mem);
  return ok;
}

TEST test_block_filter(void)
{
  int step = -1;
  for (int k = 0; k < 4; ++k) {
    if (strcmp(cur->type, block_types[k]) == 0) step = k;
  }
  ASSERT(step >= 0);
  memset(&enc, 0, sizeof(enc));

  int variant = 0;
  for (int mode = 0; mode < NUM_FILLS; ++mode) {
    for (int wi = 0; wi < NUM_FME_SIZES; ++wi) {
      for (int hi = 0; hi < NUM_FME_SIZES; ++hi) {
        // fme_level only matters for the first step (whether the top row of
        // the half-pel intermediate is needed); the later steps need a level
        // of at least step + 1.
        for (int fme_level = step + 1; fme_level <= (step == 0 ? 2 : step + 1); ++fme_level) {
          for (int off = 0; off < 9; ++off) {
            // The half-pel steps are called with zero offsets.
            if (step < 2 && off != 4) continue;
            const int8_t off_x = (int8_t)(off % 3 - 1), off_y = (int8_t)(off / 3 - 1);
            if (!run_fme(fme_sizes[wi], fme_sizes[hi], fme_level, off_x, off_y, mode, (variant++ % 3) * 7)) {
              FAILm("mismatch");
            }
          }
        }
      }
    }
  }
  PASS();
}

SUITE(simd_ipol_tests)
{
  for (volatile int i = 0; i < strategies.count; ++i) {
    cur = &strategies.strategies[i];
    if (strcmp(cur->strategy_name, "generic") == 0) continue;
    ref_fn = find_impl(cur->type, "generic");
    if (!ref_fn) continue;

    if (strcmp(cur->type, "sample_quarterpel_luma") == 0) {
      RUN_TEST(test_sample_quarterpel_luma);
    } else if (strcmp(cur->type, "sample_quarterpel_luma_hi") == 0) {
      RUN_TEST(test_sample_quarterpel_luma_hi);
    } else if (strcmp(cur->type, "sample_octpel_chroma") == 0) {
      RUN_TEST(test_sample_octpel_chroma);
    } else if (strcmp(cur->type, "sample_octpel_chroma_hi") == 0) {
      RUN_TEST(test_sample_octpel_chroma_hi);
    } else {
      for (int k = 0; k < 4; ++k) {
        if (strcmp(cur->type, block_types[k]) == 0) {
          RUN_TEST(test_block_filter);
          break;
        }
      }
    }
  }
}
