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

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "encoder.h"
#include "strategies/generic/ipol-generic.h"
#include "strategies/strategies-ipol.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

// Source plane: the sampling origin sits 16 rows/columns into the buffer so the
// filters can read their 3 (luma) or 1 (chroma) samples of left/top padding.
#define SRC_STRIDE 256
#define SRC_ORIGIN (16 * SRC_STRIDE + 16)
#define DST_STRIDE LCU_WIDTH

static kvz_pixel src_buf[256 * 256];
static kvz_pixel dst_ref[LCU_WIDTH * LCU_WIDTH];
static kvz_pixel dst_opt[LCU_WIDTH * LCU_WIDTH];
static int16_t dst_ref_hi[LCU_WIDTH * LCU_WIDTH];
static int16_t dst_opt_hi[LCU_WIDTH * LCU_WIDTH];

static encoder_control_t enc;

static kvz_sample_quarterpel_luma_func *fn_luma = NULL;
static kvz_sample_octpel_chroma_func *fn_chroma = NULL;
static kvz_sample_quarterpel_luma_hi_func *fn_luma_hi = NULL;
static kvz_sample_octpel_chroma_hi_func *fn_chroma_hi = NULL;

static uint32_t rng_state = 0x12345678u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

// Mode 0: uniform random samples. Mode 1: extremes only (0 and 255), which
// stresses the intermediate range and the final clipping.
static void fill_src(int mode)
{
  for (int i = 0; i < 256 * 256; ++i) {
    if (mode == 0) {
      src_buf[i] = (kvz_pixel)(rnd() & 0xFF);
    } else {
      src_buf[i] = (rnd() & 1) ? 255 : 0;
    }
  }
}

static const int luma_sizes[] = { 4, 8, 12, 16, 32, 64 };
static const int luma_mvs[] = { -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 5, 6, 7, 13, 61 };
static const int chroma_sizes[] = { 2, 4, 6, 8, 16, 32 };
static const int chroma_mvs[] = { -9, -8, -5, -1, 0, 1, 2, 3, 4, 5, 7, 8, 11, 15, 23 };

// Chroma configurations: (shift_w, shift_h). 1,1 is 4:2:0; 0,0 is 4:4:4; 1,0 is 4:2:2.
static const int chroma_shifts[][2] = { { 1, 1 }, { 0, 0 }, { 1, 0 } };

static int compare_pixels(int w, int h, const char *what, int *mism_x, int *mism_y)
{
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (dst_ref[y * DST_STRIDE + x] != dst_opt[y * DST_STRIDE + x]) {
        *mism_x = x;
        *mism_y = y;
        printf("\n    %s mismatch at (%d,%d): ref=%d opt=%d (w=%d h=%d)\n", what, x, y,
               dst_ref[y * DST_STRIDE + x], dst_opt[y * DST_STRIDE + x], w, h);
        return 0;
      }
    }
  }
  return 1;
}

static int compare_hi(int w, int h, const char *what, int *mism_x, int *mism_y)
{
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (dst_ref_hi[y * DST_STRIDE + x] != dst_opt_hi[y * DST_STRIDE + x]) {
        *mism_x = x;
        *mism_y = y;
        printf("\n    %s mismatch at (%d,%d): ref=%d opt=%d (w=%d h=%d)\n", what, x, y,
               dst_ref_hi[y * DST_STRIDE + x], dst_opt_hi[y * DST_STRIDE + x], w, h);
        return 0;
      }
    }
  }
  return 1;
}

static int run_luma_pixel(void)
{
  int mx = 0, my = 0;
  for (int mode = 0; mode < 2; ++mode) {
    fill_src(mode);
    for (int wi = 0; wi < 6; ++wi) for (int hi = 0; hi < 6; ++hi) {
      for (int mxi = 0; mxi < 16; ++mxi) for (int myi = 0; myi < 16; ++myi) {
        const int w = luma_sizes[wi], h = luma_sizes[hi];
        const int16_t mv[2] = { (int16_t)luma_mvs[mxi], (int16_t)luma_mvs[myi] };
        memset(dst_ref, 0xCD, sizeof(dst_ref));
        memset(dst_opt, 0xCD, sizeof(dst_opt));
        kvz_sample_quarterpel_luma_generic(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_ref, DST_STRIDE, 0, 0, mv);
        fn_luma(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_opt, DST_STRIDE, 0, 0, mv);
        if (!compare_pixels(w, h, "luma", &mx, &my)) return 0;
      }
    }
  }
  return 1;
}

static int run_luma_hi(void)
{
  int mx = 0, my = 0;
  for (int mode = 0; mode < 2; ++mode) {
    fill_src(mode);
    for (int wi = 0; wi < 6; ++wi) for (int hi = 0; hi < 6; ++hi) {
      for (int mxi = 0; mxi < 16; ++mxi) for (int myi = 0; myi < 16; ++myi) {
        const int w = luma_sizes[wi], h = luma_sizes[hi];
        const int16_t mv[2] = { (int16_t)luma_mvs[mxi], (int16_t)luma_mvs[myi] };
        memset(dst_ref_hi, 0xCD, sizeof(dst_ref_hi));
        memset(dst_opt_hi, 0xCD, sizeof(dst_opt_hi));
        kvz_sample_quarterpel_luma_hi_generic(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_ref_hi, DST_STRIDE, 0, 0, mv);
        fn_luma_hi(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_opt_hi, DST_STRIDE, 0, 0, mv);
        if (!compare_hi(w, h, "luma_hi", &mx, &my)) return 0;
      }
    }
  }
  return 1;
}

static int run_chroma_pixel(void)
{
  int mx = 0, my = 0;
  for (unsigned s = 0; s < sizeof(chroma_shifts) / sizeof(chroma_shifts[0]); ++s) {
    enc.cfg.chroma_shift_w = (uint8_t)chroma_shifts[s][0];
    enc.cfg.chroma_shift_h = (uint8_t)chroma_shifts[s][1];
    for (int mode = 0; mode < 2; ++mode) {
      fill_src(mode);
      for (int wi = 0; wi < 6; ++wi) for (int hi = 0; hi < 6; ++hi) {
        for (int mxi = 0; mxi < 15; ++mxi) for (int myi = 0; myi < 15; ++myi) {
          const int w = chroma_sizes[wi], h = chroma_sizes[hi];
          const int16_t mv[2] = { (int16_t)chroma_mvs[mxi], (int16_t)chroma_mvs[myi] };
          memset(dst_ref, 0xCD, sizeof(dst_ref));
          memset(dst_opt, 0xCD, sizeof(dst_opt));
          kvz_sample_octpel_chroma_generic(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_ref, DST_STRIDE, 0, 0, mv);
          fn_chroma(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_opt, DST_STRIDE, 0, 0, mv);
          if (!compare_pixels(w, h, "chroma", &mx, &my)) return 0;
        }
      }
    }
  }
  return 1;
}

static int run_chroma_hi(void)
{
  int mx = 0, my = 0;
  for (unsigned s = 0; s < sizeof(chroma_shifts) / sizeof(chroma_shifts[0]); ++s) {
    enc.cfg.chroma_shift_w = (uint8_t)chroma_shifts[s][0];
    enc.cfg.chroma_shift_h = (uint8_t)chroma_shifts[s][1];
    for (int mode = 0; mode < 2; ++mode) {
      fill_src(mode);
      for (int wi = 0; wi < 6; ++wi) for (int hi = 0; hi < 6; ++hi) {
        for (int mxi = 0; mxi < 15; ++mxi) for (int myi = 0; myi < 15; ++myi) {
          const int w = chroma_sizes[wi], h = chroma_sizes[hi];
          const int16_t mv[2] = { (int16_t)chroma_mvs[mxi], (int16_t)chroma_mvs[myi] };
          memset(dst_ref_hi, 0xCD, sizeof(dst_ref_hi));
          memset(dst_opt_hi, 0xCD, sizeof(dst_opt_hi));
          kvz_sample_octpel_chroma_hi_generic(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_ref_hi, DST_STRIDE, 0, 0, mv);
          fn_chroma_hi(&enc, src_buf + SRC_ORIGIN, SRC_STRIDE, w, h, dst_opt_hi, DST_STRIDE, 0, 0, mv);
          if (!compare_hi(w, h, "chroma_hi", &mx, &my)) return 0;
        }
      }
    }
  }
  return 1;
}

TEST test_sample_quarterpel_luma(void)
{
  ASSERT(run_luma_pixel());
  PASS();
}

TEST test_sample_quarterpel_luma_hi(void)
{
  ASSERT(run_luma_hi());
  PASS();
}

TEST test_sample_octpel_chroma(void)
{
  ASSERT(run_chroma_pixel());
  PASS();
}

TEST test_sample_octpel_chroma_hi(void)
{
  ASSERT(run_chroma_hi());
  PASS();
}


// ---- Block filter pipeline (hpel/qpel, hor_ver and diag) ----
// The qpel and diag stages depend on intermediates from the earlier stages, so
// the pipeline runs hpel hor_ver -> hpel diag -> qpel hor_ver -> qpel diag and
// every stage's output is compared against the generic pipeline.

typedef struct {
  ipol_blocks_func *hpel_hv;
  ipol_blocks_func *hpel_diag;
  ipol_blocks_func *qpel_hv;
  ipol_blocks_func *qpel_diag;
} block_fns_t;

static block_fns_t block_ref;
static block_fns_t block_opt;

static int16_t blk_hor[5][KVZ_IPOL_MAX_IM_SIZE_LUMA_SIMD];
static int16_t blk_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1];
static kvz_pixel blk_ref[4][4][LCU_LUMA_SIZE];
static kvz_pixel blk_opt[4][4][LCU_LUMA_SIZE];

static void run_block_pipeline(const block_fns_t *f, int width, int height, int fme_level,
                               int hx, int hy, kvz_pixel stages[4][4][LCU_LUMA_SIZE])
{
  memset(blk_hor, 0, sizeof(blk_hor));
  memset(blk_cols, 0, sizeof(blk_cols));
  memset(stages, 0xEE, sizeof(blk_ref));
  kvz_pixel *src = src_buf + SRC_ORIGIN;
  f->hpel_hv(&enc, src, SRC_STRIDE, width, height, stages[0], blk_hor, (int8_t)fme_level, blk_cols, (int8_t)hx, (int8_t)hy);
  f->hpel_diag(&enc, src, SRC_STRIDE, width, height, stages[1], blk_hor, (int8_t)fme_level, blk_cols, (int8_t)hx, (int8_t)hy);
  f->qpel_hv(&enc, src, SRC_STRIDE, width, height, stages[2], blk_hor, (int8_t)fme_level, blk_cols, (int8_t)hx, (int8_t)hy);
  f->qpel_diag(&enc, src, SRC_STRIDE, width, height, stages[3], blk_hor, (int8_t)fme_level, blk_cols, (int8_t)hx, (int8_t)hy);
}

static int run_block_filters(void)
{
  // The encoder rounds internal_width and internal_height up to multiples of 8
  // before calling the block filters (search_inter.c), so the block contract
  // is multiples of 8 only.
  static const int bsizes[] = { 8, 16, 32, 64 };
  static const int offs[] = { -1, 0, 1 };
  static const char *names[4] = { "hpel_hor_ver", "hpel_diag", "qpel_hor_ver", "qpel_diag" };
  for (int mode = 0; mode < 2; ++mode) {
    fill_src(mode);
    for (int wi = 0; wi < 4; ++wi) for (int hi = 0; hi < 4; ++hi) {
      for (int fl = 1; fl <= 3; ++fl) for (int ox = 0; ox < 3; ++ox) for (int oy = 0; oy < 3; ++oy) {
        const int w = bsizes[wi], h = bsizes[hi];
        run_block_pipeline(&block_ref, w, h, fl, offs[ox], offs[oy], blk_ref);
        run_block_pipeline(&block_opt, w, h, fl, offs[ox], offs[oy], blk_opt);
        // The encoder runs only the first fme_level steps (search_inter.c), so
        // only those stages are part of the contract.
        for (int st = 0; st < fl; ++st) for (int fi = 0; fi < 4; ++fi) {
          for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
            if (blk_ref[st][fi][y * LCU_WIDTH + x] != blk_opt[st][fi][y * LCU_WIDTH + x]) {
              printf("\n    %s filter %d mismatch at (%d,%d): ref=%d opt=%d (w=%d h=%d fme=%d hx=%d hy=%d)\n",
                     names[st], fi, x, y, blk_ref[st][fi][y * LCU_WIDTH + x],
                     blk_opt[st][fi][y * LCU_WIDTH + x], w, h, fl, offs[ox], offs[oy]);
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

TEST test_block_filters(void)
{
  ASSERT(run_block_filters());
  PASS();
}

static void *find_generic_ipol(const char *type)
{
  for (int i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, "generic") == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

SUITE(ipol_tests)
{
  block_ref.hpel_hv = find_generic_ipol("filter_hpel_blocks_hor_ver_luma");
  block_ref.hpel_diag = find_generic_ipol("filter_hpel_blocks_diag_luma");
  block_ref.qpel_hv = find_generic_ipol("filter_qpel_blocks_hor_ver_luma");
  block_ref.qpel_diag = find_generic_ipol("filter_qpel_blocks_diag_luma");

  memset(&enc, 0, sizeof(enc));
  fill_src(0);

  for (volatile int i = 0; i < strategies.count; ++i) {
    const strategy_t *s = &strategies.strategies[i];
    if (strcmp(s->strategy_name, "generic") == 0) {
      continue;
    }
    if (strcmp(s->type, "sample_quarterpel_luma") == 0) {
      fn_luma = s->fptr;
      RUN_TEST(test_sample_quarterpel_luma);
    } else if (strcmp(s->type, "sample_quarterpel_luma_hi") == 0) {
      fn_luma_hi = s->fptr;
      RUN_TEST(test_sample_quarterpel_luma_hi);
    } else if (strcmp(s->type, "sample_octpel_chroma") == 0) {
      fn_chroma = s->fptr;
      RUN_TEST(test_sample_octpel_chroma);
    } else if (strcmp(s->type, "sample_octpel_chroma_hi") == 0) {
      fn_chroma_hi = s->fptr;
      RUN_TEST(test_sample_octpel_chroma_hi);
    }
  }

  // Block filters are tested as one pipeline, so run it once all four
  // optimized versions are present.
  if (block_ref.hpel_hv && block_ref.hpel_diag && block_ref.qpel_hv && block_ref.qpel_diag) {
    for (volatile int i = 0; i < strategies.count; ++i) {
      const strategy_t *s = &strategies.strategies[i];
      if (strcmp(s->strategy_name, "generic") == 0) {
        continue;
      }
      if (strcmp(s->type, "filter_hpel_blocks_hor_ver_luma") == 0) block_opt.hpel_hv = s->fptr;
      if (strcmp(s->type, "filter_hpel_blocks_diag_luma") == 0) block_opt.hpel_diag = s->fptr;
      if (strcmp(s->type, "filter_qpel_blocks_hor_ver_luma") == 0) block_opt.qpel_hv = s->fptr;
      if (strcmp(s->type, "filter_qpel_blocks_diag_luma") == 0) block_opt.qpel_diag = s->fptr;
    }
    if (block_opt.hpel_hv && block_opt.hpel_diag && block_opt.qpel_hv && block_opt.qpel_diag) {
      RUN_TEST(test_block_filters);
    }
  }
}
