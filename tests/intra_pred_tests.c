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

#include "strategies/strategies-intra.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

// Reference arrays use index 0 as the top-left corner, so the block code's
// "-1" index maps to offset 1 from the start of the buffer. Pointers passed to
// the functions are placed at offset 64 so there is room for any read the
// generic code might make on either side.
#define REF_OFFSET 64
#define REF_LEN 256

static kvz_pixel ref_above[REF_LEN];
static kvz_pixel ref_left[REF_LEN];
static kvz_pixel out_ref[32 * 32];
static kvz_pixel out_opt[32 * 32];

static uint32_t rng_state = 0x9E3779B9u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

static void fill_refs(int mode)
{
  for (int i = 0; i < REF_LEN; ++i) {
    if (mode == 0) {
      ref_above[i] = (kvz_pixel)(rnd() & 0xFF);
      ref_left[i] = (kvz_pixel)(rnd() & 0xFF);
    } else {
      // Extremes stress the rounding and the clip-free 16-bit intermediates.
      ref_above[i] = (rnd() & 1) ? 255 : 0;
      ref_left[i] = (rnd() & 1) ? 255 : 0;
    }
  }
}

static angular_pred_func *fn_ref_angular = NULL;
static angular_pred_func *fn_opt_angular = NULL;
static intra_pred_planar_func *fn_ref_planar = NULL;
static intra_pred_planar_func *fn_opt_planar = NULL;
static intra_pred_filtered_dc_func *fn_ref_dc = NULL;
static intra_pred_filtered_dc_func *fn_opt_dc = NULL;

static int compare_block(int width, const char *what, int mode, int log2)
{
  if (memcmp(out_ref, out_opt, (size_t)(width * width)) == 0) {
    return 1;
  }
  for (int i = 0; i < width * width; ++i) {
    if (out_ref[i] != out_opt[i]) {
      printf("\n    %s mismatch mode=%d log2=%d at (%d,%d): ref=%d opt=%d\n",
             what, mode, log2, i % width, i / width, out_ref[i], out_opt[i]);
      break;
    }
  }
  return 0;
}

static int run_angular(void)
{
  for (int mode_fill = 0; mode_fill < 2; ++mode_fill) {
    for (int trial = 0; trial < 4; ++trial) {
      fill_refs(mode_fill);
      for (int log2 = 2; log2 <= 5; ++log2) {
        for (int intra_mode = 2; intra_mode <= 34; ++intra_mode) {
          memset(out_ref, 0xAB, sizeof(out_ref));
          memset(out_opt, 0xAB, sizeof(out_opt));
          fn_ref_angular((int_fast8_t)log2, (int_fast8_t)intra_mode,
                         ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_ref);
          fn_opt_angular((int_fast8_t)log2, (int_fast8_t)intra_mode,
                         ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_opt);
          if (!compare_block(1 << log2, "angular", intra_mode, log2)) return 0;
        }
      }
    }
  }
  return 1;
}

static int run_planar(void)
{
  for (int mode_fill = 0; mode_fill < 2; ++mode_fill) {
    for (int trial = 0; trial < 8; ++trial) {
      fill_refs(mode_fill);
      for (int log2 = 2; log2 <= 5; ++log2) {
        memset(out_ref, 0xAB, sizeof(out_ref));
        memset(out_opt, 0xAB, sizeof(out_opt));
        fn_ref_planar((int_fast8_t)log2, ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_ref);
        fn_opt_planar((int_fast8_t)log2, ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_opt);
        if (!compare_block(1 << log2, "planar", -1, log2)) return 0;
      }
    }
  }
  return 1;
}

static int run_dc(void)
{
  for (int mode_fill = 0; mode_fill < 2; ++mode_fill) {
    for (int trial = 0; trial < 8; ++trial) {
      fill_refs(mode_fill);
      for (int log2 = 2; log2 <= 5; ++log2) {
        memset(out_ref, 0xAB, sizeof(out_ref));
        memset(out_opt, 0xAB, sizeof(out_opt));
        fn_ref_dc((int_fast8_t)log2, ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_ref);
        fn_opt_dc((int_fast8_t)log2, ref_above + REF_OFFSET, ref_left + REF_OFFSET, out_opt);
        if (!compare_block(1 << log2, "filtered_dc", -1, log2)) return 0;
      }
    }
  }
  return 1;
}

TEST test_angular_pred(void)
{
  ASSERT(run_angular());
  PASS();
}

TEST test_intra_pred_planar(void)
{
  ASSERT(run_planar());
  PASS();
}

TEST test_intra_pred_filtered_dc(void)
{
  ASSERT(run_dc());
  PASS();
}

// Finds the generic reference for a strategy type. Returns NULL if absent.
static void *find_generic(const char *type)
{
  for (int i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, "generic") == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

SUITE(intra_pred_tests)
{
  fn_ref_angular = find_generic("angular_pred");
  fn_ref_planar = find_generic("intra_pred_planar");
  fn_ref_dc = find_generic("intra_pred_filtered_dc");

  for (volatile int i = 0; i < strategies.count; ++i) {
    const strategy_t *s = &strategies.strategies[i];
    if (strcmp(s->strategy_name, "generic") == 0) {
      continue;
    }
    if (strcmp(s->type, "angular_pred") == 0 && fn_ref_angular) {
      fn_opt_angular = s->fptr;
      RUN_TEST(test_angular_pred);
    } else if (strcmp(s->type, "intra_pred_planar") == 0 && fn_ref_planar) {
      fn_opt_planar = s->fptr;
      RUN_TEST(test_intra_pred_planar);
    } else if (strcmp(s->type, "intra_pred_filtered_dc") == 0 && fn_ref_dc) {
      fn_opt_dc = s->fptr;
      RUN_TEST(test_intra_pred_filtered_dc);
    }
  }
}
