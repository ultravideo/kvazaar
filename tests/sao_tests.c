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
#include "encoderstate.h"
#include "sao.h"
#include "strategies/strategies-sao.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define MAX_DIM 64
// Reconstruction buffers have one row and column of margin before the block
// because the edge classes read neighbours at -1.
#define STRIDE 80
#define ORIG_OFF (2 * STRIDE + 2)

static kvz_pixel rec_c[MAX_DIM * MAX_DIM];
static kvz_pixel orig_c[MAX_DIM * MAX_DIM];
static kvz_pixel rec_buf[STRIDE * 80];
static kvz_pixel out_ref_buf[STRIDE * 80];
static kvz_pixel out_opt_buf[STRIDE * 80];

static encoder_control_t enc;
static encoder_state_t st;

static sao_edge_ddistortion_func *fn_ref_edge_dd = NULL;
static sao_edge_ddistortion_func *fn_opt_edge_dd = NULL;
static calc_sao_edge_dir_func *fn_ref_edge_dir = NULL;
static calc_sao_edge_dir_func *fn_opt_edge_dir = NULL;
static sao_reconstruct_color_func *fn_ref_recon = NULL;
static sao_reconstruct_color_func *fn_opt_recon = NULL;
static sao_band_ddistortion_func *fn_ref_band_dd = NULL;
static sao_band_ddistortion_func *fn_opt_band_dd = NULL;

static uint32_t rng_state = 0xC0FFEEu;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

static void fill_pixels(kvz_pixel *p, int n, int mode)
{
  for (int i = 0; i < n; ++i) {
    if (mode == 0) {
      p[i] = (kvz_pixel)(rnd() & 0xFF);
    } else {
      p[i] = (rnd() & 1) ? 255 : 0;
    }
  }
}

// Offsets in the SAO range for 8-bit video: [-7, 7].
static int rnd_offset(void)
{
  return (int)(rnd() % 15) - 7;
}

static const int dims[] = { 4, 8, 12, 16, 32, 64 };
#define NUM_DIMS ((int)(sizeof(dims) / sizeof(dims[0])))

static int run_edge_dd(void)
{
  for (int mode = 0; mode < 2; ++mode) {
    for (int trial = 0; trial < 6; ++trial) {
      fill_pixels(rec_c, MAX_DIM * MAX_DIM, mode);
      fill_pixels(orig_c, MAX_DIM * MAX_DIM, mode);
      for (int eo = 0; eo < 4; ++eo) {
        for (int wi = 0; wi < NUM_DIMS; ++wi) for (int hi = 0; hi < NUM_DIMS; ++hi) {
          int offsets[NUM_SAO_EDGE_CATEGORIES];
          for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES; ++k) offsets[k] = rnd_offset();
          if (trial == 0) for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES; ++k) offsets[k] = 0;
          const int w = dims[wi], h = dims[hi];
          int ref = fn_ref_edge_dd(&enc, orig_c, rec_c, w, h, eo, offsets);
          int opt = fn_opt_edge_dd(&enc, orig_c, rec_c, w, h, eo, offsets);
          if (ref != opt) {
            printf("\n    edge_ddistortion mismatch eo=%d w=%d h=%d: ref=%d opt=%d\n", eo, w, h, ref, opt);
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

static int run_edge_dir(void)
{
  for (int mode = 0; mode < 2; ++mode) {
    for (int trial = 0; trial < 6; ++trial) {
      fill_pixels(rec_c, MAX_DIM * MAX_DIM, mode);
      fill_pixels(orig_c, MAX_DIM * MAX_DIM, mode);
      for (int eo = 0; eo < 4; ++eo) {
        for (int wi = 0; wi < NUM_DIMS; ++wi) for (int hi = 0; hi < NUM_DIMS; ++hi) {
          const int w = dims[wi], h = dims[hi];
          int ref[2][NUM_SAO_EDGE_CATEGORIES];
          int opt[2][NUM_SAO_EDGE_CATEGORIES];
          memset(ref, 0, sizeof(ref));
          memset(opt, 0, sizeof(opt));
          fn_ref_edge_dir(&enc, orig_c, rec_c, eo, w, h, ref);
          fn_opt_edge_dir(&enc, orig_c, rec_c, eo, w, h, opt);
          if (memcmp(ref, opt, sizeof(ref)) != 0) {
            printf("\n    calc_sao_edge_dir mismatch eo=%d w=%d h=%d\n", eo, w, h);
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

static int run_band_dd(void)
{
  for (int mode = 0; mode < 2; ++mode) {
    for (int trial = 0; trial < 6; ++trial) {
      fill_pixels(rec_c, MAX_DIM * MAX_DIM, mode);
      fill_pixels(orig_c, MAX_DIM * MAX_DIM, mode);
      for (int band_pos = 0; band_pos <= 28; band_pos += 4) {
        for (int wi = 0; wi < NUM_DIMS; ++wi) for (int hi = 0; hi < NUM_DIMS; ++hi) {
          int sao_bands[4];
          for (int k = 0; k < 4; ++k) sao_bands[k] = rnd_offset();
          const int w = dims[wi], h = dims[hi];
          int ref = fn_ref_band_dd(&st, orig_c, rec_c, w, h, band_pos, sao_bands);
          int opt = fn_opt_band_dd(&st, orig_c, rec_c, w, h, band_pos, sao_bands);
          if (ref != opt) {
            printf("\n    band_ddistortion mismatch band_pos=%d w=%d h=%d: ref=%d opt=%d\n", band_pos, w, h, ref, opt);
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

static int run_recon(void)
{
  static const color_t colors[3] = { COLOR_Y, COLOR_U, COLOR_V };
  for (int mode = 0; mode < 2; ++mode) {
    for (int trial = 0; trial < 6; ++trial) {
      fill_pixels(rec_buf, STRIDE * 80, mode);
      for (int type = SAO_TYPE_BAND; type <= SAO_TYPE_EDGE; ++type) {
        for (int ci = 0; ci < 3; ++ci) {
          for (int wi = 0; wi < NUM_DIMS; ++wi) for (int hi = 0; hi < NUM_DIMS; ++hi) {
            sao_info_t sao;
            memset(&sao, 0, sizeof(sao));
            sao.type = (sao_type)type;
            sao.eo_class = (sao_eo_class)(rnd() % 4);
            sao.band_position[0] = (int)(rnd() % 29);
            sao.band_position[1] = (int)(rnd() % 29);
            for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES * 2; ++k) sao.offsets[k] = rnd_offset();
            const int w = dims[wi], h = dims[hi];
            memset(out_ref_buf, 0x5A, sizeof(out_ref_buf));
            memset(out_opt_buf, 0x5A, sizeof(out_opt_buf));
            kvz_pixel *rec = rec_buf + ORIG_OFF;
            kvz_pixel *out_ref = out_ref_buf + ORIG_OFF;
            kvz_pixel *out_opt = out_opt_buf + ORIG_OFF;
            fn_ref_recon(&enc, rec, out_ref, &sao, STRIDE, STRIDE, w, h, colors[ci]);
            fn_opt_recon(&enc, rec, out_opt, &sao, STRIDE, STRIDE, w, h, colors[ci]);
            for (int y = 0; y < h; ++y) {
              if (memcmp(out_ref + y * STRIDE, out_opt + y * STRIDE, (size_t)w) != 0) {
                printf("\n    sao_reconstruct_color mismatch type=%d color=%d w=%d h=%d row=%d\n", type, ci, w, h, y);
                return 0;
              }
            }
          }
        }
      }
    }
  }
  return 1;
}

TEST test_edge_ddistortion(void)
{
  ASSERT(run_edge_dd());
  PASS();
}

TEST test_calc_sao_edge_dir(void)
{
  ASSERT(run_edge_dir());
  PASS();
}

TEST test_band_ddistortion(void)
{
  ASSERT(run_band_dd());
  PASS();
}

TEST test_sao_reconstruct_color(void)
{
  ASSERT(run_recon());
  PASS();
}

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

SUITE(sao_tests)
{
  memset(&enc, 0, sizeof(enc));
  enc.bitdepth = 8;
  memset(&st, 0, sizeof(st));
  st.encoder_control = &enc;

  fn_ref_edge_dd = find_generic("sao_edge_ddistortion");
  fn_ref_edge_dir = find_generic("calc_sao_edge_dir");
  fn_ref_recon = find_generic("sao_reconstruct_color");
  fn_ref_band_dd = find_generic("sao_band_ddistortion");

  for (volatile int i = 0; i < strategies.count; ++i) {
    const strategy_t *s = &strategies.strategies[i];
    if (strcmp(s->strategy_name, "generic") == 0) {
      continue;
    }
    if (strcmp(s->type, "sao_edge_ddistortion") == 0 && fn_ref_edge_dd) {
      fn_opt_edge_dd = s->fptr;
      RUN_TEST(test_edge_ddistortion);
    } else if (strcmp(s->type, "calc_sao_edge_dir") == 0 && fn_ref_edge_dir) {
      fn_opt_edge_dir = s->fptr;
      RUN_TEST(test_calc_sao_edge_dir);
    } else if (strcmp(s->type, "sao_reconstruct_color") == 0 && fn_ref_recon) {
      fn_opt_recon = s->fptr;
      RUN_TEST(test_sao_reconstruct_color);
    } else if (strcmp(s->type, "sao_band_ddistortion") == 0 && fn_ref_band_dd) {
      fn_opt_band_dd = s->fptr;
      RUN_TEST(test_band_ddistortion);
    }
  }
}
