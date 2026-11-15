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
 * Differential tests for the intra prediction and SAO strategies at the bit
 * depth the tests are compiled for (8, 10 or 12): every non-generic
 * implementation is compared with the generic one on random, extreme and
 * edge-case inputs, for all block sizes, intra modes, SAO classes, band
 * positions and offsets.
 *
 * Values outside the range a function may read (the reference arrays beyond
 * index 2 * width, the samples around a SAO block) are filled with different
 * garbage for two runs of the optimized function, so a read outside the range
 * shows up as a mismatch if it affects the result.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "encoder.h"
#include "encoderstate.h"
#include "sao.h"
#include "strategies/strategies-intra.h"
#include "strategies/strategies-sao.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

static uint32_t rng_state = 0x6A09E667u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

static kvz_pixel rnd_pixel(void)
{
  return (kvz_pixel)(rnd() & PIXEL_MAX);
}

// Pixel patterns. 0: random, 1: random 0/PIXEL_MAX, 2: all PIXEL_MAX,
// 3: all 0, 4: random near PIXEL_MAX, 5: random near 0, 6: alternating
// 0/PIXEL_MAX, 7: random in a narrow range (many equal neighbours).
#define NUM_PATTERNS 8
static kvz_pixel pattern_pixel(int pattern, int i)
{
  switch (pattern) {
    case 0: return rnd_pixel();
    case 1: return (rnd() & 1) ? PIXEL_MAX : 0;
    case 2: return PIXEL_MAX;
    case 3: return 0;
    case 4: return (kvz_pixel)(PIXEL_MAX - (rnd() & 7));
    case 5: return (kvz_pixel)(rnd() & 7);
    case 6: return (i & 1) ? PIXEL_MAX : 0;
    default: return (kvz_pixel)((PIXEL_MAX >> 1) + (rnd() % 3));
  }
}

static const strategy_t *cur;

static void *find_generic(const char *type)
{
  for (unsigned i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, type) == 0 &&
        strcmp(strategies.strategies[i].strategy_name, "generic") == 0) {
      return strategies.strategies[i].fptr;
    }
  }
  return NULL;
}

/*
 * Intra prediction
 */

// The functions get pointers to index -1 of references of length
// 2 * width + 1, like kvz_intra_ref. GUARD values on either side are garbage.
#define GUARD 64
#define REF_LEN (2 * 32 + 1)
static kvz_pixel ref_top_buf[GUARD + REF_LEN + GUARD];
static kvz_pixel ref_left_buf[GUARD + REF_LEN + GUARD];
static kvz_pixel out_ref[32 * 32];
static kvz_pixel out_opt[32 * 32];

static void fill_intra_refs(int pattern_top, int pattern_left, int width)
{
  for (int i = 0; i < GUARD + REF_LEN + GUARD; ++i) {
    ref_top_buf[i] = rnd_pixel();
    ref_left_buf[i] = rnd_pixel();
  }
  for (int i = 0; i < 2 * width + 1; ++i) {
    ref_top_buf[GUARD + i] = pattern_pixel(pattern_top, i);
    ref_left_buf[GUARD + i] = pattern_pixel(pattern_left, i);
  }
  // Index -1 (top-left corner) is shared by both references.
  ref_left_buf[GUARD] = ref_top_buf[GUARD];
}

// Changes the values the functions must not read.
static void scramble_intra_guards(int width)
{
  for (int i = 0; i < GUARD; ++i) {
    ref_top_buf[i] = rnd_pixel();
    ref_left_buf[i] = rnd_pixel();
  }
  for (int i = GUARD + 2 * width + 1; i < GUARD + REF_LEN + GUARD; ++i) {
    ref_top_buf[i] = rnd_pixel();
    ref_left_buf[i] = rnd_pixel();
  }
}

static int compare_blocks(int width, const char *what, int mode, int pt, int pl)
{
  for (int i = 0; i < width * width; ++i) {
    if (out_ref[i] != out_opt[i]) {
      printf("\n    %s:%s mismatch mode=%d width=%d patterns=%d/%d at (%d,%d): ref=%d opt=%d\n",
             what, cur->strategy_name, mode, width, pt, pl, i % width, i / width,
             out_ref[i], out_opt[i]);
      return 0;
    }
  }
  return 1;
}

typedef void (intra_any_func)(int_fast8_t log2_width, int_fast8_t mode,
                              const kvz_pixel *top, const kvz_pixel *left,
                              kvz_pixel *dst, void *fn);

static void call_angular(int_fast8_t log2_width, int_fast8_t mode,
                         const kvz_pixel *top, const kvz_pixel *left,
                         kvz_pixel *dst, void *fn)
{
  ((angular_pred_func *)fn)(log2_width, mode, top, left, dst);
}

static void call_planar(int_fast8_t log2_width, int_fast8_t mode,
                        const kvz_pixel *top, const kvz_pixel *left,
                        kvz_pixel *dst, void *fn)
{
  (void)mode;
  ((intra_pred_planar_func *)fn)(log2_width, top, left, dst);
}

static void call_dc(int_fast8_t log2_width, int_fast8_t mode,
                    const kvz_pixel *top, const kvz_pixel *left,
                    kvz_pixel *dst, void *fn)
{
  (void)mode;
  ((intra_pred_filtered_dc_func *)fn)(log2_width, top, left, dst);
}

static int run_intra(intra_any_func *call, void *ref_fn, void *opt_fn,
                     int first_mode, int last_mode, const char *what)
{
  for (int log2 = 2; log2 <= 5; ++log2) {
    const int width = 1 << log2;
    for (int pt = 0; pt < NUM_PATTERNS; ++pt) {
      for (int pl = 0; pl < NUM_PATTERNS; ++pl) {
        // Random patterns are repeated a few times.
        const int trials = (pt == 0 || pl == 0) ? 4 : 1;
        for (int trial = 0; trial < trials; ++trial) {
          fill_intra_refs(pt, pl, width);
          for (int mode = first_mode; mode <= last_mode; ++mode) {
            const kvz_pixel *top = ref_top_buf + GUARD;
            const kvz_pixel *left = ref_left_buf + GUARD;
            memset(out_ref, 0xAB, sizeof(out_ref));
            call((int_fast8_t)log2, (int_fast8_t)mode, top, left, out_ref, ref_fn);
            for (int run = 0; run < 2; ++run) {
              memset(out_opt, 0xCD, sizeof(out_opt));
              call((int_fast8_t)log2, (int_fast8_t)mode, top, left, out_opt, opt_fn);
              if (!compare_blocks(width, what, mode, pt, pl)) return 0;
              scramble_intra_guards(width);
            }
          }
        }
      }
    }
  }
  return 1;
}

static void *intra_ref_fn;
static void *intra_opt_fn;

TEST test_angular_pred(void)
{
  ASSERT(run_intra(call_angular, intra_ref_fn, intra_opt_fn, 2, 34, "angular_pred"));
  PASS();
}

TEST test_intra_pred_planar(void)
{
  ASSERT(run_intra(call_planar, intra_ref_fn, intra_opt_fn, 0, 0, "intra_pred_planar"));
  PASS();
}

TEST test_intra_pred_filtered_dc(void)
{
  ASSERT(run_intra(call_dc, intra_ref_fn, intra_opt_fn, 1, 1, "intra_pred_filtered_dc"));
  PASS();
}

/*
 * SAO
 */

// Contiguous blocks (block_width stride) as passed by the SAO search.
#define MAX_DIM 64
static kvz_pixel orig_blk[MAX_DIM * MAX_DIM + 64];
static kvz_pixel rec_blk[MAX_DIM * MAX_DIM + 64];

// Strided buffers for reconstruction, with a margin of two samples around the
// block.
#define REC_STRIDE 80
#define REC_ROWS (MAX_DIM + 4)
#define REC_ORIGIN (2 * REC_STRIDE + 2)
static kvz_pixel recon_in[REC_STRIDE * REC_ROWS];
static kvz_pixel recon_out_ref[REC_STRIDE * REC_ROWS];
static kvz_pixel recon_out_opt[REC_STRIDE * REC_ROWS];

static encoder_control_t enc;
static encoder_state_t st;

static void *sao_ref_fn;
static void *sao_opt_fn;

static const int sao_widths[] = {
  1, 2, 3, 4, 5, 7, 8, 9, 10, 15, 16, 17, 18, 19, 23, 24, 25, 30, 31, 32, 33,
  34, 35, 40, 47, 48, 49, 50, 56, 62, 63, 64
};
static const int sao_heights[] = { 1, 2, 3, 4, 5, 8, 15, 16, 31, 32, 33, 62, 63, 64 };
#define NUM_SAO_WIDTHS ((int)(sizeof(sao_widths) / sizeof(sao_widths[0])))
#define NUM_SAO_HEIGHTS ((int)(sizeof(sao_heights) / sizeof(sao_heights[0])))

// The SAO search works on whole LCUs, which are cut at the frame edge to a
// multiple of the minimum CU size (8 luma, 4 chroma samples). The 10 and
// 12-bit implementations are tested with any size, the 8-bit ones only with
// the sizes the encoder uses.
static int search_size_ok(int w, int h)
{
  return KVZ_BIT_DEPTH > 8 || (w % 4 == 0 && h % 4 == 0);
}

// Offset sets. 0: random in the valid range, 1: all +max, 2: all -max,
// 3: zeros, 4: random sparse, 5: beyond the valid range. The last one checks
// the robustness of the 10 and 12-bit implementations; the 8-bit ones are
// only correct for valid offsets. Kind 5 stays within +-48: the int32 sums of
// the generic ddistortion functions overflow with larger offsets on 64x64
// blocks at 12-bit (|offset| * (|offset| + 2 * 4095) * 4096 < 2^31).
#if KVZ_BIT_DEPTH > 8
#define NUM_OFFSET_KINDS 6
#else
#define NUM_OFFSET_KINDS 5
#endif
static int rnd_offset(int kind)
{
  const int max = SAO_ABS_OFFSET_MAX;
  switch (kind) {
    case 0: return (int)(rnd() % (2 * max + 1)) - max;
    case 1: return max;
    case 2: return -max;
    case 3: return 0;
    case 4: return (rnd() & 3) ? 0 : (int)(rnd() % (2 * max + 1)) - max;
    default: return (int)(rnd() % 97) - 48;
  }
}

static void fill_sao_blocks(int pattern_orig, int pattern_rec)
{
  for (int i = 0; i < MAX_DIM * MAX_DIM + 64; ++i) {
    orig_blk[i] = pattern_pixel(pattern_orig, i);
    rec_blk[i] = pattern_pixel(pattern_rec, i);
  }
}

static int run_edge_ddistortion(void)
{
  sao_edge_ddistortion_func *ref = sao_ref_fn, *opt = sao_opt_fn;
  for (int po = 0; po < NUM_PATTERNS; ++po) {
    for (int pr = 0; pr < NUM_PATTERNS; ++pr) {
      fill_sao_blocks(po, pr);
      for (int eo = 0; eo < 4; ++eo) {
        for (int wi = 0; wi < NUM_SAO_WIDTHS; ++wi) {
          for (int hi = 0; hi < NUM_SAO_HEIGHTS; ++hi) {
            const int w = sao_widths[wi], h = sao_heights[hi];
            if (!search_size_ok(w, h)) continue;
            const int kind = (int)(rnd() % NUM_OFFSET_KINDS);
            int offsets[NUM_SAO_EDGE_CATEGORIES];
            for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES; ++k) offsets[k] = rnd_offset(kind);
            const int r = ref(&enc, orig_blk, rec_blk, w, h, eo, offsets);
            const int o = opt(&enc, orig_blk, rec_blk, w, h, eo, offsets);
            if (r != o) {
              printf("\n    sao_edge_ddistortion:%s mismatch eo=%d w=%d h=%d patterns=%d/%d kind=%d: ref=%d opt=%d\n",
                     cur->strategy_name, eo, w, h, po, pr, kind, r, o);
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

static int run_calc_sao_edge_dir(void)
{
  calc_sao_edge_dir_func *ref = sao_ref_fn, *opt = sao_opt_fn;
  for (int po = 0; po < NUM_PATTERNS; ++po) {
    for (int pr = 0; pr < NUM_PATTERNS; ++pr) {
      fill_sao_blocks(po, pr);
      for (int eo = 0; eo < 4; ++eo) {
        for (int wi = 0; wi < NUM_SAO_WIDTHS; ++wi) {
          for (int hi = 0; hi < NUM_SAO_HEIGHTS; ++hi) {
            const int w = sao_widths[wi], h = sao_heights[hi];
            if (!search_size_ok(w, h)) continue;
            int r[2][NUM_SAO_EDGE_CATEGORIES], o[2][NUM_SAO_EDGE_CATEGORIES];
            // The functions accumulate into the array.
            for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES; ++k) {
              r[0][k] = o[0][k] = (int)(rnd() % 1000) - 500;
              r[1][k] = o[1][k] = (int)(rnd() % 1000);
            }
            ref(&enc, orig_blk, rec_blk, eo, w, h, r);
            opt(&enc, orig_blk, rec_blk, eo, w, h, o);
            if (memcmp(r, o, sizeof(r)) != 0) {
              printf("\n    calc_sao_edge_dir:%s mismatch eo=%d w=%d h=%d patterns=%d/%d\n",
                     cur->strategy_name, eo, w, h, po, pr);
              for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES; ++k) {
                printf("      cat %d: ref sum=%d cnt=%d, opt sum=%d cnt=%d\n",
                       k, r[0][k], r[1][k], o[0][k], o[1][k]);
              }
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

static int run_band_ddistortion(void)
{
  sao_band_ddistortion_func *ref = sao_ref_fn, *opt = sao_opt_fn;
  for (int po = 0; po < NUM_PATTERNS; ++po) {
    for (int pr = 0; pr < NUM_PATTERNS; ++pr) {
      fill_sao_blocks(po, pr);
      for (int band_pos = 0; band_pos < 32; ++band_pos) {
        for (int wi = 0; wi < NUM_SAO_WIDTHS; ++wi) {
          // Only a subset of heights per band position to keep the test fast.
          for (int hi = band_pos % 3; hi < NUM_SAO_HEIGHTS; hi += 3) {
            const int w = sao_widths[wi], h = sao_heights[hi];
            if (!search_size_ok(w, h)) continue;
            const int kind = (int)(rnd() % NUM_OFFSET_KINDS);
            int bands[4];
            for (int k = 0; k < 4; ++k) bands[k] = rnd_offset(kind);
            const int r = ref(&st, orig_blk, rec_blk, w, h, band_pos, bands);
            const int o = opt(&st, orig_blk, rec_blk, w, h, band_pos, bands);
            if (r != o) {
              printf("\n    sao_band_ddistortion:%s mismatch band_pos=%d w=%d h=%d patterns=%d/%d kind=%d: ref=%d opt=%d\n",
                     cur->strategy_name, band_pos, w, h, po, pr, kind, r, o);
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

static int run_sao_reconstruct_color(void)
{
  sao_reconstruct_color_func *ref = sao_ref_fn, *opt = sao_opt_fn;
  static const color_t colors[3] = { COLOR_Y, COLOR_U, COLOR_V };
  for (int pattern = 0; pattern < NUM_PATTERNS; ++pattern) {
    for (int i = 0; i < REC_STRIDE * REC_ROWS; ++i) {
      recon_in[i] = pattern_pixel(pattern, i);
    }
    for (int type = SAO_TYPE_BAND; type <= SAO_TYPE_EDGE; ++type) {
      for (int ci = 0; ci < 3; ++ci) {
        for (int wi = 0; wi < NUM_SAO_WIDTHS; ++wi) {
          for (int hi = 0; hi < NUM_SAO_HEIGHTS; ++hi) {
            const int w = sao_widths[wi], h = sao_heights[hi];
            sao_info_t sao;
            memset(&sao, 0, sizeof(sao));
            sao.type = (sao_type)type;
            sao.eo_class = (sao_eo_class)(rnd() % 4);
            sao.band_position[0] = (int)(rnd() % 32);
            sao.band_position[1] = (int)(rnd() % 32);
            const int kind = (int)(rnd() % (NUM_OFFSET_KINDS + (KVZ_BIT_DEPTH > 8)));
            for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES * 2; ++k) {
              // Kind 6 (10 and 12-bit): offsets that exceed the sample range.
              sao.offsets[k] = kind < NUM_OFFSET_KINDS ? rnd_offset(kind)
                                                      : (int)(rnd() % 200001) - 100000;
            }
            // Output stride different from the input stride.
            const int out_stride = REC_STRIDE - 8;
            memset(recon_out_ref, 0x5A, sizeof(recon_out_ref));
            memset(recon_out_opt, 0x5A, sizeof(recon_out_opt));
            ref(&enc, recon_in + REC_ORIGIN, recon_out_ref + REC_ORIGIN, &sao,
                REC_STRIDE, out_stride, w, h, colors[ci]);
            opt(&enc, recon_in + REC_ORIGIN, recon_out_opt + REC_ORIGIN, &sao,
                REC_STRIDE, out_stride, w, h, colors[ci]);
            // The whole buffer must match: nothing outside the block is written.
            if (memcmp(recon_out_ref, recon_out_opt, sizeof(recon_out_ref)) != 0) {
              int pos = 0;
              while (recon_out_ref[pos] == recon_out_opt[pos]) ++pos;
              pos -= REC_ORIGIN;
              printf("\n    sao_reconstruct_color:%s mismatch type=%d color=%d w=%d h=%d eo=%d kind=%d pattern=%d at (%d,%d)\n",
                     cur->strategy_name, type, ci, w, h, sao.eo_class, kind, pattern,
                     pos % out_stride, pos / out_stride);
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

// The original kvz_calc_sao_offset_array, which went through every value.
static void sao_offset_array_reference(const sao_info_t *sao, int *offset, color_t color_i)
{
  const int32_t values = 1 << KVZ_BIT_DEPTH;
  const int32_t shift = KVZ_BIT_DEPTH - 5;
  const int32_t band_pos = (color_i == COLOR_V) ? 1 : 0;
  const int32_t cur_bp = sao->band_position[band_pos];
  for (int32_t val = 0; val < values; val++) {
    const int32_t cb_minus_cbp = (val >> shift) - cur_bp;
    if (cb_minus_cbp >= 0 && cb_minus_cbp <= 3) {
      const int32_t val_unclipped = val + sao->offsets[cb_minus_cbp + 1 + 5 * band_pos];
      offset[val] = CLIP(0, values - 1, val_unclipped);
    } else {
      offset[val] = val;
    }
  }
}

TEST test_sao_offset_array(void)
{
  static int ref[1 << KVZ_BIT_DEPTH];
  static int out[1 << KVZ_BIT_DEPTH];
  static const color_t colors[3] = { COLOR_Y, COLOR_U, COLOR_V };
  for (int ci = 0; ci < 3; ++ci) {
    for (int bp = 0; bp < 32; ++bp) {
      for (int kind = 0; kind < 7; ++kind) {
        sao_info_t sao;
        memset(&sao, 0, sizeof(sao));
        sao.type = SAO_TYPE_BAND;
        sao.band_position[0] = bp;
        sao.band_position[1] = 31 - bp;
        for (int k = 0; k < NUM_SAO_EDGE_CATEGORIES * 2; ++k) {
          sao.offsets[k] = kind < 5 ? rnd_offset(kind) : (int)(rnd() % 20001) - 10000;
        }
        sao_offset_array_reference(&sao, ref, colors[ci]);
        memset(out, 0x5A, sizeof(out));
        kvz_calc_sao_offset_array(&enc, &sao, out, colors[ci]);
        if (memcmp(ref, out, sizeof(ref)) != 0) {
          printf("\n    kvz_calc_sao_offset_array mismatch color=%d band_pos=%d kind=%d\n", ci, bp, kind);
          FAILm("mismatch");
        }
      }
    }
  }
  PASS();
}

TEST test_sao_edge_ddistortion(void)
{
  ASSERT(run_edge_ddistortion());
  PASS();
}

TEST test_calc_sao_edge_dir(void)
{
  ASSERT(run_calc_sao_edge_dir());
  PASS();
}

TEST test_sao_band_ddistortion(void)
{
  ASSERT(run_band_ddistortion());
  PASS();
}

TEST test_sao_reconstruct_color(void)
{
  ASSERT(run_sao_reconstruct_color());
  PASS();
}

SUITE(simd_intra_sao_tests)
{
  memset(&enc, 0, sizeof(enc));
  enc.bitdepth = KVZ_BIT_DEPTH;
  memset(&st, 0, sizeof(st));
  st.encoder_control = &enc;

  RUN_TEST(test_sao_offset_array);

  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const strategy_t *s = &strategies.strategies[i];
    if (strcmp(s->strategy_name, "generic") == 0) continue;
    void *generic = find_generic(s->type);
    if (!generic) continue;
    cur = s;
    intra_ref_fn = sao_ref_fn = generic;
    intra_opt_fn = sao_opt_fn = s->fptr;
    if (strcmp(s->type, "angular_pred") == 0) {
      RUN_TEST(test_angular_pred);
    } else if (strcmp(s->type, "intra_pred_planar") == 0) {
      RUN_TEST(test_intra_pred_planar);
    } else if (strcmp(s->type, "intra_pred_filtered_dc") == 0) {
      RUN_TEST(test_intra_pred_filtered_dc);
    } else if (strcmp(s->type, "sao_edge_ddistortion") == 0) {
      RUN_TEST(test_sao_edge_ddistortion);
    } else if (strcmp(s->type, "calc_sao_edge_dir") == 0) {
      RUN_TEST(test_calc_sao_edge_dir);
    } else if (strcmp(s->type, "sao_band_ddistortion") == 0) {
      RUN_TEST(test_sao_band_ddistortion);
    } else if (strcmp(s->type, "sao_reconstruct_color") == 0) {
      RUN_TEST(test_sao_reconstruct_color);
    }
  }
}
