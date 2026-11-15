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
 * Differential tests for the quantization, residual and transform
 * strategies: every non-generic implementation of quant, dequant,
 * quantize_residual, find_last_scanpos, coeff_abs_sum, fast_coeff_cost and
 * the transforms is compared with the generic one at KVZ_BIT_DEPTH, on
 * random, extreme and edge-case inputs, over the whole QP range of the bit
 * depth, with and without scaling lists, sign hiding, RDOQ, transform skip
 * and cross-component prediction.
 */

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "context.h"
#include "encoder.h"
#include "encoderstate.h"
#include "rdo.h"
#include "scalinglist.h"
#include "tables.h"
#include "transform.h"
#include "strategies/strategies-dct.h"
#include "strategies/strategies-quant.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define QP_BD_OFFSET (6 * (KVZ_BIT_DEPTH - 8))

static uint32_t rng_state = 0x6A09E667u;
static uint32_t rnd(void)
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state >> 8;
}

static const strategy_t *cur;   // strategy under test
static void *ref_fn;            // generic implementation of the same type

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

/* ---------------- Encoder state ---------------- */

// Scaling lists: flat (off), the default lists, custom lists with entries of
// 13 to 255 and custom lists with entries of 1 to 255. With entries below 13
// the quantization scales (16 * 26214 / entry) do not fit in 16 bits, which
// the SIMD functions handle separately.
#define NUM_SL 4
static encoder_control_t enc[NUM_SL];
static encoder_state_t st;
static encoder_state_config_frame_t frame_cfg;
static encoder_state_config_tile_t tile_cfg;
static videoframe_t vframe;

static void init_encoders(void)
{
  for (int s = 0; s < NUM_SL; ++s) {
    encoder_control_t *e = &enc[s];
    memset(e, 0, sizeof(*e));
    e->bitdepth = KVZ_BIT_DEPTH;
    e->cfg.chroma_format = KVZ_CSP_420;
    kvz_scalinglist_init(&e->scaling_list);
    if (s == 1) {
      e->scaling_list.enable = 1;
      e->scaling_list.use_default_list = 1;
    } else if (s >= 2) {
      const int min_entry = s == 2 ? 13 : 1;
      e->scaling_list.enable = 1;
      for (int size = 0; size < SCALING_LIST_SIZE_NUM; ++size) {
        for (int list = 0; list < kvz_g_scaling_list_num[size]; ++list) {
          coeff_t *c = (coeff_t *)e->scaling_list.scaling_list_coeff[size][list];
          const int n = MIN(MAX_MATRIX_COEF_NUM, kvz_g_scaling_list_size[size]);
          for (int i = 0; i < n; ++i) {
            c[i] = (coeff_t)(min_entry + rnd() % (256 - min_entry));
          }
          e->scaling_list.scaling_list_dc[size][list] = min_entry + rnd() % (256 - min_entry);
        }
      }
    }
    kvz_scalinglist_process(&e->scaling_list, KVZ_BIT_DEPTH);
  }
  memset(&st, 0, sizeof(st));
  memset(&frame_cfg, 0, sizeof(frame_cfg));
  memset(&tile_cfg, 0, sizeof(tile_cfg));
  memset(&vframe, 0, sizeof(vframe));
  vframe.width = 64;
  vframe.height = 64;
  tile_cfg.frame = &vframe;
  st.frame = &frame_cfg;
  st.tile = &tile_cfg;
  st.encoder_control = &enc[0];
}

static void destroy_encoders(void)
{
  for (int s = 0; s < NUM_SL; ++s) {
    kvz_scalinglist_destroy(&enc[s].scaling_list);
  }
}

static void set_state(int sl, int qp, enum kvz_slice_type slicetype)
{
  st.encoder_control = &enc[sl];
  st.qp = (int8_t)qp;
  st.lambda = 0.57 * pow(2.0, (CLIP(0, 51, qp) - 12) / 3.0);
  frame_cfg.slicetype = slicetype;
  kvz_init_contexts(&st, (int8_t)CLIP(0, 51, qp), slicetype);
  // As kvz_search_lcu: the search counts bits with a copy of a started CABAC.
  kvz_cabac_start(&st.cabac);
  st.search_cabac = st.cabac;
  st.search_cabac.only_count = 1;
}

// The encoder clips the QP to [0, 51]; at high bit depth the quantizers add
// QP_BD_OFFSET to it.
static int random_qp(void)
{
  return (int)(rnd() % 52);
}

/* ---------------- Coefficient generators ---------------- */

static ALIGNED(64) coeff_t c_in[32 * 32];
static ALIGNED(64) coeff_t c_ref[32 * 32];
static ALIGNED(64) coeff_t c_opt[32 * 32];

// mode 0: full int16 range, 1: typical transform output, 2: small, 3: sparse,
// 4: extremes, 5: zero
#define NUM_COEFF_MODES 6
static coeff_t gen_coeff(int mode, int i, int n)
{
  switch (mode) {
    case 0: return (coeff_t)((int)(rnd() & 0xFFFF) - 32768);
    case 1: return (coeff_t)((int)(rnd() % 4001) - 2000) * ((i < n / 8) ? 8 : 1);
    case 2: return (coeff_t)((int)(rnd() % 9) - 4);
    case 3: return (rnd() % 23 == 0) ? (coeff_t)((int)(rnd() % 1001) - 500) : 0;
    case 4: return (rnd() & 1) ? 32767 : -32768;
    default: return 0;
  }
}

static void fill_coeffs(int mode, int n)
{
  for (int i = 0; i < n; ++i) c_in[i] = gen_coeff(mode, i, n);
}

#define FAIL_DIFF(fmt, ...) do { \
  printf("\n    %s:%s " fmt "\n", cur->type, cur->strategy_name, __VA_ARGS__); \
  FAILm("mismatch"); } while (0)

/* ---------------- quant ---------------- */

TEST test_quant(void)
{
  quant_func *ref = ref_fn, *opt = cur->fptr;
  for (int trial = 0; trial < 1500; ++trial) {
    const int sl = trial % NUM_SL;
    const int width = 4 << (rnd() % 4);
    const int n = width * width;
    const int mode = rnd() % NUM_COEFF_MODES;
    const int qp = (trial < 2 * 52) ? (trial / 2) : random_qp();
    const int8_t type = (rnd() & 1) ? 2 : 0;
    const int8_t scan_idx = rnd() % 3;
    const int8_t block_type = (rnd() & 1) ? CU_INTRA : CU_INTER;
    enc[sl].cfg.signhide_enable = rnd() & 1;
    enc[sl].cfg.chroma_format = (rnd() & 1) ? KVZ_CSP_420 : KVZ_CSP_444;
    set_state(sl, qp, (rnd() & 1) ? KVZ_SLICE_I : KVZ_SLICE_P);
    fill_coeffs(mode, n);
    memset(c_ref, 0x55, sizeof(c_ref));
    memset(c_opt, 0x55, sizeof(c_opt));
    ref(&st, c_in, c_ref, width, width, type, scan_idx, block_type);
    opt(&st, c_in, c_opt, width, width, type, scan_idx, block_type);
    for (int i = 0; i < 32 * 32; ++i) {
      if (c_ref[i] != c_opt[i]) {
        FAIL_DIFF("w=%d qp=%d sl=%d mode=%d type=%d scan=%d sh=%d idx=%d ref=%d opt=%d", width, qp, sl, mode,
                  type, scan_idx, enc[sl].cfg.signhide_enable, i, c_ref[i], c_opt[i]);
      }
    }
  }
  PASS();
}

/* ---------------- dequant ---------------- */

TEST test_dequant(void)
{
  dequant_func *ref = ref_fn, *opt = cur->fptr;
  for (int trial = 0; trial < 1500; ++trial) {
    const int sl = trial % NUM_SL;
    const int width = 4 << (rnd() % 4);
    const int n = width * width;
    const int mode = rnd() % NUM_COEFF_MODES;
    const int qp = (trial < 2 * 52) ? (trial / 2) : random_qp();
    static const int8_t types[3] = { 0, 2, 3 };
    const int8_t type = types[rnd() % 3];
    const int8_t block_type = (rnd() & 1) ? CU_INTRA : CU_INTER;
    enc[sl].cfg.chroma_format = (rnd() & 1) ? KVZ_CSP_420 : KVZ_CSP_422;
    set_state(sl, qp, KVZ_SLICE_P);
    // Quantized levels: the full range, or the small levels quant produces.
    // Without scaling lists dequant multiplies the level by
    // inv_quant_scale << (qp_scaled / 6) in 32 bits. At 12-bit and high QPs
    // full-range levels would overflow, which quant never produces (its levels
    // dequantize back to the int16 transform range), so the levels are limited
    // to what fits.
    const int max_scaled_qp = qp + (KVZ_BIT_DEPTH - 8) * 6;
    const int level_max = MIN(32767, 0x7FFF0000 / (72 << (max_scaled_qp / 6)));
    for (int i = 0; i < n; ++i) {
      const int level = mode < 2 ? gen_coeff(mode == 0 ? 0 : 4, i, n) : gen_coeff(mode, i, n);
      c_in[i] = (coeff_t)CLIP(-level_max, level_max, level);
    }
    memset(c_ref, 0x55, sizeof(c_ref));
    memset(c_opt, 0x55, sizeof(c_opt));
    ref(&st, c_in, c_ref, width, width, type, block_type);
    opt(&st, c_in, c_opt, width, width, type, block_type);
    for (int i = 0; i < 32 * 32; ++i) {
      if (c_ref[i] != c_opt[i]) {
        FAIL_DIFF("w=%d qp=%d sl=%d mode=%d type=%d idx=%d in=%d ref=%d opt=%d", width, qp, sl, mode,
                  type, i, i < n ? c_in[i] : 0, c_ref[i], c_opt[i]);
      }
    }
  }
  PASS();
}

/* ---------------- coeff_abs_sum and fast_coeff_cost ---------------- */

TEST test_coeff_abs_sum(void)
{
  coeff_abs_sum_func *ref = ref_fn, *opt = cur->fptr;
  static const size_t lens[] = { 8, 16, 32, 64, 128, 256, 512, 1024 };
  for (int trial = 0; trial < 200; ++trial) {
    const int mode = trial % NUM_COEFF_MODES;
    fill_coeffs(mode, 32 * 32);
    for (unsigned li = 0; li < sizeof(lens) / sizeof(lens[0]); ++li) {
      const uint32_t r = ref(c_in, lens[li]), o = opt(c_in, lens[li]);
      if (r != o) FAIL_DIFF("len=%u mode=%d ref=%u opt=%u", (unsigned)lens[li], mode, r, o);
    }
  }
  PASS();
}

TEST test_fast_coeff_cost(void)
{
  fast_coeff_cost_func *ref = ref_fn, *opt = cur->fptr;
  for (int trial = 0; trial < 200; ++trial) {
    const int mode = trial % NUM_COEFF_MODES;
    const uint64_t weights = ((uint64_t)(rnd() & 0xFFFF) << 48) | ((uint64_t)(rnd() & 0xFFFF) << 32) |
                             ((uint64_t)(rnd() & 0xFFFF) << 16) | (uint64_t)(rnd() & 0xFFFF);
    fill_coeffs(mode, 32 * 32);
    for (int width = 4; width <= 32; width *= 2) {
      const double r = ref(c_in, width, weights), o = opt(c_in, width, weights);
      if (r != o) FAIL_DIFF("w=%d mode=%d ref=%f opt=%f", width, mode, r, o);
    }
  }
  PASS();
}

/* ---------------- find_last_scanpos ---------------- */

static struct kvz_sh_rates_t sh_ref, sh_opt;

TEST test_find_last_scanpos(void)
{
  find_last_scanpos_func *ref = ref_fn, *opt = cur->fptr;
  for (int trial = 0; trial < 3000; ++trial) {
    const int sl = trial % NUM_SL;
    const int width = 4 << (rnd() % 4);
    const int n = width * width;
    const int log2_size = kvz_g_convert_to_bit[width] + 2;
    const int qp = random_qp();
    const int8_t type = (rnd() & 1) ? 2 : 0;
    const int8_t scan_mode = rnd() % 3;
    const int8_t block_type = (rnd() & 1) ? CU_INTRA : CU_INTER;
    const encoder_control_t *e = &enc[sl];

    // As in kvz_rdoq.
    const int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - KVZ_BIT_DEPTH - log2_size;
    const int32_t sl_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
    const int32_t qp_scaled = kvz_get_scaled_qp(type, (int8_t)qp, QP_BD_OFFSET, true);
    const int32_t q_bits = 14 + qp_scaled / 6 + transform_shift;
    const int32_t *quant_coeff = e->scaling_list.quant_coeff[log2_size - 2][sl_type][qp_scaled % 6];
    const uint32_t *scan = kvz_g_sig_last_scan[scan_mode][log2_size - 1];
    const uint32_t cg_num = n >> 4;

    // Sparse coefficients (the usual case), with the last non-zero one at a
    // random position, or one of the coefficient modes.
    const int mode = rnd() % (NUM_COEFF_MODES + 2);
    if (mode < NUM_COEFF_MODES) {
      fill_coeffs(mode, n);
    } else {
      memset(c_in, 0, sizeof(c_in));
      const int last = rnd() % n;
      for (int i = 0; i <= last; ++i) {
        if (i == last || rnd() % 4 == 0) {
          c_in[scan[i]] = (coeff_t)((int)(rnd() % (mode == NUM_COEFF_MODES ? 7 : 401)) - (mode == NUM_COEFF_MODES ? 3 : 200));
        }
      }
    }

    int32_t last_r = -1, last_o = -1, cg_last_r = -1, cg_last_o = -1, cg_pos_r = 1234, cg_pos_o = 1234;
    uint16_t ctx_r = 0x5555, ctx_o = 0x5555;
    for (int i = 0; i < n; ++i) {
      c_ref[i] = c_opt[i] = (coeff_t)(i * 7 + 1);
      sh_ref.sig_coeff_inc[i] = sh_opt.sig_coeff_inc[i] = -i - 1;
    }
    ref(c_in, c_ref, type, q_bits, quant_coeff, &sh_ref, 16, &ctx_r, scan, &cg_last_r, &last_r, cg_num, &cg_pos_r, width, scan_mode);
    opt(c_in, c_opt, type, q_bits, quant_coeff, &sh_opt, 16, &ctx_o, scan, &cg_last_o, &last_o, cg_num, &cg_pos_o, width, scan_mode);

    if (last_r != last_o || cg_last_r != cg_last_o || cg_pos_r != cg_pos_o || ctx_r != ctx_o) {
      FAIL_DIFF("w=%d qp=%d sl=%d mode=%d scan=%d last %d/%d cg_last %d/%d cg_pos %d/%d ctx %d/%d", width, qp, sl,
                mode, scan_mode, last_r, last_o, cg_last_r, cg_last_o, cg_pos_r, cg_pos_o, ctx_r, ctx_o);
    }
    // kvz_rdoq overwrites the levels at and before the last position, so only
    // the zeroed coefficients after it are part of the result.
    for (int i = last_r + 1; i < n; ++i) {
      if (c_ref[scan[i]] != c_opt[scan[i]]) {
        FAIL_DIFF("w=%d qp=%d mode=%d scanpos=%d last=%d dest ref=%d opt=%d", width, qp, mode, i, last_r,
                  c_ref[scan[i]], c_opt[scan[i]]);
      }
    }
    for (int i = 0; i < n; ++i) {
      if (sh_ref.sig_coeff_inc[i] != sh_opt.sig_coeff_inc[i]) {
        FAIL_DIFF("w=%d qp=%d mode=%d sig_coeff_inc[%d] ref=%d opt=%d", width, qp, mode, i,
                  sh_ref.sig_coeff_inc[i], sh_opt.sig_coeff_inc[i]);
      }
    }
  }
  PASS();
}

/* ---------------- quantize_residual ---------------- */

#define PIX_STRIDE 64
#define OUT_STRIDE 40
static ALIGNED(64) kvz_pixel ref_px[PIX_STRIDE * 32];
static ALIGNED(64) kvz_pixel pred_px[PIX_STRIDE * 32];
static ALIGNED(64) kvz_pixel rec_r[PIX_STRIDE * 32], rec_o[PIX_STRIDE * 32];
static ALIGNED(64) int16_t ccp_r[2][64 * 64], ccp_o[2][64 * 64];

// mode 0: random, 1: prediction close to the reference, 2: random 0/PIXEL_MAX,
// 3: maximum positive residual, 4: maximum negative residual, 5: flat
#define NUM_PIXEL_MODES 6
static void fill_pixels(int mode)
{
  const int noise = 1 + rnd() % 64;
  for (int i = 0; i < PIX_STRIDE * 32; ++i) {
    switch (mode) {
      case 0:
        ref_px[i] = (kvz_pixel)(rnd() & PIXEL_MAX);
        pred_px[i] = (kvz_pixel)(rnd() & PIXEL_MAX);
        break;
      case 1:
        pred_px[i] = (kvz_pixel)(rnd() & PIXEL_MAX);
        ref_px[i] = (kvz_pixel)CLIP(0, PIXEL_MAX, pred_px[i] + (int)(rnd() % (2 * noise + 1)) - noise);
        break;
      case 2:
        ref_px[i] = (rnd() & 1) ? PIXEL_MAX : 0;
        pred_px[i] = (rnd() & 1) ? PIXEL_MAX : 0;
        break;
      case 3:
        ref_px[i] = PIXEL_MAX;
        pred_px[i] = 0;
        break;
      case 4:
        ref_px[i] = 0;
        pred_px[i] = PIXEL_MAX;
        break;
      default:
        ref_px[i] = (kvz_pixel)(PIXEL_MAX / 3);
        pred_px[i] = (kvz_pixel)(PIXEL_MAX / 3);
        break;
    }
  }
}

static int check_residual_call(quant_residual_func *ref, quant_residual_func *opt, cu_info_t *cu_r, cu_info_t *cu_o,
  int width, color_t color, coeff_scan_order_t scan_order, int use_trskip, bool early_skip, bool in_place)
{
  int16_t *ccp_ptr_r[2] = { ccp_r[0], ccp_r[1] };
  int16_t *ccp_ptr_o[2] = { ccp_o[0], ccp_o[1] };
  const kvz_pixel *pred_r = pred_px, *pred_o = pred_px;
  int out_stride = OUT_STRIDE;
  if (in_place) {
    // The reconstruction is written over the prediction.
    memcpy(rec_r, pred_px, sizeof(pred_px));
    memcpy(rec_o, pred_px, sizeof(pred_px));
    pred_r = rec_r;
    pred_o = rec_o;
    out_stride = PIX_STRIDE;
  } else {
    memset(rec_r, 0x5a, sizeof(rec_r));
    memset(rec_o, 0x5a, sizeof(rec_o));
  }
  memset(c_ref, 0x55, sizeof(c_ref));
  memset(c_opt, 0x55, sizeof(c_opt));

  const int has_r = ref(&st, cu_r, width, color, scan_order, use_trskip, PIX_STRIDE, out_stride,
                        ref_px, pred_r, rec_r, c_ref, early_skip, ccp_ptr_r);
  const int has_o = opt(&st, cu_o, width, color, scan_order, use_trskip, PIX_STRIDE, out_stride,
                        ref_px, pred_o, rec_o, c_opt, early_skip, ccp_ptr_o);
  if (has_r != has_o) return 1;
  if (memcmp(c_ref, c_opt, sizeof(c_ref))) return 2;
  if (memcmp(rec_r, rec_o, sizeof(rec_r))) return 3;
  if (memcmp(cu_r, cu_o, sizeof(*cu_r))) return 4;
  if (memcmp(ccp_r, ccp_o, sizeof(ccp_r))) return 5;
  return 0;
}

TEST test_quantize_residual(void)
{
  quant_residual_func *ref = ref_fn, *opt = cur->fptr;
  for (int trial = 0; trial < 1200; ++trial) {
    const int sl = trial % NUM_SL;
    encoder_control_t *e = &enc[sl];
    const int width = 4 << (rnd() % 4);
    const int mode = rnd() % NUM_PIXEL_MODES;
    const int qp = (trial < 2 * 52) ? (trial / 2) : random_qp();
    const bool ccp = (trial % 4) == 3;
    e->cfg.rdoq_enable = rnd() & 1;
    e->cfg.rdoq_skip = rnd() & 1;
    e->cfg.signhide_enable = rnd() & 1;
    e->cfg.trskip_enable = 1;
    e->cfg.enable_cross_component_prediction = ccp;
    e->cfg.chroma_format = ccp ? KVZ_CSP_444 : (enum kvz_chroma_format)(1 + rnd() % 3);
    set_state(sl, qp, (rnd() & 1) ? KVZ_SLICE_I : KVZ_SLICE_P);
    fill_pixels(mode);

    cu_info_t cu;
    memset(&cu, 0, sizeof(cu));
    cu.type = (rnd() & 1) ? CU_INTRA : CU_INTER;
    cu.depth = 3 - kvz_g_convert_to_bit[width] + (width == 4 ? (rnd() & 1) : 0);
    cu.tr_depth = (ccp && (rnd() % 4)) ? cu.depth : cu.depth + (rnd() & 1);
    cu.part_size = (rnd() & 1) ? SIZE_NxN : SIZE_2Nx2N;
    cu.cbf = (uint16_t)rnd();
    cu.alpha_u = rnd() & 7;
    cu.alpha_v = rnd() & 7;
    cu.alpha_u_s = rnd() & 1;
    cu.alpha_v_s = rnd() & 1;
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 64 * 64; ++i) {
        ccp_r[k][i] = ccp_o[k][i] = (int16_t)((int)(rnd() % (2 * PIXEL_MAX + 1)) - PIXEL_MAX);
      }
    }
    cu_info_t cu_r = cu, cu_o = cu;
    // Horizontal and vertical scans are only used for 4x4 and 8x8 blocks.
    const coeff_scan_order_t scan_order = width <= 8 ? rnd() % 3 : SCAN_DIAG;
    const int use_trskip = width == 4 && (rnd() & 1);
    const bool early_skip = rnd() % 8 == 0;
    const bool in_place = rnd() % 4 == 0;

    // With cross-component prediction the chroma blocks use the luma residuals
    // stored by the luma call, so the planes are coded in order.
    for (int c = 0; c < 3; ++c) {
      const color_t color = (color_t)c;
      if (c > 0 && !ccp && (rnd() & 1)) continue;
      if (c > 0) fill_pixels(mode == 1 ? 1 : rnd() % NUM_PIXEL_MODES);
      const int err = check_residual_call(ref, opt, &cu_r, &cu_o, width, color, scan_order, use_trskip,
                                          early_skip, in_place);
      if (err) {
        FAIL_DIFF("trial=%d w=%d color=%d qp=%d sl=%d mode=%d rdoq=%d sh=%d trskip=%d early=%d in_place=%d "
                  "ccp=%d csp=%d err=%d", trial, width, c, qp, sl, mode, e->cfg.rdoq_enable,
                  e->cfg.signhide_enable, use_trskip, early_skip, in_place, ccp, e->cfg.chroma_format, err);
      }
    }
  }
  PASS();
}

/* ---------------- Transforms ---------------- */

static ALIGNED(64) int16_t tr_in[32 * 32];
static ALIGNED(64) int16_t tr_ref[32 * 32];
static ALIGNED(64) int16_t tr_opt[32 * 32];

// Forward transforms take residuals of KVZ_BIT_DEPTH bits, inverse transforms
// dequantized coefficients in the full int16 range.
TEST test_transform(void)
{
  dct_func *ref = ref_fn, *opt = cur->fptr;
  const bool inverse = strstr(cur->type, "idct") || strstr(cur->type, "inverse");
  int n = 4;
  if (strstr(cur->type, "8x8")) n = 8;
  if (strstr(cur->type, "16x16")) n = 16;
  if (strstr(cur->type, "32x32")) n = 32;
  for (int mode = 0; mode < 7; ++mode) {
    for (int trial = 0; trial < 64; ++trial) {
      for (int i = 0; i < n * n; ++i) {
        int v;
        if (!inverse) {
          v = mode == 0 ? (int)(rnd() % (2 * PIXEL_MAX + 1)) - PIXEL_MAX :
              mode == 1 ? ((rnd() & 1) ? PIXEL_MAX : -PIXEL_MAX) :
              mode == 2 ? PIXEL_MAX : mode == 3 ? -PIXEL_MAX :
              mode == 4 ? ((i ^ (i >> 3)) & 1 ? PIXEL_MAX : -PIXEL_MAX) :
              mode == 5 ? (int)(rnd() % 65) - 32 :
              ((rnd() >> 3) & 1) ? (((i % n) + (i / n) + trial) & 1 ? PIXEL_MAX : -PIXEL_MAX) : 0;
        } else {
          v = mode == 0 ? (int)(rnd() & 0xFFFF) - 32768 : mode == 1 ? ((rnd() & 1) ? 32767 : -32768) :
              mode == 2 ? 32767 : mode == 3 ? -32768 : mode == 4 ? (int)(rnd() % 2001) - 1000 :
              mode == 5 ? (rnd() % 7 == 0 ? (int)(rnd() & 0xFFFF) - 32768 : 0) :
              ((i % n) + (i / n)) & 1 ? 32767 : -32768;
        }
        tr_in[i] = (int16_t)v;
      }
      memset(tr_ref, 0x55, sizeof(tr_ref));
      memset(tr_opt, 0x55, sizeof(tr_opt));
      ref(KVZ_BIT_DEPTH, tr_in, tr_ref);
      opt(KVZ_BIT_DEPTH, tr_in, tr_opt);
      for (int i = 0; i < 32 * 32; ++i) {
        if (tr_ref[i] != tr_opt[i]) {
          FAIL_DIFF("mode=%d trial=%d idx=%d ref=%d opt=%d", mode, trial, i, tr_ref[i], tr_opt[i]);
        }
      }
    }
  }
  PASS();
}

SUITE(simd_quant_tests)
{
  // KVZ_TEST_SEED selects another random sequence.
  if (getenv("KVZ_TEST_SEED")) rng_state = (uint32_t)strtoul(getenv("KVZ_TEST_SEED"), NULL, 10);
  init_encoders();

  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    cur = &strategies.strategies[i];
    if (strcmp(cur->strategy_name, "generic") == 0) continue;
    ref_fn = find_generic(cur->type);
    if (!ref_fn) continue;

    if (strcmp(cur->type, "quant") == 0) {
      RUN_TEST(test_quant);
    } else if (strcmp(cur->type, "dequant") == 0) {
      RUN_TEST(test_dequant);
    } else if (strcmp(cur->type, "coeff_abs_sum") == 0) {
      RUN_TEST(test_coeff_abs_sum);
    } else if (strcmp(cur->type, "fast_coeff_cost") == 0) {
      RUN_TEST(test_fast_coeff_cost);
    } else if (strcmp(cur->type, "find_last_scanpos") == 0) {
      RUN_TEST(test_find_last_scanpos);
    } else if (strcmp(cur->type, "quantize_residual") == 0) {
      RUN_TEST(test_quantize_residual);
    } else if (!strncmp(cur->type, "dct_", 4) || !strncmp(cur->type, "idct_", 5) ||
               !strcmp(cur->type, "fast_forward_dst_4x4") || !strcmp(cur->type, "fast_inverse_dst_4x4")) {
      RUN_TEST(test_transform);
    }
  }

  destroy_encoders();
}
