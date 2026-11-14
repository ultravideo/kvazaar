/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2026, project contributors
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
 * \ingroup Reconstruction
 * \file
 * Per-LCU chroma reconstruction from the final CU tree and coefficients.
 *
 * After the search of an LCU, kvazaar's frame->rec chroma can differ from
 * what the decoder reconstructs for some 4:2:2 blocks (the search's recon can
 * read a stale work-tree reference or leave sub-TU areas without a residual).
 * This pass re-runs ONLY the chroma reconstruction in decode order, reading
 * references from the progressively-rebuilt frame->rec and using the final
 * quantized coefficients, so frame->rec chroma matches the decoder. Luma is
 * left untouched (it is already correct).
 */

#include "global.h"

#include "cu.h"
#include "encoderstate.h"
#include "inter.h"
#include "intra.h"
#include "strategies/strategies-picture.h"
#include "transform.h"
#include "videoframe.h"

/**
 * Reconstruct only the chroma of a leaf CU into lcu->rec using the stored
 * coefficients and reading references from lcu->rec (which mirrors the
 * progressively-rebuilt frame->rec).
 */
static void reconstruct_cu_chroma(encoder_state_t * const state,
                                  int x, int y, int depth,
                                  lcu_t *lcu,
                                  const cu_info_t *cur_cu)
{
  const int cu_width = LCU_WIDTH >> depth;
  const bool has_chroma = state->encoder_control->cfg.chroma_format != KVZ_CSP_400;

  // The dequant/inverse-transform in the recon path uses state->qp, which
  // after the per-LCU search loop is stale. Use the block's own QP so the
  // residual reconstruction matches the decoder.
  state->qp = cur_cu->qp;

  if (cur_cu->type == CU_INTRA) {
    // Intra transform coding always runs, so the per-sub-TU CBFs are always
    // signalled; no overall-CBF gate here.
    kvz_intra_recon_cu(state, x, y, depth,
                       -1, cur_cu->intra.mode_chroma, // skip luma
                       NULL, lcu, true, false);
  } else if (cur_cu->type == CU_INTER) {
    kvz_inter_recon_cu(state, lcu, x, y, cu_width, false, has_chroma);
    // For inter CUs the transform coding (and hence all sub-TU chroma CBFs)
    // is only invoked when the CU's overall CBF is set. The cu_array can carry
    // a stale per-sub-position CBF from a rejected candidate, so gate the whole
    // residual on the CU's own overall CBF to match the decoder.
    if (cbf_is_set_any(cur_cu->cbf, depth)) {
      kvz_quantize_lcu_residual(state, false, has_chroma, x, y, depth,
                                NULL, lcu, false, KVZ_SUBTU_ALL, true);
    }
  }
}

/**
 * Walk the CU tree of an LCU in decode order and reconstruct every leaf's
 * chroma.
 */
static void reconstruct_coding_tree_chroma(encoder_state_t * const state,
                                           int x, int y, int depth,
                                           lcu_t *lcu)
{
  const videoframe_t * const frame = state->tile->frame;
  const cu_info_t *cur_cu = kvz_cu_array_at_const(frame->cu_array, x, y);
  const int cu_width = LCU_WIDTH >> depth;

  if (GET_SPLITDATA(cur_cu, depth)) {
    const int half_cu = cu_width >> 1;
    reconstruct_coding_tree_chroma(state, x,          y,          depth + 1, lcu);
    reconstruct_coding_tree_chroma(state, x + half_cu, y,          depth + 1, lcu);
    reconstruct_coding_tree_chroma(state, x,          y + half_cu, depth + 1, lcu);
    reconstruct_coding_tree_chroma(state, x + half_cu, y + half_cu, depth + 1, lcu);
    return;
  }

  reconstruct_cu_chroma(state, x, y, depth, lcu, cur_cu);
}

/**
 * Copy lcu->rec chroma back to frame->rec.
 */
static void copy_lcu_chroma_to_frame_rec(const encoder_state_t * const state,
                                         const int x, const int y, const lcu_t *lcu)
{
  const videoframe_t * const frame = state->tile->frame;
  const int shift_w = SHIFT_W;
  const int shift_h = SHIFT_H;
  const int lcu_w_c = LCU_WIDTH >> shift_w;
  const int pic_width = frame->width;
  const int x_max = MIN(x + LCU_WIDTH, pic_width) - x;
  const int y_max = MIN(y + LCU_WIDTH, frame->height) - y;

  if (state->encoder_control->cfg.chroma_format != KVZ_CSP_400) {
    kvz_pixels_blit(lcu->rec.u, &frame->rec->u[(x >> shift_w) + (y >> shift_h) * (frame->rec->stride_c)],
                    x_max >> shift_w, y_max >> shift_h, lcu_w_c, frame->rec->stride_c);
    kvz_pixels_blit(lcu->rec.v, &frame->rec->v[(x >> shift_w) + (y >> shift_h) * (frame->rec->stride_c)],
                    x_max >> shift_w, y_max >> shift_h, lcu_w_c, frame->rec->stride_c);
  }
}

/**
 * Rebuild one LCU's frame->rec chroma from the final CU tree and the final
 * quantized coefficients so that it matches the decoder. Luma is left
 * untouched.
 *
 * This is called at the end of the LCU's search (kvz_search_lcu), before the
 * worker's normal per-LCU deblock/SAO, so the filters operate on the
 * corrected chroma and no frame-level pass is needed. The search's own
 * per-candidate recon can leave stale pixels for some 4:2:2 blocks (rejected
 * candidates, sub-TU areas with no residual), so the chroma is re-built in
 * decode order.
 *
 * The references are taken from the search's work tree (search_lcu), whose
 * top_ref/left_ref were captured from the reference buffers (recdata_to_bufs)
 * before the worker's per-LCU deblocking - i.e. the pre-deblock pixels the
 * decoder uses for its intra prediction. They must not be re-copied from
 * frame->rec, which by this point contains the neighbours' post-deblock
 * pixels. The lcu->rec starts empty and is filled leaf by leaf in decode
 * order, so each block's intra references read either the LCU border refs or
 * the already-rebuilt reconstruction of its top/left neighbours.
 *
 * The coefficients must be COPIED: the recon-from-coeffs quantize path
 * dequantizes its coefficient buffer in place, and the search's final
 * coefficients are still needed for the bitstream coding (state->coeff).
 *
 * \param search_lcu   depth 0 of the search work tree
 * \param scratch_lcu  another lcu_t of the same work tree, which is
 *                     overwritten (except for the shared source pixels,
 *                     which are not used)
 */
void kvz_reconstruct_lcu_chroma(encoder_state_t * const state,
                                const int lcu_px_x, const int lcu_px_y,
                                const lcu_t *search_lcu, lcu_t *scratch_lcu)
{
  const videoframe_t * const frame = state->tile->frame;

  const enum kvz_chroma_format csp = state->encoder_control->cfg.chroma_format;
  const int chroma_size = KVZ_LCU_CHROMA_SIZE(csp);
  if (chroma_size == 0) return;

  // Only the chroma is reconstructed, so only the chroma buffers need to be
  // initialized.
  lcu_t *const lcu = scratch_lcu;
  FILL(lcu->cu, 0);
  FILL_ARRAY(lcu->rec.u, 0, KVZ_LCU_CHROMA_REC_SIZE(csp));
  FILL_ARRAY(lcu->rec.v, 0, KVZ_LCU_CHROMA_REC_SIZE(csp));
  memcpy(lcu->coeff.u, search_lcu->coeff.u, chroma_size * sizeof(coeff_t));
  memcpy(lcu->coeff.v, search_lcu->coeff.v, chroma_size * sizeof(coeff_t));
  lcu->top_ref  = search_lcu->top_ref;
  lcu->left_ref = search_lcu->left_ref;

  // Load the final CU tree for this LCU.
  kvz_cu_array_copy_to_lcu(lcu, lcu_px_x, lcu_px_y, frame->cu_array);

  reconstruct_coding_tree_chroma(state, lcu_px_x, lcu_px_y, 0, lcu);
  copy_lcu_chroma_to_frame_rec(state, lcu_px_x, lcu_px_y, lcu);
}
