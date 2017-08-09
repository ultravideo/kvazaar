/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "encoderstate.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cabac.h"
#include "context.h"
#include "encode_coding_tree.h"
#include "encoder_state-bitstream.h"
#include "filter.h"
#include "image.h"
#include "rate_control.h"
#include "sao.h"
#include "search.h"
#include "tables.h"
#include "threadqueue.h"


int kvz_encoder_state_match_children_of_previous_frame(encoder_state_t * const state) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    //Child should also exist for previous encoder
    assert(state->previous_encoder_state->children[i].encoder_control);
    state->children[i].previous_encoder_state = &state->previous_encoder_state->children[i];
    kvz_encoder_state_match_children_of_previous_frame(&state->children[i]);
  }
  return 1;
}

static void encoder_state_recdata_to_bufs(encoder_state_t * const state, const lcu_order_element_t * const lcu, yuv_t * const hor_buf, yuv_t * const ver_buf) {
  videoframe_t* const frame = state->tile->frame;
  
  if (hor_buf) {
    //Copy the bottom row of this LCU to the horizontal buffer
    vector2d_t bottom = { lcu->position_px.x, lcu->position_px.y + lcu->size.y - 1 };
    const int lcu_row = lcu->position.y;

    unsigned from_index = bottom.y * frame->rec->stride + bottom.x;
    unsigned to_index = lcu->position_px.x + lcu_row * frame->width;
    
    kvz_pixels_blit(&frame->rec->y[from_index],
                    &hor_buf->y[to_index],
                    lcu->size.x, 1,
                    frame->rec->stride, frame->width);

    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      unsigned from_index_c = (bottom.y / 2) * frame->rec->stride / 2 + (bottom.x / 2);
      unsigned to_index_c = lcu->position_px.x / 2 + lcu_row * frame->width / 2;

      kvz_pixels_blit(&frame->rec->u[from_index_c],
                      &hor_buf->u[to_index_c],
                      lcu->size.x / 2, 1, 
                      frame->rec->stride / 2, frame->width / 2);
      kvz_pixels_blit(&frame->rec->v[from_index_c],
                      &hor_buf->v[to_index_c],
                      lcu->size.x / 2, 1,
                      frame->rec->stride / 2, frame->width / 2);
    }
  }
  
  if (ver_buf) {
    //Copy the right row of this LCU to the vertical buffer.
    
    const int lcu_col = lcu->position.x;
    vector2d_t left = { lcu->position_px.x + lcu->size.x - 1, lcu->position_px.y };
    
    kvz_pixels_blit(&frame->rec->y[left.y * frame->rec->stride + left.x],
                    &ver_buf->y[lcu->position_px.y + lcu_col * frame->height],
                    1, lcu->size.y,
                    frame->rec->stride, 1);

    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      unsigned from_index = (left.y / 2) * frame->rec->stride / 2 + (left.x / 2);
      unsigned to_index = lcu->position_px.y / 2 + lcu_col * frame->height / 2;

      kvz_pixels_blit(&frame->rec->u[from_index],
                      &ver_buf->u[to_index],
                      1, lcu->size.y / 2,
                      frame->rec->stride / 2, 1);
      kvz_pixels_blit(&frame->rec->v[from_index],
                      &ver_buf->v[to_index],
                      1, lcu->size.y / 2,
                      frame->rec->stride / 2, 1);
    }
  }
  
}


static void encode_sao_color(encoder_state_t * const state, sao_info_t *sao,
                             color_t color_i)
{
  cabac_data_t * const cabac = &state->cabac;
  sao_eo_cat i;
  int offset_index = (color_i == COLOR_V) ? 5 : 0;

  // Skip colors with no SAO.
  //FIXME: for now, we always have SAO for all channels
  if (color_i == COLOR_Y && 0) return;
  if (color_i != COLOR_Y && 0) return;

  /// sao_type_idx_luma:   TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  /// sao_type_idx_chroma: TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  // Encode sao_type_idx for Y and U+V.
  if (color_i != COLOR_V) {
    cabac->cur_ctx = &(cabac->ctx.sao_type_idx_model);
    CABAC_BIN(cabac, sao->type != SAO_TYPE_NONE, "sao_type_idx");
    if (sao->type == SAO_TYPE_BAND) {
      CABAC_BIN_EP(cabac, 0, "sao_type_idx_ep");
    } else if (sao->type == SAO_TYPE_EDGE) {
      CABAC_BIN_EP(cabac, 1, "sao_type_idx_ep");
    }
  }

  if (sao->type == SAO_TYPE_NONE) return;

  /// sao_offset_abs[][][][]: TR, cMax = (1 << (Min(bitDepth, 10) - 5)) - 1,
  ///                         cRiceParam = 0, bins = {bypass x N}
  for (i = SAO_EO_CAT1; i <= SAO_EO_CAT4; ++i) {
    kvz_cabac_write_unary_max_symbol_ep(cabac, abs(sao->offsets[i + offset_index]), SAO_ABS_OFFSET_MAX);
  }

  /// sao_offset_sign[][][][]: FL, cMax = 1, bins = {bypass}
  /// sao_band_position[][][]: FL, cMax = 31, bins = {bypass x N}
  /// sao_eo_class_luma:       FL, cMax = 3, bins = {bypass x 3}
  /// sao_eo_class_chroma:     FL, cMax = 3, bins = {bypass x 3}
  if (sao->type == SAO_TYPE_BAND) {
    for (i = SAO_EO_CAT1; i <= SAO_EO_CAT4; ++i) {
      // Positive sign is coded as 0.
      if (sao->offsets[i + offset_index] != 0) {
        CABAC_BIN_EP(cabac, sao->offsets[i + offset_index] < 0 ? 1 : 0, "sao_offset_sign");
      }
    }
    // TODO: sao_band_position
    // FL cMax=31 (5 bits)
    CABAC_BINS_EP(cabac, sao->band_position[color_i == COLOR_V ? 1:0], 5, "sao_band_position");
  } else if (color_i != COLOR_V) {
    CABAC_BINS_EP(cabac, sao->eo_class, 2, "sao_eo_class");
  }
}

static void encode_sao_merge_flags(encoder_state_t * const state, sao_info_t *sao, unsigned x_ctb, unsigned y_ctb)
{
  cabac_data_t * const cabac = &state->cabac;
  // SAO merge flags are not present for the first row and column.
  if (x_ctb > 0) {
    cabac->cur_ctx = &(cabac->ctx.sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_left_flag, "sao_merge_left_flag");
  }
  if (y_ctb > 0 && !sao->merge_left_flag) {
    cabac->cur_ctx = &(cabac->ctx.sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_up_flag, "sao_merge_up_flag");
  }
}


/**
 * \brief Encode SAO information.
 */
static void encode_sao(encoder_state_t * const state,
                       unsigned x_lcu, uint16_t y_lcu,
                       sao_info_t *sao_luma, sao_info_t *sao_chroma)
{
  // TODO: transmit merge flags outside sao_info
  encode_sao_merge_flags(state, sao_luma, x_lcu, y_lcu);

  // If SAO is merged, nothing else needs to be coded.
  if (!sao_luma->merge_left_flag && !sao_luma->merge_up_flag) {
    encode_sao_color(state, sao_luma, COLOR_Y);
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      encode_sao_color(state, sao_chroma, COLOR_U);
      encode_sao_color(state, sao_chroma, COLOR_V);
    }
  }
}


/**
 * \brief Sets the QP for each CU in state->tile->frame->cu_array.
 *
 * The QPs are used in deblocking.
 *
 * The delta QP for an LCU is coded when the first CU with coded block flag
 * set is encountered. Hence, for the purposes of deblocking, all CUs
 * before the first one with cbf set use state->ref_qp and all CUs after
 * that use state->qp.
 *
 * \param state           encoder state
 * \param x               x-coordinate of the left edge of the root CU
 * \param y               y-coordinate of the top edge of the root CU
 * \param depth           depth in the CU quadtree
 * \param coeffs_coded    Used for tracking whether a CU with a residual
 *                        has been encountered. Should be set to false at
 *                        the top level.
 * \return Whether there were any CUs with residual or not.
 */
static bool set_cu_qps(encoder_state_t *state, int x, int y, int depth, bool coeffs_coded)
{
  if (state->qp == state->ref_qp) {
    // If the QPs are equal there is no need to care about the residuals.
    coeffs_coded = true;
  }

  cu_info_t *cu = kvz_cu_array_at(state->tile->frame->cu_array, x, y);
  const int cu_width = LCU_WIDTH >> depth;
  coeffs_coded = coeffs_coded || cbf_is_set_any(cu->cbf, cu->depth);

  if (!coeffs_coded && cu->depth > depth) {
    // Recursively process sub-CUs.
    const int d = cu_width >> 1;
    coeffs_coded = set_cu_qps(state, x,     y,     depth + 1, coeffs_coded);
    coeffs_coded = set_cu_qps(state, x + d, y,     depth + 1, coeffs_coded);
    coeffs_coded = set_cu_qps(state, x,     y + d, depth + 1, coeffs_coded);
    coeffs_coded = set_cu_qps(state, x + d, y + d, depth + 1, coeffs_coded);

  } else {
    if (!coeffs_coded && cu->tr_depth > depth) {
      // The CU is split into smaller transform units. Check whether coded
      // block flag is set for any of the TUs.
      const int tu_width = LCU_WIDTH >> cu->tr_depth;
      for (int y_scu = y; y_scu < y + cu_width; y_scu += tu_width) {
        for (int x_scu = x; x_scu < x + cu_width; x_scu += tu_width) {
          cu_info_t *tu = kvz_cu_array_at(state->tile->frame->cu_array, x_scu, y_scu);
          if (cbf_is_set_any(tu->cbf, cu->depth)) {
            coeffs_coded = true;
          }
        }
      }
    }

    // Set the correct QP for all state->tile->frame->cu_array elements in
    // the area covered by the CU.
    const int8_t qp = coeffs_coded ? state->qp : state->ref_qp;

    for (int y_scu = y; y_scu < y + cu_width; y_scu += SCU_WIDTH) {
      for (int x_scu = x; x_scu < x + cu_width; x_scu += SCU_WIDTH) {
        kvz_cu_array_at(state->tile->frame->cu_array, x_scu, y_scu)->qp = qp;
      }
    }
  }

  return coeffs_coded;
}


static void encoder_state_worker_encode_lcu(void * opaque)
{
  const lcu_order_element_t * const lcu = opaque;
  encoder_state_t *state = lcu->encoder_state;
  const encoder_control_t * const encoder = state->encoder_control;
  videoframe_t* const frame = state->tile->frame;

  kvz_set_lcu_lambda_and_qp(state, lcu->position);

  lcu_coeff_t coeff;
  state->coeff = &coeff;

  //This part doesn't write to bitstream, it's only search, deblock and sao
  kvz_search_lcu(state, lcu->position_px.x, lcu->position_px.y, state->tile->hor_buf_search, state->tile->ver_buf_search);
    
  encoder_state_recdata_to_bufs(state, lcu, state->tile->hor_buf_search, state->tile->ver_buf_search);

  if (encoder->cfg.deblock_enable) {
    if (encoder->lcu_dqp_enabled) {
      set_cu_qps(state, lcu->position_px.x, lcu->position_px.y, 0, false);
    }

    kvz_filter_deblock_lcu(state, lcu->position_px.x, lcu->position_px.y);
  }

  if (encoder->cfg.sao_enable) {
    kvz_sao_search_lcu(state, lcu->position.x, lcu->position.y);
  }

  // Copy LCU cu_array to main states cu_array, because that is the only one
  // which is given to the next frame through image_list_t.
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);

    encoder_state_t *main_state = state;
    while (main_state->parent) main_state = main_state->parent;
    assert(main_state != state);

    const unsigned tile_x_px = state->tile->lcu_offset_x << LOG2_LCU_WIDTH;
    const unsigned tile_y_px = state->tile->lcu_offset_y << LOG2_LCU_WIDTH;
    const unsigned x_px = lcu->position_px.x;
    const unsigned y_px = lcu->position_px.y;
    kvz_cu_array_copy(main_state->tile->frame->cu_array,
                      x_px + tile_x_px, y_px + tile_y_px,
                      state->tile->frame->cu_array,
                      x_px, y_px,
                      LCU_WIDTH, LCU_WIDTH);

    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=copy_cuinfo,frame=%d,tile=%d", state->frame->num, state->tile->id);
  }
  
  //Now write data to bitstream (required to have a correct CABAC state)
  const uint64_t existing_bits = kvz_bitstream_tell(&state->stream);
  
  //Encode SAO
  if (encoder->cfg.sao_enable) {
    encode_sao(state, lcu->position.x, lcu->position.y, &frame->sao_luma[lcu->position.y * frame->width_in_lcu + lcu->position.x], &frame->sao_chroma[lcu->position.y * frame->width_in_lcu + lcu->position.x]);
  }
  

  // QP delta is not used when rate control is turned off.
  state->must_code_qp_delta = encoder->lcu_dqp_enabled;

  //Encode coding tree
  kvz_encode_coding_tree(state, lcu->position.x * LCU_WIDTH, lcu->position.y * LCU_WIDTH, 0);

  // Coeffs are not needed anymore.
  state->coeff = NULL;

  bool end_of_slice_segment_flag;
  if (state->encoder_control->cfg.slices & KVZ_SLICES_WPP) {
    // Slice segments end after each WPP row.
    end_of_slice_segment_flag = lcu->last_column;
  } else if (state->encoder_control->cfg.slices & KVZ_SLICES_TILES) {
    // Slices end after each tile.
    end_of_slice_segment_flag = lcu->last_column && lcu->last_row;
  } else {
    // Slice ends after the last row of the last tile.
    int last_tile_id = -1 + encoder->cfg.tiles_width_count * encoder->cfg.tiles_height_count;
    bool is_last_tile = state->tile->id == last_tile_id;
    end_of_slice_segment_flag = is_last_tile && lcu->last_column && lcu->last_row;
  }
  kvz_cabac_encode_bin_trm(&state->cabac, end_of_slice_segment_flag);

  {
    const bool end_of_tile = lcu->last_column && lcu->last_row;
    const bool end_of_wpp_row = encoder->cfg.wpp && lcu->last_column;


    if (end_of_tile || end_of_wpp_row) {
      if (!end_of_slice_segment_flag) {
        // end_of_sub_stream_one_bit
        kvz_cabac_encode_bin_trm(&state->cabac, 1);
      }

      // Finish the substream by writing out remaining state.
      kvz_cabac_finish(&state->cabac);

      // Write a rbsp_trailing_bits or a byte_alignment. The first one is used
      // for ending a slice_segment_layer_rbsp and the second one for ending
      // a substream. They are identical and align the byte stream.
      kvz_bitstream_put(state->cabac.stream, 1, 1);
      kvz_bitstream_align_zero(state->cabac.stream);

      kvz_cabac_start(&state->cabac);

      kvz_crypto_delete(&state->crypto_hdl);
    }
  }

  const uint32_t bits = kvz_bitstream_tell(&state->stream) - existing_bits;
  kvz_get_lcu_stats(state, lcu->position.x, lcu->position.y)->bits = bits;

  //Wavefronts need the context to be copied to the next row
  if (state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW && lcu->index == 1) {
    int j;
    //Find next encoder (next row)
    for (j=0; state->parent->children[j].encoder_control; ++j) {
      if (state->parent->children[j].wfrow->lcu_offset_y == state->wfrow->lcu_offset_y + 1) {
        //And copy context
        kvz_context_copy(&state->parent->children[j], state);
      }
    }
  }
  
  if (encoder->cfg.sao_enable && lcu->above) {
    // Add the post-deblocking but pre-SAO pixels of the LCU row above this
    // row to a buffer so this row can use them on it's own SAO
    // reconstruction.

    // The pixels need to be taken to from the LCU to the top-left, because
    // not all of the pixels could be deblocked before prediction of this
    // LCU was reconstructed.
    if (lcu->above->left) {
      encoder_state_recdata_to_bufs(state, lcu->above->left, state->tile->hor_buf_before_sao, NULL);
    }
    // If this is the last LCU in the row, we can save the pixels from the top
    // also, as they have been fully deblocked.
    if (!lcu->right) {
      encoder_state_recdata_to_bufs(state, lcu->above, state->tile->hor_buf_before_sao, NULL);
    }
  }
}

static void encoder_state_encode_leaf(encoder_state_t * const state)
{
  assert(state->is_leaf);
  assert(state->lcu_order_count > 0);

  const kvz_config *cfg = &state->encoder_control->cfg;

  state->ref_qp = state->frame->QP;

  if (cfg->crypto_features) {
    state->crypto_hdl = kvz_crypto_create();
    state->crypto_prev_pos = 0;
  }

  // Select whether to encode the frame/tile in current thread or to define
  // wavefront jobs for other threads to handle.
  bool wavefront = state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW;
  bool use_parallel_encoding = (wavefront && state->parent->children[1].encoder_control);
  if (!use_parallel_encoding) {
    // Encode every LCU in order and perform SAO reconstruction after every
    // frame is encoded. Deblocking and SAO search is done during LCU encoding.

    for (int i = 0; i < state->lcu_order_count; ++i) {
      PERFORMANCE_MEASURE_START(KVZ_PERF_LCU);

      encoder_state_worker_encode_lcu(&state->lcu_order[i]);

#ifdef KVZ_DEBUG
      {
        const lcu_order_element_t * const lcu = &state->lcu_order[i];
        PERFORMANCE_MEASURE_END(KVZ_PERF_LCU, state->encoder_control->threadqueue, "type=encode_lcu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", state->frame->num, state->tile->id, state->slice->id, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH + lcu->size.x - 1, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH + lcu->size.y - 1);
      }
#endif //KVZ_DEBUG
    }

    if (state->encoder_control->cfg.sao_enable) {
      PERFORMANCE_MEASURE_START(KVZ_PERF_SAOREC);
      kvz_sao_reconstruct_frame(state);
      PERFORMANCE_MEASURE_END(KVZ_PERF_SAOREC, state->encoder_control->threadqueue, "type=kvz_sao_reconstruct_frame,frame=%d,tile=%d,slice=%d,row=%d-%d,px_x=%d-%d,px_y=%d-%d", state->frame->num, state->tile->id, state->slice->id, state->lcu_order[0].position.y + state->tile->lcu_offset_y, state->lcu_order[state->lcu_order_count - 1].position.y + state->tile->lcu_offset_y,
        state->tile->lcu_offset_x * LCU_WIDTH, state->tile->frame->width + state->tile->lcu_offset_x * LCU_WIDTH - 1,
        state->tile->lcu_offset_y * LCU_WIDTH, state->tile->frame->height + state->tile->lcu_offset_y * LCU_WIDTH - 1
      );
    }
  } else {
    // Add each LCU in the wavefront row as it's own job to the queue.

    // Select which frame dependancies should be set to.
    const encoder_state_t * ref_state = NULL;
    if (cfg->gop_lowdelay &&
        cfg->gop_len > 0 &&
        state->previous_encoder_state != state)
    {
      // For LP-gop, depend on the state of the first reference.
      int ref_neg = cfg->gop[(state->frame->poc - 1) % cfg->gop_len].ref_neg[0];
      if (ref_neg > state->encoder_control->cfg.owf) {
        // If frame is not within OWF range, it's already done.
        ref_state = NULL;
      } else {
        ref_state = state->previous_encoder_state;
        while (ref_neg > 1) {
          ref_neg -= 1;
          ref_state = ref_state->previous_encoder_state;
        }
      }
    } else {
      // Otherwise, depend on the previous frame.
      ref_state = state->previous_encoder_state;
    }

    for (int i = 0; i < state->lcu_order_count; ++i) {
      const lcu_order_element_t * const lcu = &state->lcu_order[i];

#ifdef KVZ_DEBUG
      char job_description[256];
      sprintf(job_description, "type=encode_lcu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", state->frame->num, state->tile->id, state->slice->id, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH + lcu->size.x - 1, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH + lcu->size.y - 1);
#else
      char* job_description = NULL;
#endif
      kvz_threadqueue_free_job(&state->tile->wf_jobs[lcu->id]);
      state->tile->wf_jobs[lcu->id] = kvz_threadqueue_submit(state->encoder_control->threadqueue, encoder_state_worker_encode_lcu, (void*)lcu, 1, job_description);
      
      // If job object was returned, add dependancies and allow it to run.
      if (state->tile->wf_jobs[lcu->id]) {
        // Add inter frame dependancies when ecoding more than one frame at
        // once. The added dependancy is for the first LCU of each wavefront
        // row to depend on the reconstruction status of the row below in the
        // previous frame.
        if (ref_state != NULL &&
            state->previous_encoder_state->tqj_recon_done &&
            state->frame->slicetype != KVZ_SLICE_I)
        {
          if (!lcu->left) {
            const lcu_order_element_t * const ref_lcu = &ref_state->lcu_order[i];
            if (lcu->below) {
              kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], ref_lcu->below->encoder_state->tqj_recon_done);
            } else {
              kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], ref_lcu->encoder_state->tqj_recon_done);
            }
          }
        }

        // Add local WPP dependancy to the LCU on the left.
        if (lcu->left) {
          kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], state->tile->wf_jobs[lcu->id - 1]);
        }
        // Add local WPP dependancy to the LCU on the top right.
        if (lcu->above) {
          if (lcu->above->right) {
            kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], state->tile->wf_jobs[lcu->id - state->tile->frame->width_in_lcu + 1]);
          } else {
            kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], state->tile->wf_jobs[lcu->id - state->tile->frame->width_in_lcu]);
          }
        }

        kvz_threadqueue_job_unwait_job(state->encoder_control->threadqueue, state->tile->wf_jobs[lcu->id]);

        // In the case where SAO is not enabled, the wavefront row is
        // done when the last LCU in the row is done.
        if (!state->encoder_control->cfg.sao_enable && i + 1 == state->lcu_order_count) {
          assert(!state->tqj_recon_done);
          state->tqj_recon_done =
            kvz_threadqueue_copy_ref(state->tile->wf_jobs[lcu->id]);
        }
      }
    }
  }
}

static void encoder_state_encode(encoder_state_t * const main_state);

static void encoder_state_worker_encode_children(void * opaque)
{
  encoder_state_t *sub_state = opaque;
  encoder_state_encode(sub_state);

  if (sub_state->is_leaf && sub_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
    // Set the last wavefront job of this row as the job that completes
    // the bitstream for this wavefront row state.

    int wpp_row = sub_state->wfrow->lcu_offset_y;
    int tile_width = sub_state->tile->frame->width_in_lcu;
    int end_of_row = (wpp_row + 1) * tile_width - 1;
    assert(!sub_state->tqj_bitstream_written);
    if (sub_state->tile->wf_jobs[end_of_row]) {
      sub_state->tqj_bitstream_written =
        kvz_threadqueue_copy_ref(sub_state->tile->wf_jobs[end_of_row]);
    }
  }
}

typedef struct {
  int y;
  const encoder_state_t * encoder_state;
} worker_sao_reconstruct_lcu_data;

static void encoder_state_worker_sao_reconstruct_lcu(void *opaque) {
  worker_sao_reconstruct_lcu_data *data = opaque;
  videoframe_t * const frame = data->encoder_state->tile->frame;
  unsigned stride = frame->width_in_lcu;
  int x;
  
  //TODO: copy only needed data
  kvz_pixel *new_y_data = MALLOC(kvz_pixel, frame->width * frame->height);
  kvz_pixel *new_u_data = NULL;
  kvz_pixel *new_v_data = NULL;
  if (frame->rec->chroma_format != KVZ_CSP_400) {
    new_u_data = MALLOC(kvz_pixel, (frame->width * frame->height) >> 2);
    new_v_data = MALLOC(kvz_pixel, (frame->width * frame->height) >> 2);
  }
  
  const int offset = frame->width * (data->y*LCU_WIDTH);
  const int offset_c = frame->width/2 * (data->y*LCU_WIDTH_C);
  int num_pixels = frame->width * (LCU_WIDTH + 2);
  
  if (num_pixels + offset > frame->width * frame->height) {
    num_pixels = frame->width * frame->height - offset;
  }
  
  memcpy(&new_y_data[offset], &frame->rec->y[offset], sizeof(kvz_pixel) * num_pixels);
  if (frame->rec->chroma_format != KVZ_CSP_400) {
    memcpy(&new_u_data[offset_c], &frame->rec->u[offset_c], sizeof(kvz_pixel) * num_pixels >> 2);
    memcpy(&new_v_data[offset_c], &frame->rec->v[offset_c], sizeof(kvz_pixel) * num_pixels >> 2);
  }
  
  if (data->y>0) {
    //copy first row from buffer
    memcpy(&new_y_data[frame->width * (data->y*LCU_WIDTH-1)], &data->encoder_state->tile->hor_buf_before_sao->y[frame->width * (data->y-1)], frame->width * sizeof(kvz_pixel));
    if (frame->rec->chroma_format != KVZ_CSP_400) {
      memcpy(&new_u_data[frame->width / 2 * (data->y*LCU_WIDTH_C - 1)], &data->encoder_state->tile->hor_buf_before_sao->u[frame->width / 2 * (data->y - 1)], frame->width / 2 * sizeof(kvz_pixel));
      memcpy(&new_v_data[frame->width / 2 * (data->y*LCU_WIDTH_C - 1)], &data->encoder_state->tile->hor_buf_before_sao->v[frame->width / 2 * (data->y - 1)], frame->width / 2 * sizeof(kvz_pixel));
    }
  }

  for (x = 0; x < frame->width_in_lcu; x++) {
  // sao_do_rdo(encoder, lcu.x, lcu.y, sao_luma, sao_chroma);
    sao_info_t *sao_luma = &frame->sao_luma[data->y * stride + x];
    sao_info_t *sao_chroma = &frame->sao_chroma[data->y * stride + x];
    kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_y_data, x, data->y, sao_luma, COLOR_Y);
    if (frame->rec->chroma_format != KVZ_CSP_400) {
      kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_u_data, x, data->y, sao_chroma, COLOR_U);
      kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_v_data, x, data->y, sao_chroma, COLOR_V);
    }
  }
  
  free(new_y_data);
  free(new_u_data);
  free(new_v_data);

  free(opaque);
}


static int encoder_state_tree_is_a_chain(const encoder_state_t * const state) {
  if (!state->children[0].encoder_control) return 1;
  if (state->children[1].encoder_control) return 0;
  return encoder_state_tree_is_a_chain(&state->children[0]);
}

static void encoder_state_encode(encoder_state_t * const main_state) {
  //If we have children, encode at child level
  if (main_state->children[0].encoder_control) {
    int i=0;
    //If we have only one child, than it cannot be the last split in tree
    int node_is_the_last_split_in_tree = (main_state->children[1].encoder_control != 0);
    
    for (i=0; main_state->children[i].encoder_control; ++i) {
      encoder_state_t *sub_state = &(main_state->children[i]);
      
      if (sub_state->tile != main_state->tile) {
        const int offset_x = sub_state->tile->lcu_offset_x * LCU_WIDTH;
        const int offset_y = sub_state->tile->lcu_offset_y * LCU_WIDTH;
        const int width = MIN(sub_state->tile->frame->width_in_lcu * LCU_WIDTH, main_state->tile->frame->width - offset_x);
        const int height = MIN(sub_state->tile->frame->height_in_lcu * LCU_WIDTH, main_state->tile->frame->height - offset_y);
        
        if (sub_state->tile->frame->source) {
          kvz_image_free(sub_state->tile->frame->source);
          sub_state->tile->frame->source = NULL;
        }
        if (sub_state->tile->frame->rec) {
          kvz_image_free(sub_state->tile->frame->rec);
          sub_state->tile->frame->rec = NULL;
        }
        
        assert(!sub_state->tile->frame->source);
        assert(!sub_state->tile->frame->rec);
        sub_state->tile->frame->source = kvz_image_make_subimage(main_state->tile->frame->source, offset_x, offset_y, width, height);
        sub_state->tile->frame->rec = kvz_image_make_subimage(main_state->tile->frame->rec, offset_x, offset_y, width, height);
      }
      
      //To be the last split, we require that every child is a chain
      node_is_the_last_split_in_tree = node_is_the_last_split_in_tree && encoder_state_tree_is_a_chain(&main_state->children[i]);
    }
    //If it's the latest split point
    if (node_is_the_last_split_in_tree) {
      for (i=0; main_state->children[i].encoder_control; ++i) {
        //If we don't have wavefronts, parallelize encoding of children.
        if (main_state->children[i].type != ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
#ifdef KVZ_DEBUG
          char job_description[256];
          switch (main_state->children[i].type) {
            case ENCODER_STATE_TYPE_TILE: 
              sprintf(job_description, "type=encode_child,frame=%d,tile=%d,row=%d-%d,px_x=%d-%d,px_y=%d-%d", main_state->children[i].frame->num, main_state->children[i].tile->id, main_state->children[i].lcu_order[0].position.y + main_state->children[i].tile->lcu_offset_y, main_state->children[i].lcu_order[0].position.y + main_state->children[i].tile->lcu_offset_y, 
                      main_state->children[i].lcu_order[0].position_px.x + main_state->children[i].tile->lcu_offset_x * LCU_WIDTH, main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].position_px.x + main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].size.x + main_state->children[i].tile->lcu_offset_x * LCU_WIDTH - 1,
                      main_state->children[i].lcu_order[0].position_px.y + main_state->children[i].tile->lcu_offset_y * LCU_WIDTH, main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].position_px.y + main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].size.y + main_state->children[i].tile->lcu_offset_y * LCU_WIDTH - 1);
              break;
            case ENCODER_STATE_TYPE_SLICE:
              sprintf(job_description, "type=encode_child,frame=%d,slice=%d,start_in_ts=%d", main_state->children[i].frame->num, main_state->children[i].slice->id, main_state->children[i].slice->start_in_ts);
              break;
            default:
              sprintf(job_description, "type=encode_child,frame=%d,invalid", main_state->children[i].frame->num);
              break;
          }
#else
          char* job_description = NULL;
#endif
          kvz_threadqueue_free_job(&main_state->children[i].tqj_recon_done);
          main_state->children[i].tqj_recon_done = kvz_threadqueue_submit(main_state->encoder_control->threadqueue, encoder_state_worker_encode_children, &(main_state->children[i]), 1, job_description);
          if (main_state->children[i].previous_encoder_state != &main_state->children[i] && main_state->children[i].previous_encoder_state->tqj_recon_done && !main_state->children[i].frame->is_idr_frame) {
#if 0
            // Disabled due to non-determinism.
            if (main_state->encoder_control->cfg->mv_constraint == KVZ_MV_CONSTRAIN_FRAME_AND_TILE_MARGIN)
            {
              // When MV's don't cross tile boundaries, add dependancy only to the same tile.
              kvz_threadqueue_job_dep_add(main_state->children[i].tqj_recon_done, main_state->children[i].previous_encoder_state->tqj_recon_done);
            } else 
#endif      
            {
              // Add dependancy to each child in the previous frame.
              for (int child_id = 0; main_state->children[child_id].encoder_control; ++child_id) {
                kvz_threadqueue_job_dep_add(main_state->children[i].tqj_recon_done, main_state->children[child_id].previous_encoder_state->tqj_recon_done);
              }
            }
          }
          kvz_threadqueue_job_unwait_job(main_state->encoder_control->threadqueue, main_state->children[i].tqj_recon_done);
        } else {
          //Wavefront rows have parallelism at LCU level, so we should not launch multiple threads here!
          //FIXME: add an assert: we can only have wavefront children
          encoder_state_worker_encode_children(&(main_state->children[i]));
        }
      }
      
      // Add SAO reconstruction jobs and their dependancies when using WPP coding.
      if (main_state->encoder_control->cfg.sao_enable && 
          main_state->children[0].type == ENCODER_STATE_TYPE_WAVEFRONT_ROW)
      {
        videoframe_t * const frame = main_state->tile->frame;
        threadqueue_job_t *previous_job = NULL;
        
        for (int y = 0; y < frame->height_in_lcu; ++y) {
          // Queue a single job performing SAO reconstruction for the whole wavefront row.

          worker_sao_reconstruct_lcu_data *data = MALLOC(worker_sao_reconstruct_lcu_data, 1);
          threadqueue_job_t *job;
#ifdef KVZ_DEBUG
          char job_description[256];
          sprintf(job_description, "type=sao,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d", main_state->frame->num, main_state->tile->id, main_state->tile->lcu_offset_x * LCU_WIDTH, main_state->tile->lcu_offset_x * LCU_WIDTH + main_state->tile->frame->width - 1, (main_state->tile->lcu_offset_y + y) * LCU_WIDTH, MIN(main_state->tile->lcu_offset_y * LCU_WIDTH + main_state->tile->frame->height, (main_state->tile->lcu_offset_y + y + 1) * LCU_WIDTH)-1);
#else
          char* job_description = NULL;
#endif
          data->y = y;
          data->encoder_state = main_state;
          
          job = kvz_threadqueue_submit(main_state->encoder_control->threadqueue, encoder_state_worker_sao_reconstruct_lcu, data, 1, job_description);

          // If job object was returned, add dependancies and allow it to run.
          if (job) {
            // This dependancy is needed, because the pre-SAO pixels from the LCU row
            // below this one are read straigh from the frame.
            if (previous_job) {
              kvz_threadqueue_job_dep_add(job, previous_job);
            }
            previous_job = job;
            
            // This depepndancy ensures that the bottom edge of this LCU row
            // has been fully deblocked.
            if (y < frame->height_in_lcu - 1) {
              // Not last row: depend on the last LCU of the row below.
              kvz_threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 1) * frame->width_in_lcu + frame->width_in_lcu - 1]);
            } else {
              // Last row: depend on the last LCU of the row
              kvz_threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 0) * frame->width_in_lcu + frame->width_in_lcu - 1]);
            }
            kvz_threadqueue_job_unwait_job(main_state->encoder_control->threadqueue, job);
            
            // The wavefront row is finished, when the SAO-reconstruction is
            // finished.
            kvz_threadqueue_free_job(&main_state->children[y].tqj_recon_done);
            main_state->children[y].tqj_recon_done = job;
            
            if (y == frame->height_in_lcu - 1) {
              // This tile is finished, when the reconstruction of the last
              // WPP-row is finished.
              assert(!main_state->tqj_recon_done);
              main_state->tqj_recon_done = kvz_threadqueue_copy_ref(job);
            }
          }
        }
      }
    } else {
      for (i=0; main_state->children[i].encoder_control; ++i) {
        encoder_state_worker_encode_children(&(main_state->children[i]));
      }
    }
  } else {
    switch (main_state->type) {
      case ENCODER_STATE_TYPE_TILE:
      case ENCODER_STATE_TYPE_SLICE:
      case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
        encoder_state_encode_leaf(main_state);
        break;
      default:
        fprintf(stderr, "Unsupported leaf type %c!\n", main_state->type);
        assert(0);
    }
  }
}


static void encoder_ref_insertion_sort(const encoder_state_t *const state, uint8_t reflist[16], uint8_t length) {

  for (uint8_t i = 1; i < length; ++i) {
    const uint8_t cur_idx = reflist[i];
    const int32_t cur_poc = state->frame->ref->pocs[cur_idx];
    int8_t j = i;
    while (j > 0 && cur_poc < state->frame->ref->pocs[reflist[j - 1]]) {
      reflist[j] = reflist[j - 1];
      --j;
    }
    reflist[j] = cur_idx;
  }
}

/**
 * \brief Generate reference picture lists.
 *
 * \param state             main encoder state
 */
void kvz_encoder_create_ref_lists(const encoder_state_t *const state)
{
  // TODO check possibility to add L0 references to L1 list also
  
  FILL_ARRAY(state->frame->ref_LX_size, 0, 2);

  // List all pocs of lists
  int j = 0;
  for (j = 0; j < state->frame->ref->used_size; j++) {
    if (state->frame->ref->pocs[j] < state->frame->poc) {
      state->frame->ref_LX[0][state->frame->ref_LX_size[0]] = j;
      state->frame->ref_LX_size[0] += 1;
    } else {
      state->frame->ref_LX[1][state->frame->ref_LX_size[1]] = j;
      state->frame->ref_LX_size[1] += 1;
    }
  }

  // Fill the rest with -1s.
  for (; j < 16; j++) {
    state->frame->ref_LX[0][j] = (uint8_t) -1;
    state->frame->ref_LX[1][j] = (uint8_t) -1;
  }

  //encoder_ref_insertion_sort(state, state->frame->ref_LX[0], state->frame->ref_LX_size[0]);
  //encoder_ref_insertion_sort(state, state->frame->ref_LX[1], state->frame->ref_LX_size[1]);
}

/**
 * \brief Remove any references that should no longer be used.
 */
static void encoder_state_remove_refs(encoder_state_t *state) {
  const encoder_control_t * const encoder = state->encoder_control;
  
  int neg_refs = encoder->cfg.gop[state->frame->gop_offset].ref_neg_count;
  int pos_refs = encoder->cfg.gop[state->frame->gop_offset].ref_pos_count;

  unsigned target_ref_num;
  if (encoder->cfg.gop_len) {
    target_ref_num = neg_refs + pos_refs;
  } else {
    target_ref_num = encoder->cfg.ref_frames;
  }
  if (state->frame->slicetype == KVZ_SLICE_I) {
    target_ref_num = 0;
  }

  if (encoder->cfg.gop_len && target_ref_num > 0) {
    // With GOP in use, go through all the existing reference pictures and
    // remove any picture that is not referenced by the current picture.

    for (int ref = state->frame->ref->used_size - 1; ref >= 0; --ref) {
      bool is_referenced = false;

      int ref_poc = state->frame->ref->pocs[ref];
      
      for (int i = 0; i < neg_refs; i++) {
        int ref_relative_poc = -encoder->cfg.gop[state->frame->gop_offset].ref_neg[i];
        if (ref_poc == state->frame->poc + ref_relative_poc) {
          is_referenced = true;
          break;
        }
      }

      
      for (int i = 0; i < pos_refs; i++) {
        int ref_relative_poc = encoder->cfg.gop[state->frame->gop_offset].ref_pos[i];
        if (ref_poc == state->frame->poc + ref_relative_poc) {
          is_referenced = true;
          break;
        }
      }

      if (!is_referenced) {
        // This reference is not referred to by this frame, it must be removed.
        kvz_image_list_rem(state->frame->ref, ref);
      }
    }
  } else {
    // Without GOP, remove the oldest picture.
    while (state->frame->ref->used_size > target_ref_num) {
      int8_t oldest_ref = state->frame->ref->used_size - 1;
      kvz_image_list_rem(state->frame->ref, oldest_ref);
    }
  }

  assert(state->frame->ref->used_size <= target_ref_num);
}

static void encoder_state_reset_poc(encoder_state_t *state) {
  state->frame->poc = 0;
  kvz_videoframe_set_poc(state->tile->frame, 0);

  for (int i = 0; state->children[i].encoder_control; ++i) {
    encoder_state_t *sub_state = &(state->children[i]);
    encoder_state_reset_poc(sub_state);
  }
}

static void encoder_set_source_picture(encoder_state_t * const state, kvz_picture* frame)
{
  assert(!state->tile->frame->source);
  assert(!state->tile->frame->rec);

  state->tile->frame->source = frame;
  if (state->encoder_control->cfg.lossless) {
    // In lossless mode, the reconstruction is equal to the source frame.
    state->tile->frame->rec = kvz_image_copy_ref(frame);
  } else {
    state->tile->frame->rec = kvz_image_alloc(state->encoder_control->chroma_format, frame->width, frame->height);
    state->tile->frame->rec->dts = frame->dts;
    state->tile->frame->rec->pts = frame->pts;
  }

  kvz_videoframe_set_poc(state->tile->frame, state->frame->poc);
}

static void encoder_state_init_children(encoder_state_t * const state) {
  kvz_bitstream_clear(&state->stream);

  if (state->is_leaf) {
    //Leaf states have cabac and context
    kvz_cabac_start(&state->cabac);
    kvz_init_contexts(state, state->frame->QP, state->frame->slicetype);
  }

  //Clear the jobs
  kvz_threadqueue_free_job(&state->tqj_bitstream_written);
  kvz_threadqueue_free_job(&state->tqj_recon_done);

  for (int i = 0; state->children[i].encoder_control; ++i) {
    encoder_state_init_children(&state->children[i]);
  }
}

static void normalize_lcu_weights(encoder_state_t * const state)
{
  if (state->frame->num == 0) return;

  const uint32_t num_lcus = state->encoder_control->in.width_in_lcu *
                            state->encoder_control->in.height_in_lcu;
  double sum = 0.0;
  for (uint32_t i = 0; i < num_lcus; i++) {
    sum += state->frame->lcu_stats[i].weight;
  }

  for (uint32_t i = 0; i < num_lcus; i++) {
    state->frame->lcu_stats[i].weight /= sum;
  }
}

static void encoder_state_init_new_frame(encoder_state_t * const state, kvz_picture* frame) {
  assert(state->type == ENCODER_STATE_TYPE_MAIN);

  const kvz_config * const cfg = &state->encoder_control->cfg;

  encoder_set_source_picture(state, frame);

  // Check whether the frame is a keyframe or not.
  if (state->frame->num == 0) {
    state->frame->is_idr_frame = true;
  } else {
    bool is_i_idr = (cfg->intra_period == 1 && state->frame->num % 2 == 0);
    bool is_p_idr = (cfg->intra_period > 1 && (state->frame->num % cfg->intra_period) == 0);
    state->frame->is_idr_frame = is_i_idr || is_p_idr;
  }

  // Set pictype.
  if (state->frame->is_idr_frame) {
    state->frame->pictype = KVZ_NAL_IDR_W_RADL;
  } else {
    state->frame->pictype = KVZ_NAL_TRAIL_R;
  }

  // Set slicetype.
  if (state->frame->is_idr_frame || cfg->intra_period == 1) {
    state->frame->slicetype = KVZ_SLICE_I;
  } else if (cfg->gop_len > 0 && !cfg->gop_lowdelay) {
    state->frame->slicetype = KVZ_SLICE_B;
  } else {
    state->frame->slicetype = KVZ_SLICE_P;
  }

  // Set POC.
  if (state->frame->is_idr_frame) {
    encoder_state_reset_poc(state);
  } else if (cfg->intra_period != 1 && cfg->gop_len > 0) {
    // Calculate POC according to the global frame counter and GOP
    // structure.
    int32_t poc;
    if (cfg->intra_period > 0) {
      poc = state->frame->num % cfg->intra_period - 1;
    } else {
      poc = state->frame->num - 1;
    }
    int32_t poc_offset = cfg->gop[state->frame->gop_offset].poc_offset;
    state->frame->poc = poc - poc % cfg->gop_len + poc_offset;
    kvz_videoframe_set_poc(state->tile->frame, state->frame->poc);
  }

  encoder_state_remove_refs(state);
  kvz_encoder_create_ref_lists(state);

  normalize_lcu_weights(state);
  kvz_set_picture_lambda_and_qp(state);

  encoder_state_init_children(state);
}

static void _encode_one_frame_add_bitstream_deps(const encoder_state_t * const state, threadqueue_job_t * const job) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    _encode_one_frame_add_bitstream_deps(&state->children[i], job);
  }
  if (state->tqj_bitstream_written) {
    kvz_threadqueue_job_dep_add(job, state->tqj_bitstream_written);
  }
  if (state->tqj_recon_done) {
    kvz_threadqueue_job_dep_add(job, state->tqj_recon_done);
  }
}


void kvz_encode_one_frame(encoder_state_t * const state, kvz_picture* frame)
{
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    encoder_state_init_new_frame(state, frame);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=init_new_frame,frame=%d,poc=%d", state->frame->num, state->frame->poc);
  }
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    encoder_state_encode(state);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=encode,frame=%d", state->frame->num);
  }
  //kvz_threadqueue_flush(main_state->encoder_control->threadqueue);
  {
    threadqueue_job_t *job;
#ifdef KVZ_DEBUG
    char job_description[256];
    sprintf(job_description, "type=write_bitstream,frame=%d", state->frame->num);
#else
    char* job_description = NULL;
#endif

    job = kvz_threadqueue_submit(state->encoder_control->threadqueue, kvz_encoder_state_worker_write_bitstream, (void*) state, 1, job_description);
    
    _encode_one_frame_add_bitstream_deps(state, job);
    if (state->previous_encoder_state != state && state->previous_encoder_state->tqj_bitstream_written) {
      //We need to depend on previous bitstream generation
      kvz_threadqueue_job_dep_add(job, state->previous_encoder_state->tqj_bitstream_written);
    }
    kvz_threadqueue_job_unwait_job(state->encoder_control->threadqueue, job);
    assert(!state->tqj_bitstream_written);
    state->tqj_bitstream_written = job;
  }
  state->frame->done = 0;
  //kvz_threadqueue_flush(main_state->encoder_control->threadqueue);
}


/**
 * Prepare the encoder state for encoding the next frame.
 *
 * - Add the previous reconstructed picture as a reference, if needed.
 * - Free the previous reconstructed and source pictures.
 * - Create a new cu array, if needed.
 * - Update frame count and POC.
 */
void kvz_encoder_prepare(encoder_state_t *state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  // The previous frame must be done before the next one is started.
  assert(state->frame->done);

  if (state->frame->num == -1) {
    // We're at the first frame, so don't care about all this stuff.
    state->frame->num = 0;
    state->frame->poc   = 0;
    assert(!state->tile->frame->source);
    assert(!state->tile->frame->rec);
    state->frame->prepared = 1;
    return;
  }

  // NOTE: prev_state is equal to state when OWF is zero
  encoder_state_t *prev_state = state->previous_encoder_state;

  if (state->previous_encoder_state != state) {
    kvz_cu_array_free(state->tile->frame->cu_array);
    state->tile->frame->cu_array = NULL;
    unsigned width  = state->tile->frame->width_in_lcu  * LCU_WIDTH;
    unsigned height = state->tile->frame->height_in_lcu * LCU_WIDTH;
    state->tile->frame->cu_array = kvz_cu_array_alloc(width, height);

    kvz_image_list_copy_contents(state->frame->ref, prev_state->frame->ref);
  }

  if (!encoder->cfg.gop_len ||
      !prev_state->frame->poc ||
      encoder->cfg.gop[prev_state->frame->gop_offset].is_ref) {

    // Store current list of POCs for use in TMVP derivation
    memcpy(prev_state->tile->frame->rec->ref_pocs, state->frame->ref->pocs, sizeof(int32_t)*state->frame->ref->used_size);

    // Add previous reconstructed picture as a reference
    kvz_image_list_add(state->frame->ref,
                   prev_state->tile->frame->rec,
                   prev_state->tile->frame->cu_array,
                   prev_state->frame->poc);
    kvz_cu_array_free(state->tile->frame->cu_array);
    unsigned height = state->tile->frame->height_in_lcu * LCU_WIDTH;
    unsigned width  = state->tile->frame->width_in_lcu  * LCU_WIDTH;
    state->tile->frame->cu_array = kvz_cu_array_alloc(width, height);
  }

  // Remove source and reconstructed picture.
  kvz_image_free(state->tile->frame->source);
  state->tile->frame->source = NULL;
  kvz_image_free(state->tile->frame->rec);
  state->tile->frame->rec = NULL;

  // Update POC and frame count.
  state->frame->num = prev_state->frame->num + 1;
  state->frame->poc   = prev_state->frame->poc   + 1;

  state->frame->prepared = 1;
}

coeff_scan_order_t kvz_get_scan_order(int8_t cu_type, int intra_mode, int depth)
{
  // Scan mode is diagonal, except for 4x4+8x8 luma and 4x4 chroma, where:
  // - angular 6-14 = vertical
  // - angular 22-30 = horizontal
  if (cu_type == CU_INTRA && depth >= 3) {
    if (intra_mode >= 6 && intra_mode <= 14) {
      return SCAN_VER;
    } else if (intra_mode >= 22 && intra_mode <= 30) {
      return SCAN_HOR;
    }
  }

  return SCAN_DIAG;
}

lcu_stats_t* kvz_get_lcu_stats(encoder_state_t *state, int lcu_x, int lcu_y)
{
  const int index = lcu_x + state->tile->lcu_offset_x +
                    (lcu_y + state->tile->lcu_offset_y) *
                    state->encoder_control->in.width_in_lcu;
  return &state->frame->lcu_stats[index];
}
