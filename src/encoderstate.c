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

/*
 * \file
 */

#include "encoderstate.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tables.h"
#include "config.h"
#include "cabac.h"
#include "image.h"
#include "nal.h"
#include "context.h"
#include "transform.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "search.h"
#include "sao.h"
#include "rdo.h"
#include "rate_control.h"
#include "strategies/strategies-picture.h"

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
    const int rdpx = lcu->position_px.x;
    const int rdpy = lcu->position_px.y + lcu->size.y - 1;
    const int by = lcu->position.y;
    
    //Copy the bottom row of this LCU to the horizontal buffer
    kvz_pixels_blit(&frame->rec->y[rdpy * frame->rec->stride + rdpx],
                        &hor_buf->y[lcu->position_px.x + by * frame->width],
                        lcu->size.x, 1, frame->rec->stride, frame->width);
    kvz_pixels_blit(&frame->rec->u[(rdpy/2) * frame->rec->stride/2 + (rdpx/2)],
                        &hor_buf->u[lcu->position_px.x / 2 + by * frame->width / 2],
                        lcu->size.x / 2, 1, frame->rec->stride / 2, frame->width / 2);
    kvz_pixels_blit(&frame->rec->v[(rdpy/2) * frame->rec->stride/2 + (rdpx/2)],
                        &hor_buf->v[lcu->position_px.x / 2 + by * frame->width / 2],
                        lcu->size.x / 2, 1, frame->rec->stride / 2, frame->width / 2);
  }
  
  if (ver_buf) {
    const int rdpx = lcu->position_px.x + lcu->size.x - 1;
    const int rdpy = lcu->position_px.y;
    const int bx = lcu->position.x;
    
    
    //Copy the right row of this LCU to the vertical buffer.
    kvz_pixels_blit(&frame->rec->y[rdpy * frame->rec->stride + rdpx],
                        &ver_buf->y[lcu->position_px.y + bx * frame->height],
                        1, lcu->size.y, frame->rec->stride, 1);
    kvz_pixels_blit(&frame->rec->u[(rdpy/2) * frame->rec->stride/2 + (rdpx/2)],
                        &ver_buf->u[lcu->position_px.y / 2 + bx * frame->height / 2],
                        1, lcu->size.y / 2, frame->rec->stride / 2, 1);
    kvz_pixels_blit(&frame->rec->v[(rdpy/2) * frame->rec->stride/2 + (rdpx/2)],
                        &ver_buf->v[lcu->position_px.y / 2 + bx * frame->height / 2],
                        1, lcu->size.y / 2, frame->rec->stride / 2, 1);
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
    encode_sao_color(state, sao_chroma, COLOR_U);
    encode_sao_color(state, sao_chroma, COLOR_V);
  }
}


static void encoder_state_worker_encode_lcu(void * opaque) {
  const lcu_order_element_t * const lcu = opaque;
  encoder_state_t *state = lcu->encoder_state;
  const encoder_control_t * const encoder = state->encoder_control;
  videoframe_t* const frame = state->tile->frame;
  
  //This part doesn't write to bitstream, it's only search, deblock and sao
  
  kvz_search_lcu(state, lcu->position_px.x, lcu->position_px.y, state->tile->hor_buf_search, state->tile->ver_buf_search);
    
  encoder_state_recdata_to_bufs(state, lcu, state->tile->hor_buf_search, state->tile->ver_buf_search);

  if (encoder->deblock_enable) {
    kvz_filter_deblock_lcu(state, lcu->position_px.x, lcu->position_px.y);
  }

  if (encoder->sao_enable) {
    const int stride = frame->width_in_lcu;
    int32_t merge_cost_luma[3] = { INT32_MAX };
    int32_t merge_cost_chroma[3] = { INT32_MAX };
    sao_info_t *sao_luma = &frame->sao_luma[lcu->position.y * stride + lcu->position.x];
    sao_info_t *sao_chroma = &frame->sao_chroma[lcu->position.y * stride + lcu->position.x];

    // Merge candidates
    sao_info_t *sao_top_luma = lcu->position.y != 0 ? &frame->sao_luma[(lcu->position.y - 1) * stride + lcu->position.x] : NULL;
    sao_info_t *sao_left_luma = lcu->position.x != 0 ? &frame->sao_luma[lcu->position.y * stride + lcu->position.x - 1] : NULL;
    sao_info_t *sao_top_chroma = lcu->position.y != 0 ? &frame->sao_chroma[(lcu->position.y - 1) * stride + lcu->position.x] : NULL;
    sao_info_t *sao_left_chroma = lcu->position.x != 0 ? &frame->sao_chroma[lcu->position.y * stride + lcu->position.x - 1] : NULL;

    kvz_sao_search_luma(state, frame, lcu->position.x, lcu->position.y, sao_luma, sao_top_luma, sao_left_luma, merge_cost_luma);
    kvz_sao_search_chroma(state, frame, lcu->position.x, lcu->position.y, sao_chroma, sao_top_chroma, sao_left_chroma, merge_cost_chroma);

    sao_luma->merge_up_flag = sao_luma->merge_left_flag = 0;
    // Check merge costs
    if (sao_top_luma) {
      // Merge up if cost is equal or smaller to the searched mode cost
      if (merge_cost_luma[2] + merge_cost_chroma[2] <= merge_cost_luma[0] + merge_cost_chroma[0]) {        
        *sao_luma = *sao_top_luma;
        *sao_chroma = *sao_top_chroma;
        sao_luma->merge_up_flag = 1;
        sao_luma->merge_left_flag = 0;
      }
    }
    if (sao_left_luma) {
      // Merge left if cost is equal or smaller to the searched mode cost 
      // AND smaller than merge up cost, if merge up was already chosen
      if (merge_cost_luma[1] + merge_cost_chroma[1] <= merge_cost_luma[0] + merge_cost_chroma[0]) {
        if (!sao_luma->merge_up_flag || merge_cost_luma[1] + merge_cost_chroma[1] < merge_cost_luma[2] + merge_cost_chroma[2]) {      
          *sao_luma = *sao_left_luma;
          *sao_chroma = *sao_left_chroma;
          sao_luma->merge_left_flag = 1;
          sao_luma->merge_up_flag = 0;
        }
      }
    }
    assert(sao_luma->eo_class < SAO_NUM_EO);
    assert(sao_chroma->eo_class < SAO_NUM_EO);
    
    CHECKPOINT_SAO_INFO("sao_luma", *sao_luma);
    CHECKPOINT_SAO_INFO("sao_chroma", *sao_chroma);
  }
  
  // Copy LCU cu_array to main states cu_array, because that is the only one
  // which is given to the next frame through image_list_t.
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);

    encoder_state_t *main_state = state;
    while (main_state->parent) main_state = main_state->parent;
    assert(main_state != state);

    unsigned child_width_in_scu = state->tile->frame->width_in_lcu << MAX_DEPTH;
    unsigned main_width_in_scu = main_state->tile->frame->width_in_lcu << MAX_DEPTH;
    unsigned tile_x = state->tile->lcu_offset_x;
    unsigned tile_y = state->tile->lcu_offset_y;

    unsigned x = lcu->position.x << MAX_DEPTH;
    unsigned y = lcu->position.y << MAX_DEPTH;

    for (unsigned lcu_row = 0; lcu_row < 8; ++lcu_row) {
      cu_info_t *main_row = &main_state->tile->frame->cu_array->data[x + tile_x + (y + tile_y + lcu_row) * main_width_in_scu];
      cu_info_t *child_row = &state->tile->frame->cu_array->data[x + (y + lcu_row) * child_width_in_scu];
      memcpy(main_row, child_row, sizeof(cu_info_t) * 8);
    }

    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=copy_cuinfo,frame=%d,tile=%d", state->global->frame, state->tile->id);
  }
  
  //Now write data to bitstream (required to have a correct CABAC state)
  
  //First LCU, and we are in a slice. We need a slice header
  if (state->type == ENCODER_STATE_TYPE_SLICE && lcu->index == 0) {
    kvz_encoder_state_write_bitstream_slice_header(state);
    kvz_bitstream_add_rbsp_trailing_bits(&state->stream); 
  }
  
  //Encode SAO
  if (encoder->sao_enable) {
    encode_sao(state, lcu->position.x, lcu->position.y, &frame->sao_luma[lcu->position.y * frame->width_in_lcu + lcu->position.x], &frame->sao_chroma[lcu->position.y * frame->width_in_lcu + lcu->position.x]);
  }
  
  //Encode coding tree
  kvz_encode_coding_tree(state, lcu->position.x << MAX_DEPTH, lcu->position.y << MAX_DEPTH, 0);

  //Terminator
  if (lcu->index < state->lcu_order_count - 1) {
    //Since we don't handle slice segments, end of slice segment == end of slice
    //Always 0 since otherwise it would be split
    kvz_cabac_encode_bin_trm(&state->cabac, 0);  // end_of_slice_segment_flag
  }
  
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
  
  if (encoder->sao_enable && lcu->above) {
    //If we're not the first in the row
    if (lcu->above->left) {
      encoder_state_recdata_to_bufs(state, lcu->above->left, state->tile->hor_buf_before_sao, NULL);
    }
    //Latest LCU in the row, copy the data from the one above also
    if (!lcu->right) {
      encoder_state_recdata_to_bufs(state, lcu->above, state->tile->hor_buf_before_sao, NULL);
    }
  }
}

static void encoder_state_encode_leaf(encoder_state_t * const state) {
  assert(state->is_leaf);
  assert(state->lcu_order_count > 0);
  
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
        PERFORMANCE_MEASURE_END(KVZ_PERF_LCU, state->encoder_control->threadqueue, "type=encode_lcu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", state->global->frame, state->tile->id, state->slice->id, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH + lcu->size.x - 1, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH + lcu->size.y - 1);
      }
#endif //KVZ_DEBUG
    }
    
    if (state->encoder_control->sao_enable) {
      PERFORMANCE_MEASURE_START(KVZ_PERF_SAOREC);
      kvz_sao_reconstruct_frame(state);
      PERFORMANCE_MEASURE_END(KVZ_PERF_SAOREC, state->encoder_control->threadqueue, "type=kvz_sao_reconstruct_frame,frame=%d,tile=%d,slice=%d,row=%d-%d,px_x=%d-%d,px_y=%d-%d", state->global->frame, state->tile->id, state->slice->id, state->lcu_order[0].position.y + state->tile->lcu_offset_y, state->lcu_order[state->lcu_order_count - 1].position.y + state->tile->lcu_offset_y,
        state->tile->lcu_offset_x * LCU_WIDTH, state->tile->frame->width + state->tile->lcu_offset_x * LCU_WIDTH - 1,
        state->tile->lcu_offset_y * LCU_WIDTH, state->tile->frame->height + state->tile->lcu_offset_y * LCU_WIDTH - 1
      );
    }
  } else {
    // Add each LCU in the wavefront row as it's own job to the queue.

    for (int i = 0; i < state->lcu_order_count; ++i) {
      const lcu_order_element_t * const lcu = &state->lcu_order[i];
#ifdef KVZ_DEBUG
      char job_description[256];
      sprintf(job_description, "type=encode_lcu,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", state->global->frame, state->tile->id, state->slice->id, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH, lcu->position_px.x + state->tile->lcu_offset_x * LCU_WIDTH + lcu->size.x - 1, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH, lcu->position_px.y + state->tile->lcu_offset_y * LCU_WIDTH + lcu->size.y - 1);
#else
      char* job_description = NULL;
#endif
      state->tile->wf_jobs[lcu->id] = kvz_threadqueue_submit(state->encoder_control->threadqueue, encoder_state_worker_encode_lcu, (void*)lcu, 1, job_description);
      
      // If job object was returned, add dependancies and allow it to run.
      if (state->tile->wf_jobs[lcu->id]) {
        // Add inter frame dependancies when ecoding more than one frame at
        // once. The added dependancy is for the first LCU of each wavefront
        // row to depend on the reconstruction status of the row below in the
        // previous frame.
        if (state->previous_encoder_state != state && state->previous_encoder_state->tqj_recon_done && state->global->slicetype != KVZ_SLICE_I) {
          if (!lcu->left) {
            if (lcu->below) {
              kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], lcu->below->encoder_state->previous_encoder_state->tqj_recon_done);
            } else {
              kvz_threadqueue_job_dep_add(state->tile->wf_jobs[lcu->id], lcu->encoder_state->previous_encoder_state->tqj_recon_done);
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
      }

      // In the case where SAO is not enabled, the wavefront row is
      // done when the last LCU in the row is done.
      if (!state->encoder_control->sao_enable && i + 1 == state->lcu_order_count) {
        assert(!state->tqj_recon_done);
        state->tqj_recon_done = state->tile->wf_jobs[lcu->id];
      }
    }
  }
}

static void encoder_state_encode(encoder_state_t * const main_state);

static void encoder_state_worker_encode_children(void * opaque) {
  encoder_state_t *sub_state = opaque;
  encoder_state_encode(sub_state);
  if (sub_state->is_leaf) {
    if (sub_state->type != ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
      PERFORMANCE_MEASURE_START(KVZ_PERF_BSLEAF);
      kvz_encoder_state_write_bitstream_leaf(sub_state);
      PERFORMANCE_MEASURE_END(KVZ_PERF_BSLEAF, sub_state->encoder_control->threadqueue, "type=encoder_state_write_bitstream_leaf,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", sub_state->global->frame, sub_state->tile->id, sub_state->slice->id, sub_state->lcu_order[0].position_px.x + sub_state->tile->lcu_offset_x * LCU_WIDTH, sub_state->lcu_order[sub_state->lcu_order_count - 1].position_px.x + sub_state->lcu_order[sub_state->lcu_order_count - 1].size.x + sub_state->tile->lcu_offset_x * LCU_WIDTH - 1, sub_state->lcu_order[0].position_px.y + sub_state->tile->lcu_offset_y * LCU_WIDTH, sub_state->lcu_order[sub_state->lcu_order_count - 1].position_px.y + sub_state->lcu_order[sub_state->lcu_order_count - 1].size.y + sub_state->tile->lcu_offset_y * LCU_WIDTH - 1);
    } else {
      threadqueue_job_t *job;
#ifdef KVZ_DEBUG
      char job_description[256];
      sprintf(job_description, "type=encoder_state_write_bitstream_leaf,frame=%d,tile=%d,slice=%d,px_x=%d-%d,px_y=%d-%d", sub_state->global->frame, sub_state->tile->id, sub_state->slice->id, sub_state->lcu_order[0].position_px.x + sub_state->tile->lcu_offset_x * LCU_WIDTH, sub_state->lcu_order[sub_state->lcu_order_count-1].position_px.x + sub_state->lcu_order[sub_state->lcu_order_count-1].size.x + sub_state->tile->lcu_offset_x * LCU_WIDTH - 1, sub_state->lcu_order[0].position_px.y + sub_state->tile->lcu_offset_y * LCU_WIDTH, sub_state->lcu_order[sub_state->lcu_order_count-1].position_px.y + sub_state->lcu_order[sub_state->lcu_order_count-1].size.y + sub_state->tile->lcu_offset_y * LCU_WIDTH - 1);
#else
      char* job_description = NULL;
#endif
      job = kvz_threadqueue_submit(sub_state->encoder_control->threadqueue, kvz_encoder_state_worker_write_bitstream_leaf, sub_state, 1, job_description);
      kvz_threadqueue_job_dep_add(job, sub_state->tile->wf_jobs[sub_state->wfrow->lcu_offset_y * sub_state->tile->frame->width_in_lcu + sub_state->lcu_order_count - 1]);
      kvz_threadqueue_job_unwait_job(sub_state->encoder_control->threadqueue, job);
      
      assert(!sub_state->tqj_bitstream_written);
      //Bitstream is written for the row, if we're at the last LCU
      sub_state->tqj_bitstream_written = job;
      return;
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
  kvz_pixel *new_u_data = MALLOC(kvz_pixel, (frame->width * frame->height) >> 2);
  kvz_pixel *new_v_data = MALLOC(kvz_pixel, (frame->width * frame->height) >> 2);
  
  const int offset = frame->width * (data->y*LCU_WIDTH);
  const int offset_c = frame->width/2 * (data->y*LCU_WIDTH_C);
  int num_pixels = frame->width * (LCU_WIDTH + 2);
  
  if (num_pixels + offset > frame->width * frame->height) {
    num_pixels = frame->width * frame->height - offset;
  }
  
  memcpy(&new_y_data[offset], &frame->rec->y[offset], sizeof(kvz_pixel) * num_pixels);
  memcpy(&new_u_data[offset_c], &frame->rec->u[offset_c], sizeof(kvz_pixel) * num_pixels >> 2);
  memcpy(&new_v_data[offset_c], &frame->rec->v[offset_c], sizeof(kvz_pixel) * num_pixels >> 2);
  
  if (data->y>0) {
    //copy first row from buffer
    memcpy(&new_y_data[frame->width * (data->y*LCU_WIDTH-1)], &data->encoder_state->tile->hor_buf_before_sao->y[frame->width * (data->y-1)], frame->width * sizeof(kvz_pixel));
    memcpy(&new_u_data[frame->width/2 * (data->y*LCU_WIDTH_C-1)], &data->encoder_state->tile->hor_buf_before_sao->u[frame->width/2 * (data->y-1)], frame->width/2 * sizeof(kvz_pixel));
    memcpy(&new_v_data[frame->width/2 * (data->y*LCU_WIDTH_C-1)], &data->encoder_state->tile->hor_buf_before_sao->v[frame->width/2 * (data->y-1)], frame->width/2 * sizeof(kvz_pixel));
  }

  for (x = 0; x < frame->width_in_lcu; x++) {
  // sao_do_rdo(encoder, lcu.x, lcu.y, sao_luma, sao_chroma);
    sao_info_t *sao_luma = &frame->sao_luma[data->y * stride + x];
    sao_info_t *sao_chroma = &frame->sao_chroma[data->y * stride + x];
    kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_y_data, x, data->y, sao_luma, COLOR_Y);
    kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_u_data, x, data->y, sao_chroma, COLOR_U);
    kvz_sao_reconstruct(data->encoder_state->encoder_control, frame, new_v_data, x, data->y, sao_chroma, COLOR_V);
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
              sprintf(job_description, "type=encode_child,frame=%d,tile=%d,row=%d-%d,px_x=%d-%d,px_y=%d-%d", main_state->children[i].global->frame, main_state->children[i].tile->id, main_state->children[i].lcu_order[0].position.y + main_state->children[i].tile->lcu_offset_y, main_state->children[i].lcu_order[0].position.y + main_state->children[i].tile->lcu_offset_y, 
                      main_state->children[i].lcu_order[0].position_px.x + main_state->children[i].tile->lcu_offset_x * LCU_WIDTH, main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].position_px.x + main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].size.x + main_state->children[i].tile->lcu_offset_x * LCU_WIDTH - 1,
                      main_state->children[i].lcu_order[0].position_px.y + main_state->children[i].tile->lcu_offset_y * LCU_WIDTH, main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].position_px.y + main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].size.y + main_state->children[i].tile->lcu_offset_y * LCU_WIDTH - 1);
              break;
            case ENCODER_STATE_TYPE_SLICE:
              sprintf(job_description, "type=encode_child,frame=%d,slice=%d,start_in_ts=%d", main_state->children[i].global->frame, main_state->children[i].slice->id, main_state->children[i].slice->start_in_ts);
              break;
            default:
              sprintf(job_description, "type=encode_child,frame=%d,invalid", main_state->children[i].global->frame);
              break;
          }
#else
          char* job_description = NULL;
#endif
          main_state->children[i].tqj_recon_done = kvz_threadqueue_submit(main_state->encoder_control->threadqueue, encoder_state_worker_encode_children, &(main_state->children[i]), 1, job_description);
          if (main_state->children[i].previous_encoder_state != &main_state->children[i] && main_state->children[i].previous_encoder_state->tqj_recon_done && !main_state->children[i].global->is_idr_frame) {
            // Add dependancy to each child in the previous frame.
            // TODO: Make it so that only adjacent tiles are dependet upon and search is constrained to those?
            for (int child_id = 0; main_state->children[child_id].encoder_control; ++child_id) {
              kvz_threadqueue_job_dep_add(main_state->children[i].tqj_recon_done, main_state->children[child_id].previous_encoder_state->tqj_recon_done);
            }
          }
          kvz_threadqueue_job_unwait_job(main_state->encoder_control->threadqueue, main_state->children[i].tqj_recon_done);
        } else {
          //Wavefront rows have parallelism at LCU level, so we should not launch multiple threads here!
          //FIXME: add an assert: we can only have wavefront children
          encoder_state_worker_encode_children(&(main_state->children[i]));
        }
      }
      
      //If children are wavefront, we need to reconstruct SAO
      if (main_state->encoder_control->sao_enable && main_state->children[0].type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
        int y;
        videoframe_t * const frame = main_state->tile->frame;
        threadqueue_job_t *previous_job = NULL;
        
        for (y = 0; y < frame->height_in_lcu; ++y) {
          worker_sao_reconstruct_lcu_data *data = MALLOC(worker_sao_reconstruct_lcu_data, 1);
          threadqueue_job_t *job;
#ifdef KVZ_DEBUG
          char job_description[256];
          sprintf(job_description, "type=sao,frame=%d,tile=%d,px_x=%d-%d,px_y=%d-%d", main_state->global->frame, main_state->tile->id, main_state->tile->lcu_offset_x * LCU_WIDTH, main_state->tile->lcu_offset_x * LCU_WIDTH + main_state->tile->frame->width - 1, (main_state->tile->lcu_offset_y + y) * LCU_WIDTH, MIN(main_state->tile->lcu_offset_y * LCU_WIDTH + main_state->tile->frame->height, (main_state->tile->lcu_offset_y + y + 1) * LCU_WIDTH)-1);
#else
          char* job_description = NULL;
#endif
          data->y = y;
          data->encoder_state = main_state;
          
          job = kvz_threadqueue_submit(main_state->encoder_control->threadqueue, encoder_state_worker_sao_reconstruct_lcu, data, 1, job_description);
          
          if (previous_job) {
            kvz_threadqueue_job_dep_add(job, previous_job);
          }
          previous_job = job;
          
          if (y < frame->height_in_lcu - 1) {
            //Not last row: depend on the last LCU of the row below
            kvz_threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 1) * frame->width_in_lcu + frame->width_in_lcu - 1]);
          } else {
            //Last row: depend on the last LCU of the row
            kvz_threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 0) * frame->width_in_lcu + frame->width_in_lcu - 1]);
          }
          kvz_threadqueue_job_unwait_job(main_state->encoder_control->threadqueue, job);
          
          //Set wfrow recon job
          main_state->children[y].tqj_recon_done = job;
          
          if (y == frame->height_in_lcu - 1) {
            assert(!main_state->tqj_recon_done);
            main_state->tqj_recon_done = job;
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


static void encoder_ref_insertion_sort(int reflist[16], int length) {

  for (uint8_t i = 1; i < length; ++i) {
    const int16_t cur_poc = reflist[i];
    int16_t j = i;
    while (j > 0 && cur_poc < reflist[j - 1]) {
      reflist[j] = reflist[j - 1];
      --j;
    }
    reflist[j] = cur_poc;
  }
}

/**
 * \brief Return reference picture lists.
 *
 * \param state             main encoder state
 * \param ref_list_len_out  Returns the lengths of the reference lists.
 * \param ref_list_poc_out  Returns two lists of POCs of the reference pictures.
 */
void kvz_encoder_get_ref_lists(const encoder_state_t *const state,
                               int ref_list_len_out[2],
                               int ref_list_poc_out[2][16])
{
  FILL_ARRAY(ref_list_len_out, 0, 2);

  // List all pocs of lists
  int j = 0;
  for (j = 0; j < state->global->ref->used_size; j++) {
    if (state->global->ref->pocs[j] < state->global->poc) {
      ref_list_poc_out[0][ref_list_len_out[0]] = state->global->ref->pocs[j];
      ref_list_len_out[0]++;
    } else {
      ref_list_poc_out[1][ref_list_len_out[1]] = state->global->ref->pocs[j];
      ref_list_len_out[1]++;
    }
  }

  // Fill the rest of ref_list_poc_out array with -1s.
  for (; j < 16; j++) {
    ref_list_poc_out[0][j] = -1;
    ref_list_poc_out[1][j] = -1;
  }

  encoder_ref_insertion_sort(ref_list_poc_out[0], ref_list_len_out[0]);
  encoder_ref_insertion_sort(ref_list_poc_out[1], ref_list_len_out[1]);
}

static void encoder_state_ref_sort(encoder_state_t *state) {
  int ref_list_len[2];
  int ref_list_poc[2][16];

  kvz_encoder_get_ref_lists(state, ref_list_len, ref_list_poc);

  for (int j = 0; j < state->global->ref->used_size; j++) {
    if (state->global->ref->pocs[j] < state->global->poc) {
      for (int ref_idx = 0; ref_idx < ref_list_len[0]; ref_idx++) {
        if (ref_list_poc[0][ref_idx] == state->global->ref->pocs[j]) {
          state->global->refmap[j].idx = ref_list_len[0] - ref_idx - 1;
          break;
        }
      }
      state->global->refmap[j].list = 1;

    } else {
      for (int ref_idx = 0; ref_idx < ref_list_len[1]; ref_idx++) {
        if (ref_list_poc[1][ref_idx] == state->global->ref->pocs[j]) {
          state->global->refmap[j].idx = ref_idx;
          break;
        }
      }
      state->global->refmap[j].list = 2;
    }
    state->global->refmap[j].poc = state->global->ref->pocs[j];
  }
}

static void encoder_state_remove_refs(encoder_state_t *state) {
  const encoder_control_t * const encoder = state->encoder_control;
  int8_t refnumber = encoder->cfg->ref_frames;
  int8_t check_refs = 0;
  if (encoder->cfg->gop_len) {
    refnumber = encoder->cfg->gop[state->global->gop_offset].ref_neg_count + encoder->cfg->gop[state->global->gop_offset].ref_pos_count;
    check_refs = 1;
  } else if (state->global->slicetype == KVZ_SLICE_I) {
    refnumber = 0;
  }
  // Remove the ref pic (if present)
  while (check_refs || state->global->ref->used_size > (uint32_t)refnumber) {
    int8_t ref_to_remove = state->global->ref->used_size - 1;
    if (encoder->cfg->gop_len) {
      for (int ref = 0; ref < state->global->ref->used_size; ref++) {
        uint8_t found = 0;
        for (int i = 0; i < encoder->cfg->gop[state->global->gop_offset].ref_neg_count; i++) {
          if (state->global->ref->pocs[ref] == state->global->poc - encoder->cfg->gop[state->global->gop_offset].ref_neg[i]) {
            found = 1;
            break;
          }
        }
        if (found) continue;
        for (int i = 0; i < encoder->cfg->gop[state->global->gop_offset].ref_pos_count; i++) {
          if (state->global->ref->pocs[ref] == state->global->poc + encoder->cfg->gop[state->global->gop_offset].ref_pos[i]) {
            found = 1;
            break;
          }
        }
        if (!found) {
          kvz_image_list_rem(state->global->ref, ref);
          ref--;
        }
      }
      check_refs = 0;
    } else kvz_image_list_rem(state->global->ref, ref_to_remove);
  }

}

static void encoder_state_reset_poc(encoder_state_t *state) {
  int i;

  state->global->poc = 0;
  kvz_videoframe_set_poc(state->tile->frame, 0);
  
  for (i=0; state->children[i].encoder_control; ++i) {
    encoder_state_t *sub_state = &(state->children[i]);
    encoder_state_reset_poc(sub_state);
  }
}

static void encoder_state_new_frame(encoder_state_t * const state) {
  int i;
  //FIXME Move this somewhere else!
  if (state->type == ENCODER_STATE_TYPE_MAIN) {
    const encoder_control_t * const encoder = state->encoder_control;

    if (state->global->frame == 0) {
      state->global->is_idr_frame = true;
    }  else if (encoder->cfg->gop_len) {
      // Closed GOP / CRA is not yet supported.
      state->global->is_idr_frame = false;
    
      // Calculate POC according to the global frame counter and GOP structure
      int32_t poc = state->global->frame - 1;
      int32_t poc_offset = encoder->cfg->gop[state->global->gop_offset].poc_offset;
      state->global->poc = poc - poc % encoder->cfg->gop_len + poc_offset;
      kvz_videoframe_set_poc(state->tile->frame, state->global->poc);
    } else {
      bool is_i_idr = (encoder->cfg->intra_period == 1 && state->global->frame % 2 == 0);
      bool is_p_idr = (encoder->cfg->intra_period > 1 && (state->global->frame % encoder->cfg->intra_period) == 0);
      state->global->is_idr_frame = is_i_idr || is_p_idr;
    }
   
    if (state->global->is_idr_frame) {
      encoder_state_reset_poc(state);
      state->global->slicetype = KVZ_SLICE_I;
      state->global->pictype = KVZ_NAL_IDR_W_RADL;
    } else {
      state->global->slicetype = encoder->cfg->intra_period==1 ? KVZ_SLICE_I : (state->encoder_control->cfg->gop_len?KVZ_SLICE_B:KVZ_SLICE_P);

      // Use P-slice for lowdelay.
      if (state->global->slicetype == KVZ_SLICE_B && encoder->cfg->gop_lowdelay) {
        state->global->slicetype = KVZ_SLICE_P;
      }

      state->global->pictype = KVZ_NAL_TRAIL_R;
      if (state->encoder_control->cfg->gop_len) {
        if (encoder->cfg->intra_period > 1 && (state->global->poc % encoder->cfg->intra_period) == 0) {
          state->global->slicetype = KVZ_SLICE_I;
        }
      }

    }

    encoder_state_remove_refs(state);
    encoder_state_ref_sort(state);
    double lambda;
    if (encoder->cfg->target_bitrate > 0) {
      // Rate control enabled.
      lambda = kvz_select_picture_lambda(state);
      state->global->QP = kvz_lambda_to_QP(lambda);
    } else {
      if (encoder->cfg->gop_len > 0 && state->global->slicetype != KVZ_SLICE_I) {
        kvz_gop_config const * const gop =
          encoder->cfg->gop + state->global->gop_offset;
        state->global->QP = encoder->cfg->qp + gop->qp_offset;
        state->global->QP_factor = gop->qp_factor;
      } else {
        state->global->QP = encoder->cfg->qp;
      }
      lambda = kvz_select_picture_lambda_from_qp(state);
    }
    state->global->cur_lambda_cost = lambda;
    state->global->cur_lambda_cost_sqrt = sqrt(lambda);

  }
  kvz_bitstream_clear(&state->stream);
  
  if (state->is_leaf) {
    //Leaf states have cabac and context
    kvz_cabac_start(&state->cabac);
    kvz_init_contexts(state, state->global->QP, state->global->slicetype);
  }
  
  //Clear the jobs
  state->tqj_bitstream_written = NULL;
  state->tqj_recon_done = NULL;
  
  for (i = 0; state->children[i].encoder_control; ++i) {
    encoder_state_new_frame(&state->children[i]);
  }
  

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


void kvz_encode_one_frame(encoder_state_t * const state)
{
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    encoder_state_new_frame(state);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=new_frame,frame=%d,poc=%d", state->global->frame, state->global->poc);
  }
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    encoder_state_encode(state);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, state->encoder_control->threadqueue, "type=encode,frame=%d", state->global->frame);
  }
  //kvz_threadqueue_flush(main_state->encoder_control->threadqueue);
  {
    threadqueue_job_t *job;
#ifdef KVZ_DEBUG
    char job_description[256];
    sprintf(job_description, "type=write_bitstream,frame=%d", state->global->frame);
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
  state->frame_done = 0;
  //kvz_threadqueue_flush(main_state->encoder_control->threadqueue);
}


void kvz_encoder_next_frame(encoder_state_t *state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  // The previous frame must be done before the next one is started.
  assert(state->frame_done);

  if (state->global->frame == -1) {
    //We're at the first frame, so don't care about all this stuff;
    state->global->frame = 0;
    state->global->poc = 0;
    assert(!state->tile->frame->source);
    assert(!state->tile->frame->rec);
    state->tile->frame->rec = kvz_image_alloc(state->tile->frame->width, state->tile->frame->height);
    assert(state->tile->frame->rec);
    state->prepared = 1;
    return;
  }
  
  if (state->previous_encoder_state != state) {
    encoder_state_t *prev_state = state->previous_encoder_state;

    //We have a "real" previous encoder
    state->global->frame = prev_state->global->frame + 1;
    state->global->poc = prev_state->global->poc + 1;

    kvz_cu_array_free(state->tile->frame->cu_array);
    kvz_image_free(state->tile->frame->source);
    state->tile->frame->source = NULL;
    kvz_image_free(state->tile->frame->rec);
    state->tile->frame->rec = kvz_image_alloc(state->tile->frame->width, state->tile->frame->height);
    assert(state->tile->frame->rec);
    {
      // Allocate height_in_scu x width_in_scu x sizeof(CU_info)
      unsigned height_in_scu = state->tile->frame->height_in_lcu << MAX_DEPTH;
      unsigned width_in_scu = state->tile->frame->width_in_lcu << MAX_DEPTH;
      state->tile->frame->cu_array = kvz_cu_array_alloc(width_in_scu, height_in_scu);
    }
    kvz_videoframe_set_poc(state->tile->frame, state->global->poc);
    kvz_image_list_copy_contents(state->global->ref, prev_state->global->ref);
    if (!encoder->cfg->gop_len ||
        !prev_state->global->poc ||
        encoder->cfg->gop[prev_state->global->gop_offset].is_ref) {
      kvz_image_list_add(state->global->ref,
                     prev_state->tile->frame->rec,
                     prev_state->tile->frame->cu_array,
                     prev_state->global->poc);
    }

    state->prepared = 1;
    return;
  }


  if (!encoder->cfg->gop_len ||
      !state->global->poc ||
      encoder->cfg->gop[state->global->gop_offset].is_ref) {
    // Add current reconstructed picture as reference
    kvz_image_list_add(state->global->ref,
                   state->tile->frame->rec,
                   state->tile->frame->cu_array,
                   state->global->poc);
  }


  state->global->frame++;
  state->global->poc++;

  // Remove current source picture.
  kvz_image_free(state->tile->frame->source);
  state->tile->frame->source = NULL;

  // Remove current reconstructed picture, and alloc a new one
  kvz_image_free(state->tile->frame->rec);

  state->tile->frame->rec = kvz_image_alloc(state->tile->frame->width, state->tile->frame->height);
  assert(state->tile->frame->rec);
  kvz_videoframe_set_poc(state->tile->frame, state->global->poc);
  state->prepared = 1;
}

static void encode_part_mode(cabac_data_t * const cabac,
                             const cu_info_t * const cur_cu,
                             int depth)
{
  if (cur_cu->type == CU_INTRA) {
    if (depth == MAX_DEPTH) {
      cabac->cur_ctx = &(cabac->ctx.part_size_model[0]);
      if (cur_cu->part_size == SIZE_2Nx2N) {
        CABAC_BIN(cabac, 1, "part_mode 2Nx2N");
      } else {
        CABAC_BIN(cabac, 0, "part_mode NxN");
      }
    }
  } else {
    // TODO: Handle inter sizes other than 2Nx2N
    cabac->cur_ctx = &(cabac->ctx.part_size_model[0]);
    CABAC_BIN(cabac, 1, "part_mode 2Nx2N");
  }
}

static void encode_inter_prediction_unit(encoder_state_t * const state,
                                         cabac_data_t * const cabac,
                                         const cu_info_t * const cur_cu,
                                         int x_ctb, int y_ctb, int depth)
{
  // Mergeflag
  int16_t num_cand = 0;
  cabac->cur_ctx = &(cabac->ctx.cu_merge_flag_ext_model);
  CABAC_BIN(cabac, cur_cu->merged, "MergeFlag");
  num_cand = MRG_MAX_NUM_CANDS;
  if (cur_cu->merged) { //merge
    if (num_cand > 1) {
      int32_t ui;
      for (ui = 0; ui < num_cand - 1; ui++) {
        int32_t symbol = (ui != cur_cu->merge_idx);
        if (ui == 0) {
          cabac->cur_ctx = &(cabac->ctx.cu_merge_idx_ext_model);
          CABAC_BIN(cabac, symbol, "MergeIndex");
        } else {
          CABAC_BIN_EP(cabac,symbol,"MergeIndex");
        }
        if (symbol == 0) break;
      }
    }
  } else {
    uint32_t ref_list_idx;
    uint32_t j;
    int ref_list[2] = { 0, 0 };
    for (j = 0; j < state->global->ref->used_size; j++) {
      if (state->global->ref->pocs[j] < state->global->poc) {
        ref_list[0]++;
      } else {
        ref_list[1]++;
      }
    }

    // Void TEncSbac::codeInterDir( TComDataCU* pcCU, UInt uiAbsPartIdx )
    if (state->global->slicetype == KVZ_SLICE_B)
    {
      // Code Inter Dir
      uint8_t inter_dir = cur_cu->inter.mv_dir-1;
      uint8_t ctx = depth;
      

      if (cur_cu->part_size == SIZE_2Nx2N || (LCU_WIDTH >> depth) != 8)
      {
        cabac->cur_ctx = &(cabac->ctx.inter_dir[ctx]);
        CABAC_BIN(cabac, (inter_dir == 2), "inter_pred_idc");
      }
      if (inter_dir < 2)
      {
        cabac->cur_ctx = &(cabac->ctx.inter_dir[4]);
        CABAC_BIN(cabac, inter_dir, "inter_pred_idc");
      }
    }

    for (ref_list_idx = 0; ref_list_idx < 2; ref_list_idx++) {
      if (cur_cu->inter.mv_dir & (1 << ref_list_idx)) {
        if (ref_list[ref_list_idx] > 1) {
          // parseRefFrmIdx
          int32_t ref_frame = cur_cu->inter.mv_ref_coded[ref_list_idx];

          cabac->cur_ctx = &(cabac->ctx.cu_ref_pic_model[0]);
          CABAC_BIN(cabac, (ref_frame != 0), "ref_idx_lX");

          if (ref_frame > 0) {
            int32_t i;
            int32_t ref_num = ref_list[ref_list_idx] - 2;

            cabac->cur_ctx = &(cabac->ctx.cu_ref_pic_model[1]);
            ref_frame--;

            for (i = 0; i < ref_num; ++i) {
              const uint32_t symbol = (i == ref_frame) ? 0 : 1;

              if (i == 0) {
                CABAC_BIN(cabac, symbol, "ref_idx_lX");
              } else {
                CABAC_BIN_EP(cabac, symbol, "ref_idx_lX");
              }
              if (symbol == 0) break;
            }
          }
        }

        if (!(/*pcCU->getSlice()->getMvdL1ZeroFlag() &&*/ state->global->ref_list == REF_PIC_LIST_1 && cur_cu->inter.mv_dir == 3)) {
          const int32_t mvd_hor = cur_cu->inter.mvd[ref_list_idx][0];
          const int32_t mvd_ver = cur_cu->inter.mvd[ref_list_idx][1];
          const int8_t hor_abs_gr0 = mvd_hor != 0;
          const int8_t ver_abs_gr0 = mvd_ver != 0;
          const uint32_t mvd_hor_abs = abs(mvd_hor);
          const uint32_t mvd_ver_abs = abs(mvd_ver);

          cabac->cur_ctx = &(cabac->ctx.cu_mvd_model[0]);
          CABAC_BIN(cabac, (mvd_hor != 0), "abs_mvd_greater0_flag_hor");
          CABAC_BIN(cabac, (mvd_ver != 0), "abs_mvd_greater0_flag_ver");

          cabac->cur_ctx = &(cabac->ctx.cu_mvd_model[1]);

          if (hor_abs_gr0) {
            CABAC_BIN(cabac, (mvd_hor_abs>1), "abs_mvd_greater1_flag_hor");
          }

          if (ver_abs_gr0) {
            CABAC_BIN(cabac, (mvd_ver_abs>1), "abs_mvd_greater1_flag_ver");
          }

          if (hor_abs_gr0) {
            if (mvd_hor_abs > 1) {
              kvz_cabac_write_ep_ex_golomb(cabac,mvd_hor_abs-2, 1);
            }

            CABAC_BIN_EP(cabac, (mvd_hor>0)?0:1, "mvd_sign_flag_hor");
          }

          if (ver_abs_gr0) {
            if (mvd_ver_abs > 1) {
              kvz_cabac_write_ep_ex_golomb(cabac,mvd_ver_abs-2, 1);
            }

            CABAC_BIN_EP(cabac, (mvd_ver>0)?0:1, "mvd_sign_flag_ver");
          }
        }

        // Signal which candidate MV to use
        kvz_cabac_write_unary_max_symbol(cabac, cabac->ctx.mvp_idx_model, cur_cu->inter.mv_cand[ref_list_idx], 1,
                                    AMVP_MAX_NUM_CANDS - 1);
      }
    } // for ref_list
  } // if !merge
}

static void encode_intra_coding_unit(encoder_state_t * const state,
                                     cabac_data_t * const cabac,
                                     const cu_info_t * const cur_cu,
                                     int x_ctb, int y_ctb, int depth)
{
  const videoframe_t * const frame = state->tile->frame;
  uint8_t intra_pred_mode[4] = {
    cur_cu->intra[0].mode, cur_cu->intra[1].mode,
    cur_cu->intra[2].mode, cur_cu->intra[3].mode };
    uint8_t intra_pred_mode_chroma = cur_cu->intra[0].mode_chroma;
  int8_t intra_preds[4][3] = {{-1, -1, -1},{-1, -1, -1},{-1, -1, -1},{-1, -1, -1}};
  int8_t mpm_preds[4] = {-1, -1, -1, -1};
  int i, j;
  uint32_t flag[4];
  int num_pred_units = (cur_cu->part_size == SIZE_2Nx2N ? 1 : 4);

  #if ENABLE_PCM == 1
  // Code must start after variable initialization
  kvz_cabac_encode_bin_trm(cabac, 0); // IPCMFlag == 0
  #endif

  // PREDINFO CODING
  // If intra prediction mode is found from the predictors,
  // it can be signaled with two EP's. Otherwise we can send
  // 5 EP bins with the full predmode
  for (j = 0; j < num_pred_units; ++j) {
    static const vector2d_t offset[4] = {{0,0},{1,0},{0,1},{1,1}};
    const cu_info_t *left_cu = NULL;
    const cu_info_t *above_cu = NULL;

    if (x_ctb > 0) {
      left_cu = kvz_videoframe_get_cu_const(frame, x_ctb - 1, y_ctb);
    }
    // Don't take the above CU across the LCU boundary.
    if (y_ctb > 0 && (y_ctb & 7) != 0) {
      above_cu = kvz_videoframe_get_cu_const(frame, x_ctb, y_ctb - 1);
    }

    kvz_intra_get_dir_luma_predictor((x_ctb<<3) + (offset[j].x<<2),
                                 (y_ctb<<3) + (offset[j].y<<2),
                                 intra_preds[j], cur_cu,
                                 left_cu, above_cu);
    for (i = 0; i < 3; i++) {
      if (intra_preds[j][i] == intra_pred_mode[j]) {
        mpm_preds[j] = (int8_t)i;
        break;
      }
    }
    flag[j] = (mpm_preds[j] == -1) ? 0 : 1;
  }

  cabac->cur_ctx = &(cabac->ctx.intra_mode_model);
  for (j = 0; j < num_pred_units; ++j) {
    CABAC_BIN(cabac, flag[j], "prev_intra_luma_pred_flag");
  }

  for (j = 0; j < num_pred_units; ++j) {
    // Signal index of the prediction mode in the prediction list.
    if (flag[j]) {
      CABAC_BIN_EP(cabac, (mpm_preds[j] == 0 ? 0 : 1), "mpm_idx");
      if (mpm_preds[j] != 0) {
        CABAC_BIN_EP(cabac, (mpm_preds[j] == 1 ? 0 : 1), "mpm_idx");
      }
    } else {
      // Signal the actual prediction mode.
      int32_t tmp_pred = intra_pred_mode[j];

      // Sort prediction list from lowest to highest.
      if (intra_preds[j][0] > intra_preds[j][1]) SWAP(intra_preds[j][0], intra_preds[j][1], int8_t);
      if (intra_preds[j][0] > intra_preds[j][2]) SWAP(intra_preds[j][0], intra_preds[j][2], int8_t);
      if (intra_preds[j][1] > intra_preds[j][2]) SWAP(intra_preds[j][1], intra_preds[j][2], int8_t);

      // Reduce the index of the signaled prediction mode according to the
      // prediction list, as it has been already signaled that it's not one
      // of the prediction modes.
      for (i = 2; i >= 0; i--) {
        tmp_pred = (tmp_pred > intra_preds[j][i] ? tmp_pred - 1 : tmp_pred);
      }

      CABAC_BINS_EP(cabac, tmp_pred, 5, "rem_intra_luma_pred_mode");
    }
  }

  {  // start intra chroma pred mode coding
    unsigned pred_mode = 5;
    unsigned chroma_pred_modes[4] = {0, 26, 10, 1};

    if (intra_pred_mode_chroma == intra_pred_mode[0]) {
      pred_mode = 4;
    } else if (intra_pred_mode_chroma == 34) {
      // Angular 34 mode is possible only if intra pred mode is one of the
      // possible chroma pred modes, in which case it is signaled with that
      // duplicate mode.
      for (i = 0; i < 4; ++i) {
        if (intra_pred_mode[0] == chroma_pred_modes[i]) pred_mode = i;
      }
    } else {
      for (i = 0; i < 4; ++i) {
        if (intra_pred_mode_chroma == chroma_pred_modes[i]) pred_mode = i;
      }
    }

    // pred_mode == 5 mean intra_pred_mode_chroma is something that can't
    // be coded.
    assert(pred_mode != 5);

    /**
     * Table 9-35 - Binarization for intra_chroma_pred_mode
     *   intra_chroma_pred_mode  bin_string
     *                        4           0
     *                        0         100
     *                        1         101
     *                        2         110
     *                        3         111
     * Table 9-37 - Assignment of ctxInc to syntax elements with context coded bins
     *   intra_chroma_pred_mode[][] = 0, bypass, bypass
     */
    cabac->cur_ctx = &(cabac->ctx.chroma_pred_model[0]);
    if (pred_mode == 4) {
      CABAC_BIN(cabac, 0, "intra_chroma_pred_mode");
    } else {
      CABAC_BIN(cabac, 1, "intra_chroma_pred_mode");
      CABAC_BINS_EP(cabac, pred_mode, 2, "intra_chroma_pred_mode");
    }
  }  // end intra chroma pred mode coding

  kvz_encode_transform_coeff(state, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
}

void kvz_encode_coding_tree(encoder_state_t * const state,
                        uint16_t x_ctb, uint16_t y_ctb, uint8_t depth)
{
  cabac_data_t * const cabac = &state->cabac;
  const videoframe_t * const frame = state->tile->frame;
  const cu_info_t *cur_cu = kvz_videoframe_get_cu_const(frame, x_ctb, y_ctb);
  uint8_t split_flag = GET_SPLITDATA(cur_cu, depth);
  uint8_t split_model = 0;
  
  //Absolute ctb
  uint16_t abs_x_ctb = x_ctb + (state->tile->lcu_offset_x * LCU_WIDTH) / (LCU_WIDTH >> MAX_DEPTH);
  uint16_t abs_y_ctb = y_ctb + (state->tile->lcu_offset_y * LCU_WIDTH) / (LCU_WIDTH >> MAX_DEPTH);

  // Check for slice border FIXME
  uint8_t border_x = ((state->encoder_control->in.width) < (abs_x_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_y = ((state->encoder_control->in.height) < (abs_y_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_split_x = ((state->encoder_control->in.width)  < ((abs_x_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border_split_y = ((state->encoder_control->in.height) < ((abs_y_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (depth != MAX_DEPTH) {
    // Implisit split flag when on border
    if (!border) {
      // Get left and top block split_flags and if they are present and true, increase model number
      if (x_ctb > 0 && GET_SPLITDATA(kvz_videoframe_get_cu_const(frame, x_ctb - 1, y_ctb), depth) == 1) {
        split_model++;
      }

      if (y_ctb > 0 && GET_SPLITDATA(kvz_videoframe_get_cu_const(frame, x_ctb, y_ctb - 1), depth) == 1) {
        split_model++;
      }

      cabac->cur_ctx = &(cabac->ctx.split_flag_model[split_model]);
      CABAC_BIN(cabac, split_flag, "SplitFlag");
    }

    if (split_flag || border) {
      // Split blocks and remember to change x and y block positions
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      kvz_encode_coding_tree(state, x_ctb, y_ctb, depth + 1); // x,y

      // TODO: fix when other half of the block would not be completely over the border
      if (!border_x || border_split_x) {
        kvz_encode_coding_tree(state, x_ctb + change, y_ctb, depth + 1);
      }
      if (!border_y || border_split_y) {
        kvz_encode_coding_tree(state, x_ctb, y_ctb + change, depth + 1);
      }
      if (!border || (border_split_x && border_split_y)) {
        kvz_encode_coding_tree(state, x_ctb + change, y_ctb + change, depth + 1);
      }
      return;
    }
  }



    // Encode skip flag
  if (state->global->slicetype != KVZ_SLICE_I) {
    int8_t ctx_skip = 0; // uiCtxSkip = aboveskipped + leftskipped;
    int ui;
    int16_t num_cand = MRG_MAX_NUM_CANDS;
    // Get left and top skipped flags and if they are present and true, increase context number
    if (x_ctb > 0 && (kvz_videoframe_get_cu_const(frame, x_ctb - 1, y_ctb))->skipped) {
      ctx_skip++;
    }

    if (y_ctb > 0 && (kvz_videoframe_get_cu_const(frame, x_ctb, y_ctb - 1))->skipped) {
      ctx_skip++;
    }

    cabac->cur_ctx = &(cabac->ctx.cu_skip_flag_model[ctx_skip]);
    CABAC_BIN(cabac, cur_cu->skipped, "SkipFlag");

    // IF SKIP
    if (cur_cu->skipped) {
      if (num_cand > 1) {
        for (ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui != cur_cu->merge_idx);
          if (ui == 0) {
            cabac->cur_ctx = &(cabac->ctx.cu_merge_idx_ext_model);
            CABAC_BIN(cabac, symbol, "MergeIndex");
          } else {
            CABAC_BIN_EP(cabac,symbol,"MergeIndex");
          }
          if (symbol == 0) {
            break;
          }
        }
      }
      return;
    }
  }

  // ENDIF SKIP

  // Prediction mode
  if (state->global->slicetype != KVZ_SLICE_I) {
    cabac->cur_ctx = &(cabac->ctx.cu_pred_mode_model);
    CABAC_BIN(cabac, (cur_cu->type == CU_INTRA), "PredMode");
  }

  // part_mode
  encode_part_mode(cabac, cur_cu, depth);

  if (cur_cu->type == CU_INTER) {
    // FOR each part
    encode_inter_prediction_unit(state, cabac, cur_cu, x_ctb, y_ctb, depth);
    // END for each part

    {
      int cbf = (cbf_is_set(cur_cu->cbf.y, depth) ||
                 cbf_is_set(cur_cu->cbf.u, depth) ||
                 cbf_is_set(cur_cu->cbf.v, depth));

      // Only need to signal coded block flag if not skipped or merged
      // skip = no coded residual, merge = coded residual
      if (!cur_cu->merged) {
        cabac->cur_ctx = &(cabac->ctx.cu_qt_root_cbf_model);
        CABAC_BIN(cabac, cbf, "rqt_root_cbf");
      }
      // Code (possible) coeffs to bitstream

      if (cbf) {
        kvz_encode_transform_coeff(state, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
      }
    }
  } else if (cur_cu->type == CU_INTRA) {
    encode_intra_coding_unit(state, cabac, cur_cu, x_ctb, y_ctb, depth);
  }

    #if ENABLE_PCM == 1
  // Code IPCM block
  if (cur_cu->type == CU_PCM) {
    kvz_cabac_encode_bin_trm(cabac, 1); // IPCMFlag == 1
      kvz_cabac_finish(cabac);
      kvz_bitstream_add_rbsp_trailing_bits(cabac.stream);
    // PCM sample
      {
      unsigned y, x;

      pixel *base_y = &cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH))    + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH))) * encoder->in.width];
      pixel *base_u = &cur_pic->u_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];
      pixel *base_v = &cur_pic->v_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];

      // Luma
      for (y = 0; y < LCU_WIDTH >> depth; y++) {
        for (x = 0; x < LCU_WIDTH >> depth; x++) {
          kvz_bitstream_put(cabac.stream, base_y[x + y * encoder->in.width], 8);
          }
        }

      // Chroma
      if (encoder->in.video_format != FORMAT_400) {
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            kvz_bitstream_put(cabac.stream, base_u[x + y * (encoder->in.width >> 1)], 8);
          }
        }
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            kvz_bitstream_put(cabac.stream, base_v[x + y * (encoder->in.width >> 1)], 8);
          }
        }
      }
    }
    // end PCM sample
      kvz_cabac_start(cabac);
  } // end Code IPCM block
#endif /* END ENABLE_PCM */
  else { /* Should not happend */
    printf("UNHANDLED TYPE!\r\n");
    assert(0);
    exit(1);
  }

   /* end prediction unit */
  /* end coding_unit */
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


static void encode_transform_unit(encoder_state_t * const state,
                                  int x_pu, int y_pu, int depth)
{
  assert(depth >= 1 && depth <= MAX_PU_DEPTH);

  const videoframe_t * const frame = state->tile->frame;
  uint8_t width = LCU_WIDTH >> depth;
  uint8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  int x_cu = x_pu / 2;
  int y_cu = y_pu / 2;
  const cu_info_t *cur_cu = kvz_videoframe_get_cu_const(frame, x_cu, y_cu);

  coeff_t coeff_y[LCU_WIDTH*LCU_WIDTH+1];
  coeff_t coeff_u[LCU_WIDTH*LCU_WIDTH>>2];
  coeff_t coeff_v[LCU_WIDTH*LCU_WIDTH>>2];
  int32_t coeff_stride = frame->width;

  int8_t scan_idx = kvz_get_scan_order(cur_cu->type, cur_cu->intra[PU_INDEX(x_pu, y_pu)].mode, depth);

  int cbf_y = cbf_is_set(cur_cu->cbf.y, depth + PU_INDEX(x_pu, y_pu));

  if (cbf_y) {
    int x = x_pu * (LCU_WIDTH >> MAX_PU_DEPTH);
    int y = y_pu * (LCU_WIDTH >> MAX_PU_DEPTH);
    coeff_t *orig_pos = &frame->coeff_y[x + y * frame->width];
    for (y = 0; y < width; y++) {
      for (x = 0; x < width; x++) {
        coeff_y[x+y*width] = orig_pos[x];
      }
      orig_pos += coeff_stride;
    }
  }

  // CoeffNxN
  // Residual Coding
  if (cbf_y) {
    kvz_encode_coeff_nxn(state, coeff_y, width, 0, scan_idx, cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip);
  }

  if (depth == MAX_DEPTH + 1 && !(x_pu % 2 && y_pu % 2)) {
    // For size 4x4 luma transform the corresponding chroma transforms are
    // also of size 4x4 covering 8x8 luma pixels. The residual is coded
    // in the last transform unit so for the other ones, don't do anything.
    return;
  }

  if (cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth)) {
    int x, y;
    coeff_t *orig_pos_u, *orig_pos_v;

    if (depth <= MAX_DEPTH) {
      x = x_pu * (LCU_WIDTH >> (MAX_PU_DEPTH + 1));
      y = y_pu * (LCU_WIDTH >> (MAX_PU_DEPTH + 1));
    } else {
      // for 4x4 select top left pixel of the CU.
      x = x_cu * (LCU_WIDTH >> (MAX_DEPTH + 1));
      y = y_cu * (LCU_WIDTH >> (MAX_DEPTH + 1));
    }
    orig_pos_u = &frame->coeff_u[x + y * (frame->width >> 1)];
    orig_pos_v = &frame->coeff_v[x + y * (frame->width >> 1)];
    for (y = 0; y < (width_c); y++) {
      for (x = 0; x < (width_c); x++) {
        coeff_u[x+y*(width_c)] = orig_pos_u[x];
        coeff_v[x+y*(width_c)] = orig_pos_v[x];
      }
      orig_pos_u += coeff_stride>>1;
      orig_pos_v += coeff_stride>>1;
    }

    scan_idx = kvz_get_scan_order(cur_cu->type, cur_cu->intra[0].mode_chroma, depth);

    if (cbf_is_set(cur_cu->cbf.u, depth)) {
      kvz_encode_coeff_nxn(state, coeff_u, width_c, 2, scan_idx, 0);
    }

    if (cbf_is_set(cur_cu->cbf.v, depth)) {
      kvz_encode_coeff_nxn(state, coeff_v, width_c, 2, scan_idx, 0);
    }
  }
}

/**
 * \param encoder
 * \param x_pu            Prediction units' x coordinate.
 * \param y_pu            Prediction units' y coordinate.
 * \param depth           Depth from LCU.
 * \param tr_depth        Depth from last CU.
 * \param parent_coeff_u  What was signaled at previous level for cbf_cb.
 * \param parent_coeff_v  What was signlaed at previous level for cbf_cr.
 */
void kvz_encode_transform_coeff(encoder_state_t * const state, int32_t x_pu,int32_t y_pu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v)
{
  cabac_data_t * const cabac = &state->cabac;
  int32_t x_cu = x_pu / 2;
  int32_t y_cu = y_pu / 2;
  const videoframe_t * const frame = state->tile->frame;
  const cu_info_t *cur_cu = kvz_videoframe_get_cu_const(frame, x_cu, y_cu);

  // NxN signifies implicit transform split at the first transform level.
  // There is a similar implicit split for inter, but it is only used when
  // transform hierarchy is not in use.
  int intra_split_flag = (cur_cu->type == CU_INTRA && cur_cu->part_size == SIZE_NxN);

  // The implicit split by intra NxN is not counted towards max_tr_depth.
  int tr_depth_intra = state->encoder_control->tr_depth_intra;
  int max_tr_depth = (cur_cu->type == CU_INTRA ? tr_depth_intra + intra_split_flag : TR_DEPTH_INTER);

  int8_t split = (cur_cu->tr_depth > depth);

  const int cb_flag_y = cbf_is_set(cur_cu->cbf.y, depth + PU_INDEX(x_pu, y_pu));
  const int cb_flag_u = cbf_is_set(cur_cu->cbf.u, depth);
  const int cb_flag_v = cbf_is_set(cur_cu->cbf.v, depth);

  // The split_transform_flag is not signaled when:
  // - transform size is greater than 32 (depth == 0)
  // - transform size is 4 (depth == MAX_PU_DEPTH)
  // - transform depth is max
  // - cu is intra NxN and it's the first split
  if (depth > 0 &&
      depth < MAX_PU_DEPTH &&
      tr_depth < max_tr_depth &&
      !(intra_split_flag && tr_depth == 0))
  {
    cabac->cur_ctx = &(cabac->ctx.trans_subdiv_model[5 - ((kvz_g_convert_to_bit[LCU_WIDTH] + 2) - depth)]);
    CABAC_BIN(cabac, split, "split_transform_flag");
  }

  // Chroma cb flags are not signaled when one of the following:
  // - transform size is 4 (2x2 chroma transform doesn't exist)
  // - they have already been signaled to 0 previously
  // When they are not present they are inferred to be 0, except for size 4
  // when the flags from previous level are used.
  if (depth < MAX_PU_DEPTH) {
    cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_chroma[tr_depth]);
    if (tr_depth == 0 || parent_coeff_u) {
      CABAC_BIN(cabac, cb_flag_u, "cbf_cb");
    }
    if (tr_depth == 0 || parent_coeff_v) {
      CABAC_BIN(cabac, cb_flag_v, "cbf_cr");
    }
  }

  if (split) {
    uint8_t pu_offset = 1 << (MAX_PU_DEPTH - (depth + 1));
    kvz_encode_transform_coeff(state, x_pu, y_pu, depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    kvz_encode_transform_coeff(state, x_pu + pu_offset, y_pu,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    kvz_encode_transform_coeff(state, x_pu, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    kvz_encode_transform_coeff(state, x_pu + pu_offset, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    return;
  }

  // Luma coded block flag is signaled when one of the following:
  // - prediction mode is intra
  // - transform depth > 0
  // - we have chroma coefficients at this level
  // When it is not present, it is inferred to be 1.
  if(cur_cu->type == CU_INTRA || tr_depth > 0 || cb_flag_u || cb_flag_v) {
      cabac->cur_ctx = &(cabac->ctx.qt_cbf_model_luma[!tr_depth]);
      CABAC_BIN(cabac, cb_flag_y, "cbf_luma");
  }

  if (cb_flag_y | cb_flag_u | cb_flag_v) {
    encode_transform_unit(state, x_pu, y_pu, depth);
  }
}

void kvz_encode_coeff_nxn(encoder_state_t * const state, coeff_t *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode, int8_t tr_skip)
{
  const encoder_control_t * const encoder = state->encoder_control;
  cabac_data_t * const cabac = &state->cabac;
  int c1 = 1;
  uint8_t last_coeff_x = 0;
  uint8_t last_coeff_y = 0;
  int32_t i;
  uint32_t sig_coeffgroup_flag[8 * 8] = { 0 };

  int8_t be_valid = encoder->sign_hiding;
  int32_t scan_pos_sig;
  uint32_t go_rice_param = 0;
  uint32_t blk_pos, pos_y, pos_x, sig, ctx_sig;

  // CONSTANTS
  const uint32_t num_blk_side    = width >> TR_MIN_LOG2_SIZE;
  const uint32_t log2_block_size = kvz_g_convert_to_bit[width] + 2;
  const uint32_t *scan           =
    kvz_g_sig_last_scan[scan_mode][log2_block_size - 1];
  const uint32_t *scan_cg = g_sig_last_scan_cg[log2_block_size - 2][scan_mode];

  // Init base contexts according to block type
  cabac_ctx_t *base_coeff_group_ctx = &(cabac->ctx.cu_sig_coeff_group_model[type]);
  cabac_ctx_t *baseCtx           = (type == 0) ? &(cabac->ctx.cu_sig_model_luma[0]) :
                                 &(cabac->ctx.cu_sig_model_chroma[0]);

  // Scan all coeff groups to find out which of them have coeffs.
  // Populate sig_coeffgroup_flag with that info.
  unsigned sig_cg_cnt = 0;
  for (int cg_y = 0; cg_y < width / 4; ++cg_y) {
    for (int cg_x = 0; cg_x < width / 4; ++cg_x) {
      unsigned cg_pos = cg_y * width * 4 + cg_x * 4;
      for (int coeff_row = 0; coeff_row < 4; ++coeff_row) {
        // Load four 16-bit coeffs and see if any of them are non-zero.
        unsigned coeff_pos = cg_pos + coeff_row * width;
        uint64_t four_coeffs = *(uint64_t*)(&coeff[coeff_pos]);
        if (four_coeffs) {
          ++sig_cg_cnt;
          unsigned cg_pos_y = (cg_pos >> log2_block_size) >> TR_MIN_LOG2_SIZE;
          unsigned cg_pos_x = (cg_pos & (width - 1)) >> TR_MIN_LOG2_SIZE;
          sig_coeffgroup_flag[cg_pos_x + cg_pos_y * num_blk_side] = 1;
          break;
        }
      }
    }
  }

  // Rest of the code assumes at least one non-zero coeff.
  assert(sig_cg_cnt > 0);

  // Find the last coeff group by going backwards in scan order.
  unsigned scan_cg_last = num_blk_side * num_blk_side - 1;
  while (!sig_coeffgroup_flag[scan_cg[scan_cg_last]]) {
    --scan_cg_last;
  }

  // Find the last coeff by going backwards in scan order.
  unsigned scan_pos_last = scan_cg_last * 16 + 15;
  while (!coeff[scan[scan_pos_last]]) {
    --scan_pos_last;
  }

  int pos_last = scan[scan_pos_last];

  // transform skip flag
  if(width == 4 && encoder->trskip_enable) {
    cabac->cur_ctx = (type == 0) ? &(cabac->ctx.transform_skip_model_luma) : &(cabac->ctx.transform_skip_model_chroma);
    CABAC_BIN(cabac, tr_skip, "transform_skip_flag");
  }

  last_coeff_x = pos_last & (width - 1);
  last_coeff_y = (uint8_t)(pos_last >> log2_block_size);

  // Code last_coeff_x and last_coeff_y
  kvz_encode_last_significant_xy(state, last_coeff_x, last_coeff_y, width, width,
                             type, scan_mode);

  scan_pos_sig  = scan_pos_last;

  // significant_coeff_flag
  for (i = scan_cg_last; i >= 0; i--) {
    int32_t sub_pos        = i << 4; // LOG2_SCAN_SET_SIZE;
    int32_t abs_coeff[16];
    int32_t cg_blk_pos     = scan_cg[i];
    int32_t cg_pos_y       = cg_blk_pos / num_blk_side;
    int32_t cg_pos_x       = cg_blk_pos - (cg_pos_y * num_blk_side);

    uint32_t coeff_signs   = 0;
    int32_t last_nz_pos_in_cg = -1;
    int32_t first_nz_pos_in_cg = 16;
    int32_t num_non_zero = 0;
    go_rice_param = 0;

    if (scan_pos_sig == scan_pos_last) {
      abs_coeff[0] = abs(coeff[pos_last]);
      coeff_signs  = (coeff[pos_last] < 0);
      num_non_zero = 1;
      last_nz_pos_in_cg  = scan_pos_sig;
      first_nz_pos_in_cg = scan_pos_sig;
      scan_pos_sig--;
    }

    if (i == scan_cg_last || i == 0) {
      sig_coeffgroup_flag[cg_blk_pos] = 1;
    } else {
      uint32_t sig_coeff_group   = (sig_coeffgroup_flag[cg_blk_pos] != 0);
      uint32_t ctx_sig  = kvz_context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                      cg_pos_y, width);
      cabac->cur_ctx = &base_coeff_group_ctx[ctx_sig];
      CABAC_BIN(cabac, sig_coeff_group, "coded_sub_block_flag");
    }

    if (sig_coeffgroup_flag[cg_blk_pos]) {
      int32_t pattern_sig_ctx = kvz_context_calc_pattern_sig_ctx(sig_coeffgroup_flag,
                                                             cg_pos_x, cg_pos_y, width);

      for (; scan_pos_sig >= sub_pos; scan_pos_sig--) {
        blk_pos = scan[scan_pos_sig];
        pos_y   = blk_pos >> log2_block_size;
        pos_x   = blk_pos - (pos_y << log2_block_size);
        sig    = (coeff[blk_pos] != 0) ? 1 : 0;

        if (scan_pos_sig > sub_pos || i == 0 || num_non_zero) {
          ctx_sig  = kvz_context_get_sig_ctx_inc(pattern_sig_ctx, scan_mode, pos_x, pos_y,
                                             log2_block_size, type);
          cabac->cur_ctx = &baseCtx[ctx_sig];
          CABAC_BIN(cabac, sig, "sig_coeff_flag");
        }

        if (sig) {
          abs_coeff[num_non_zero] = abs(coeff[blk_pos]);
          coeff_signs              = 2 * coeff_signs + (coeff[blk_pos] < 0);
          num_non_zero++;

          if (last_nz_pos_in_cg == -1) {
            last_nz_pos_in_cg = scan_pos_sig;
          }

          first_nz_pos_in_cg  = scan_pos_sig;
        }
      }
    } else {
      scan_pos_sig = sub_pos - 1;
    }

    if (num_non_zero > 0) {
      int8_t sign_hidden = (last_nz_pos_in_cg - first_nz_pos_in_cg >=
                            4 /*SBH_THRESHOLD*/) ? 1 : 0;
      uint32_t ctx_set  = (i > 0 && type == 0) ? 2 : 0;
      cabac_ctx_t *base_ctx_mod;
      int32_t num_c1_flag, first_c2_flag_idx, idx, first_coeff2;

      if (c1 == 0) {
        ctx_set++;
      }

      c1 = 1;

      base_ctx_mod     = (type == 0) ? &(cabac->ctx.cu_one_model_luma[4 * ctx_set]) :
                         &(cabac->ctx.cu_one_model_chroma[4 * ctx_set]);
      num_c1_flag      = MIN(num_non_zero, C1FLAG_NUMBER);
      first_c2_flag_idx = -1;

      for (idx = 0; idx < num_c1_flag; idx++) {
        uint32_t symbol = (abs_coeff[idx] > 1) ? 1 : 0;
        cabac->cur_ctx = &base_ctx_mod[c1];
        CABAC_BIN(cabac, symbol, "coeff_abs_level_greater1_flag");

        if (symbol) {
          c1 = 0;

          if (first_c2_flag_idx == -1) {
            first_c2_flag_idx = idx;
          }
        } else if ((c1 < 3) && (c1 > 0)) {
          c1++;
        }
      }

      if (c1 == 0) {
        base_ctx_mod = (type == 0) ? &(cabac->ctx.cu_abs_model_luma[ctx_set]) :
                       &(cabac->ctx.cu_abs_model_chroma[ctx_set]);

        if (first_c2_flag_idx != -1) {
          uint8_t symbol = (abs_coeff[first_c2_flag_idx] > 2) ? 1 : 0;
          cabac->cur_ctx      = &base_ctx_mod[0];
          CABAC_BIN(cabac, symbol, "coeff_abs_level_greater2_flag");
        }
      }

      if (be_valid && sign_hidden) {
        CABAC_BINS_EP(cabac, (coeff_signs >> 1), (num_non_zero - 1), "coeff_sign_flag");
      } else {
        CABAC_BINS_EP(cabac, coeff_signs, num_non_zero, "coeff_sign_flag");
      }

      if (c1 == 0 || num_non_zero > C1FLAG_NUMBER) {
        first_coeff2 = 1;

        for (idx = 0; idx < num_non_zero; idx++) {
          int32_t base_level  = (idx < C1FLAG_NUMBER) ? (2 + first_coeff2) : 1;

          if (abs_coeff[idx] >= base_level) {
            kvz_cabac_write_coeff_remain(cabac, abs_coeff[idx] - base_level, go_rice_param);

            if (abs_coeff[idx] > 3 * (1 << go_rice_param)) {
              go_rice_param = MIN(go_rice_param + 1, 4);
            }
          }

          if (abs_coeff[idx] >= 2) {
            first_coeff2 = 0;
          }
        }
      }
    }
  }
}

/*!
 \brief Encode (X,Y) position of the last significant coefficient
 \param lastpos_x X component of last coefficient
 \param lastpos_y Y component of last coefficient
 \param width  Block width
 \param height Block height
 \param type plane type / luminance or chrominance
 \param scan scan type (diag, hor, ver)

 This method encodes the X and Y component within a block of the last significant coefficient.
*/
void kvz_encode_last_significant_xy(encoder_state_t * const state,
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan)
{
  cabac_data_t * const cabac = &state->cabac;
  uint8_t offset_x  = type?0:((TOBITS(width)*3) + ((TOBITS(width)+1)>>2)),offset_y = offset_x;
  uint8_t shift_x   = type?(TOBITS(width)):((TOBITS(width)+3)>>2), shift_y = shift_x;
  int group_idx_x;
  int group_idx_y;
  int last_x,last_y,i;
  cabac_ctx_t *base_ctx_x = (type ? cabac->ctx.cu_ctx_last_x_chroma : cabac->ctx.cu_ctx_last_x_luma);
  cabac_ctx_t *base_ctx_y = (type ? cabac->ctx.cu_ctx_last_y_chroma : cabac->ctx.cu_ctx_last_y_luma);

  if (scan == SCAN_VER) {
    SWAP( lastpos_x, lastpos_y,uint8_t );
  }

  group_idx_x   = g_group_idx[lastpos_x];
  group_idx_y   = g_group_idx[lastpos_y];

  // Last X binarization
  for (last_x = 0; last_x < group_idx_x ; last_x++) {
    cabac->cur_ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(cabac,1,"last_sig_coeff_x_prefix");
  }

  if (group_idx_x < g_group_idx[width - 1]) {
    cabac->cur_ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(cabac,0,"last_sig_coeff_x_prefix");
  }

  // Last Y binarization
  for (last_y = 0; last_y < group_idx_y ; last_y++) {
    cabac->cur_ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,1,"last_sig_coeff_y_prefix");
  }

  if (group_idx_y < g_group_idx[height - 1]) {
    cabac->cur_ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,0,"last_sig_coeff_y_prefix");
  }

  // Last X
  if (group_idx_x > 3) {
    lastpos_x -= g_min_in_group[group_idx_x];

    for (i = ((group_idx_x - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_x>>i) & 1,"last_sig_coeff_x_suffix");
    }
  }

  // Last Y
  if (group_idx_y > 3) {
    lastpos_y -= g_min_in_group[group_idx_y];

    for (i = ((group_idx_y - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_y>>i) & 1,"last_sig_coeff_y_suffix");
    }
  }

  // end LastSignificantXY
}
