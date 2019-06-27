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

// ***********************************************
    // Modified for SHVC.
int kvz_encoder_state_match_ILR_states_of_children(encoder_state_t *const state)
{
  if (state->ILR_state == NULL) return 0; //State has no ILR_state so children can't have one either

  //Check that matched states have matching type, might be proplematic if not.
  for (int i = 0; i < state->num_ILR_states; i++) {
    assert(state->type == state->ILR_state[i].type);
  }

  int retval = 1;

  for(int i = 0; state->children[i].encoder_control; ++i) {
    //TODO: Add to childrens' ILR state from each ILR state of state (i.e. do switch for num_ILR_states
    //Handle each type of state type
    switch ( state->children[i].type ){
    
    case ENCODER_STATE_TYPE_TILE:
    case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
    {
      //Need to add each row in block scaling source height range as an ILR state
      
      //int range[2]; //Range of blocks needed for scaling
      //kvz_blockScalingSrcHeightRange(range, &state->encoder_control->layer.upscaling, state->children[i].tile->offset_y + state->children[i].lcu_order->position_px.y, state->children[i].lcu_order->size.y);

      ////Map the pixel range to LCU pos
      //range[0] = range[0] / LCU_WIDTH; //First LCU that is needed
      //range[1] = (range[1] + LCU_WIDTH - 1) / LCU_WIDTH; //First LCU that is not needed

      //Just set ILR_state to the first state of children. The correct child is calculatet later anyway
      state->children[i].ILR_state = state->ILR_state->children;//&state->ILR_state->children[range[0]];
      state->children[i].num_ILR_states = 0; //range[1] - range[0];

      //Count children
      for (int j = 0; state->ILR_state->children[j].encoder_control; ++j){
        state->children[i].num_ILR_states++;
      }      

      break;
    }
    
    default: //Assumes 1-to-1 correspondence of states
    {
      state->children[i].ILR_state = &state->ILR_state->children[i];
      state->children[i].num_ILR_states = 1;
      
      break;
    }

    }//END SWITCH

    retval &= kvz_encoder_state_match_ILR_states_of_children(&state->children[i]);
  }

  return retval;
}
// ***********************************************
// ***********************************************
// Modified for SHVC.

/**
* Prepare the encoder state for scalability related stuff.
*
* - Copy ref list info to the ilr frame for TMVP
*/
static void scalability_prepare(encoder_state_t *state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  if (encoder->cfg.ILR_frames > 0 && state->ILR_state != NULL && state->ILR_state->tile->frame->rec != NULL) {
    const encoder_state_t *ILR_state = state->ILR_state;
    // Store current ilr frames list of POCs for use in TMVP derivation
    memcpy(ILR_state->tile->frame->rec->ref_pocs, ILR_state->frame->ref->pocs, sizeof(int32_t) * ILR_state->frame->ref->used_size);
    // Also store image info
    memcpy(ILR_state->tile->frame->rec->picture_info, ILR_state->frame->ref->image_info, sizeof(kvz_picture_info_t) * ILR_state->frame->ref->used_size);

  }
}

//Populate extra rps (rps[num_rps]) if existing rps cannot be used
static void populate_local_rps(encoder_state_t *const state, kvz_rps_config *const rps)
{
  //No point in doing this if # of ref pics is 0
  if( rps->num_negative_pics + rps->num_positive_pics == 0 ){
    return;
  }

  //Populate given rps. Try using existing rps as ref.
  const encoder_control_t* const encoder = state->encoder_control;
  const kvz_config *const cfg = &encoder->cfg;
  const kvz_gop_config *const gop = encoder->cfg.gop_len ? &cfg->gop[state->frame->gop_offset] : NULL;

  //Populate rps from frame ref
  int poc_shift = 0;
  for (int j = 0; j < rps->num_negative_pics; j++) {
    int8_t delta_poc = 0;

    if (gop != NULL) {
      int8_t found = 0;
      do {
        delta_poc = gop->ref_neg[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->image_info[i].layer_id != encoder->layer.layer_id ||
            state->frame->ref->image_info[i].temporal_id > gop->tId) {
            continue;
          }
          if (state->frame->ref->pocs[i] == state->frame->poc - delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == rps->num_negative_pics) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }
    else {
      delta_poc = j > 0 ? -rps->delta_poc[j - 1] + 1 : 1; //If no gop, needs to be -last+1 (prev - delta_poc - 1 used later)
    }

    rps->delta_poc[j] = -delta_poc;
    rps->is_used[j] = !state->frame->is_irap;

    //last_poc = delta_poc;
  }
  //last_poc = 0;
  poc_shift = 0;
  for (int j = 0; j < rps->num_positive_pics; j++) {
    int8_t delta_poc = 0;

    if (gop != NULL) {
      int8_t found = 0;
      do {
        delta_poc = gop->ref_pos[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->image_info[i].layer_id != encoder->layer.layer_id ||
            state->frame->ref->image_info[i].temporal_id > gop->tId) {
            continue;
          }
          if (state->frame->ref->pocs[i] == state->frame->poc + delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == rps->num_positive_pics) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }
    else {
      delta_poc = j > 0 ? rps->delta_poc[j - 1] + 1 : 1; //If no gop needs to be 1+last (delta_poc - prev - 1 used later)
    }

    rps->delta_poc[j + rps->num_negative_pics] = delta_poc;
    rps->is_used[j + rps->num_negative_pics] = !state->frame->is_irap;

  }

  rps->inter_rps_pred_flag = 0;

  //Try to use inter rps pred if num rps > 0
  // With gop try to use the rps that matches the current (gop_offset-1)%gop_len
  // If no gop used, try the prev rps if 
  if (cfg->num_rps > 0) {
    const kvz_rps_config *ref_rps = NULL;
    int delta_rps = 0;
    if (gop != NULL) {
      if (state->frame->gop_offset != 0) {
        ref_rps = &cfg->rps[state->frame->gop_offset - 1];
        rps->delta_ridx = cfg->num_rps - (state->frame->gop_offset - 1);
        delta_rps = (ref_rps + 1)->delta_rps;
      }
      else {
        ref_rps = &cfg->rps[1];
        rps->delta_ridx = cfg->num_rps - 1;
        delta_rps = -ref_rps->delta_rps;
      }
    }
    else {
      ref_rps = &cfg->rps[MAX(rps->num_negative_pics - 1, 0)];
      rps->delta_ridx = 1;
      delta_rps = -1;
    }

    int num_ref_pic = ref_rps->num_negative_pics + ref_rps->num_positive_pics;
    rps->delta_rps = delta_rps;
    rps->num_ref_idc = num_ref_pic + 1;

    //loop through ref pics
    int count = 0;
    for (int j = 0; j < rps->num_ref_idc; j++) {
      int ref_delta_POC = (j < num_ref_pic) ? ref_rps->delta_poc[j] : 0;
      rps->ref_idc[j] = 0;
      rps->is_used[j] = 0;
      if( count != rps->num_negative_pics + rps->num_positive_pics){
        //Loop through pics in the cur rps 
        for (int k = 0; k < rps->num_negative_pics + rps->num_positive_pics; k++) {
          if (rps->delta_poc[k] == ref_delta_POC + delta_rps) {
            //Found a poc that mathces the ref poc
            rps->ref_idc[j] = 1;
            rps->is_used[j] = !state->frame->is_irap;
            count++;
            break;
          }
        }
      }
    }
    //Check that we found a match for all ref idc
    if (count == rps->num_negative_pics + rps->num_positive_pics) {
      rps->inter_rps_pred_flag = 1;
    }
  }
}

/**
* Set is used status for cur ref list
*/
static void set_cur_is_used(encoder_state_t *state){
  memset(state->local_rps->is_used, 0, sizeof(state->local_rps->is_used));
  const int cur_poc = state->frame->poc;
  const kvz_rps_config *const rps = &state->local_rps->rps;
  
  //Get delta poc from rps
  int num_negative_pics = 0;
  int num_positive_pics = 0;
  int num_delta_poc = 0;
  int16_t delta_pocs[MAX_REF_PIC_COUNT] = { 0 };
  uint8_t used_by_curr[MAX_REF_PIC_COUNT] = { 0 };
  //Find the rps is_used ind that poc-wise matches the reference picture
  if (!state->local_rps->rps.inter_rps_pred_flag) {
    //No inter pred used so can use delta poc directly to derive the rps ref poc for the given i
    num_negative_pics = rps->num_negative_pics;
    num_positive_pics = rps->num_positive_pics;
    num_delta_poc = num_negative_pics + num_positive_pics;
    memcpy(delta_pocs, rps->delta_poc, sizeof(delta_pocs));
    memcpy(used_by_curr, rps->is_used, sizeof(used_by_curr));
  } else {
    //Need to derive the rps poc from the reference rps delta pocs [Rec. ITU-T H.265 (7-59)/(7-60)
    const kvz_rps_config *const ref_rps = &state->encoder_control->cfg.rps[state->local_rps->rps_idx - rps->delta_ridx];
    int cur_ind = 0;
    //Add future references that become past references in the cur rps
    for (int i = ref_rps->num_negative_pics + ref_rps->num_positive_pics - 1; i >= ref_rps->num_negative_pics; i--) {
      const int d_poc = ref_rps->delta_poc[i] + rps->delta_rps;
      if (d_poc < 0 && rps->ref_idc[i]) {
        delta_pocs[cur_ind] = d_poc;
        used_by_curr[cur_ind++] = rps->is_used[i];
      }
    }
    //Add ref frame if it is in the past
    if (rps->delta_rps < 0 && rps->ref_idc[rps->num_ref_idc - 1]) {
      delta_pocs[cur_ind] = rps->delta_rps;
      used_by_curr[cur_ind++] = rps->is_used[rps->num_ref_idc - 1];
    }
    //Add other past references that stay as past references
    for (int i = 0; i < ref_rps->num_negative_pics; i++) {
      const int d_poc = ref_rps->delta_poc[i] + rps->delta_rps;
      if (d_poc < 0 && rps->ref_idc[i]) {
        delta_pocs[cur_ind] = d_poc;
        used_by_curr[cur_ind++] = rps->is_used[i];
      }
    }
    num_negative_pics = cur_ind;

    //Add past references that become future references in the cur rps
    for (int i = ref_rps->num_negative_pics - 1; i >= 0; i--) {
      const int d_poc = ref_rps->delta_poc[i] + rps->delta_rps;
      if (d_poc > 0 && rps->ref_idc[i]) {
        delta_pocs[cur_ind] = d_poc;
        used_by_curr[cur_ind++] = rps->is_used[i];
      }
    }
    //Add ref frame if it is in the future
    if (rps->delta_rps > 0 && rps->ref_idc[rps->num_ref_idc - 1]) {
      delta_pocs[cur_ind] = rps->delta_rps;
      used_by_curr[cur_ind++] = rps->is_used[rps->num_ref_idc - 1];
    }
    //Add other future references that stay as future references
    for (int i = ref_rps->num_negative_pics; i < ref_rps->num_positive_pics + ref_rps->num_negative_pics; i++) {
      const int d_poc = ref_rps->delta_poc[i] + rps->delta_rps;
      if (d_poc > 0 && rps->ref_idc[i]) {
        delta_pocs[cur_ind] = d_poc;
        used_by_curr[cur_ind++] = rps->is_used[i];
      }
    }
    num_positive_pics = cur_ind - num_negative_pics;
    num_delta_poc = cur_ind;
  }

  //Loop over reference frames
  state->local_rps->num_ref_idx_LX_active[0] = 0;
  state->local_rps->num_ref_idx_LX_active[1] = 0;
  int count = 0;
  for(int ref_idx = 0; ref_idx < state->frame->ref->used_size; ref_idx++){
    //If ref is an ilr always set it as used
    if(state->frame->ref->image_info[ref_idx].layer_id < state->encoder_control->layer.layer_id ){
      state->local_rps->is_used[ref_idx] = 1;
      state->local_rps->num_ref_idx_LX_active[0]++;
      count++;
      continue;
    }
    //Find the rps is_used ind that poc-wise matches the reference picture
    for (int i = 0; i < num_delta_poc; i++) {
      if( cur_poc + delta_pocs[i] == state->frame->ref->pocs[ref_idx]){
        state->local_rps->is_used[ref_idx] = used_by_curr[i];
        //Update num ridx active
        for (int reflist = 0; reflist < 2; reflist++) {
          for (int j = 0; j < state->frame->ref_LX_size[reflist]; j++) {
            if( state->frame->ref_LX[reflist][j] == ref_idx){
              if(state->local_rps->is_used[ref_idx]) state->local_rps->num_ref_idx_LX_active[reflist]++;
              break;
            }
          }
        }
        count++;
        break;
      }
    }
  }

  assert(count == state->frame->ref->used_size);
}

/**
* Set local rps for state
*/
static void encoder_state_set_rps(encoder_state_t *state)
{
  //I-Slice does not need a rps
  if (state->frame->ref->used_size != 0) {
    //Check if cfg has the correct rps or if we need to generate a new one
    const encoder_control_t *const encoder = state->encoder_control;

    int j;
    int ref_negative = 0;
    int ref_positive = 0;
    if (encoder->cfg.gop_len) {
      for (j = 0; j < state->frame->ref->used_size; j++) {
        if (state->frame->ref->image_info[j].layer_id != encoder->layer.layer_id) {
          continue;
        }
        else if (state->frame->ref->pocs[j] < state->frame->poc) {
          ref_negative++;
        }
        else {
          ref_positive++;
        }
      }
    }
    else ref_negative = state->frame->ref->used_size - encoder->cfg.ILR_frames; //TODO: Need to check actual number of ilr frames?

    uint8_t num_short_term_ref_pic_sets = encoder->cfg.num_rps;//encoder->layer.num_short_term_ref_pic_sets;
    uint8_t short_term_ref_pic_set_sps_flag = num_short_term_ref_pic_sets != 0 ? encoder->layer.short_term_ref_pic_set_sps_flag : 0;

    uint8_t selector = 0;
    selector += short_term_ref_pic_set_sps_flag ? 1 : 0;
    selector += state->frame->ref->used_size - encoder->cfg.ILR_frames == ref_negative + ref_positive ? 2 : 0;
    selector += encoder->cfg.gop_len > 0 ? 4 : 0;
    selector += num_short_term_ref_pic_sets >= 2 ? 8 : 0;
    selector += encoder->cfg.gop_len > 0 && (ref_negative != encoder->cfg.gop[state->frame->gop_offset].ref_neg_count || ref_positive != encoder->cfg.gop[state->frame->gop_offset].ref_pos_count) ? 16 : 0;
    selector += state->frame->is_irap ? 32 : 0;

    state->local_rps->local_st_rps_sps_flag = selector < 16 ? short_term_ref_pic_set_sps_flag : 0;

    switch (selector) {
    case 1:
    case 3:
    case 7:
      state->local_rps->rps_idx = 0;
      break;

    case 9:
    case 11:
      // No gop
      // There should be a valid rps for each ref->used_size - encoder->cfg.ILR_frames
      state->local_rps->rps_idx = MAX(state->frame->ref->used_size - encoder->cfg.ILR_frames - 1, 0);
      break;

    case 15:
      // Gop
      // There should be a valid rps for the cur gop structure indexed by gop offset
      state->local_rps->rps_idx = state->frame->gop_offset;
      break;

    default: {
      // No valid rps
      // Need to write a new rps (with rpsIdx == num_rps)
      // Set num pos/neg ref in the rps
      state->local_rps->rps_idx = encoder->cfg.num_rps;
      kvz_rps_config *rps = &state->local_rps->rps;
      rps->num_negative_pics = ref_negative;
      rps->num_positive_pics = ref_positive;
      populate_local_rps(state, rps);
      set_cur_is_used(state);
      return;
    }
    }//END switch

    //Copy rps to local rps
    state->local_rps->rps = state->encoder_control->cfg.rps[state->local_rps->rps_idx];
    assert(ref_negative == state->local_rps->rps.num_negative_pics);
    assert(ref_positive == state->local_rps->rps.num_positive_pics);
    set_cur_is_used(state);
  } else {
    //Set stuff to zero for good measure
    memset(&state->local_rps->rps, 0, sizeof(kvz_rps_config));
    state->local_rps->rps.num_negative_pics = 0;
    state->local_rps->rps.num_positive_pics = 0;
    state->local_rps->local_st_rps_sps_flag = 0;
    state->local_rps->rps_idx = 0;
    memset(state->local_rps->is_used, 0, sizeof(state->local_rps->is_used));
    memset(state->local_rps->num_ref_idx_LX_active, 0, sizeof(state->local_rps->num_ref_idx_LX_active));
  }
}

//TODO: disable for now, until the need for a waitfor after this function is fixed
////TODO: A better way?
//static void propagate_tqj_ilr_rec_scaling_done_to_children(const encoder_state_t *parent){
//  if (parent->encoder_control != NULL){
//    for (int i = 0; parent->children[i].encoder_control; i++){
//      parent->children[i].tqj_ilr_rec_scaling_done = kvz_threadqueue_copy_ref(parent->tqj_ilr_rec_scaling_done);
//      propagate_tqj_ilr_rec_scaling_done_to_children(&parent->children[i]);
//    }
//  }
//}
//
//
////TODO: A better way?
//static void propagate_tqj_ilr_cua_upsampling_done_to_children(const encoder_state_t *parent){
//  if (parent->encoder_control != NULL){
//    for (int i = 0; parent->children[i].encoder_control; i++){
//      parent->children[i].tqj_ilr_cua_upsampling_done = kvz_threadqueue_copy_ref(parent->tqj_ilr_cua_upsampling_done);
//      propagate_tqj_ilr_cua_upsampling_done_to_children(&parent->children[i]);
//    }
//  }
//}

//TODO: To be depricated
//TODO: Propably overkill, figure out a better way. Need to add bitsream written?
/*static void add_dep_from_children(threadqueue_job_t *job, const encoder_state_t *const state)
{
  if (state->encoder_control != NULL) {
    for (int i = 0; state->children[i].encoder_control; i++) {
      if (state->children[i].tqj_recon_done != NULL) {
        kvz_threadqueue_job_dep_add(job, state->children[i].tqj_recon_done);
      }
      if (state->children[i].tqj_bitstream_written != NULL) {
        kvz_threadqueue_job_dep_add(job, state->children[i].tqj_bitstream_written);
      }
      add_dep_from_children(job, &state->children[i]);
    }
  }
}*/
//TODO: To be depricated
//Start the per wpp scaling jobs recursively
/*static void start_block_scaling_jobs(encoder_state_t *state, kvz_image_scaling_parameter_t *base_param)
{
  //Go deeper until a leaf state is reached
  if (state->is_leaf) {
    switch (state->type) {
    case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
    {
      //Start a job for each ctu on the wavefront row
      for (int i = 0; i < state->lcu_order_count; ++i) {
        const lcu_order_element_t * const lcu = &state->lcu_order[i];

        //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
        kvz_image_scaling_parameter_t *param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
        param->param = base_param->param;
        param->pic_in = kvz_image_copy_ref(base_param->pic_in);
        param->pic_out = kvz_image_copy_ref(base_param->pic_out);
        param->block_x = state->tile->offset_x + lcu->position_px.x;
        param->block_y = state->tile->offset_y + lcu->position_px.y;
        param->block_width = lcu->size.x;
        param->block_height = lcu->size.y;

        kvz_threadqueue_free_job(&state->layer->image_ver_scaling_jobs[lcu->id]);
        state->layer->image_ver_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_block_scaler_worker, (void*)param);

        //Calculate vertical range of block scaling
        int range[2]; //Range of blocks needed for scaling
        kvz_blockScalingSrcWidthRange(range, base_param->param, param->block_x, param->block_width);

        //Map the pixel range to LCU pos
        range[0] = range[0] / LCU_WIDTH; //First LCU that is needed
        range[1] = (range[1] + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed

                                                           //Add dependencies to ilr states
        for (int j = 0; j < state->num_ILR_states; j++) {
          for (int k = range[0]; k < range[1]; k++) {
            const lcu_order_element_t * const ilr_lcu = &state->ILR_state[j].lcu_order[k];
            kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], state->ILR_state[j].tile->wf_jobs[ilr_lcu->id]);
          }
        }

        //Dependencies added so submit the job
        kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_ver_scaling_jobs[lcu->id]);
      }
      break;
    }

    default:
    {
      //TODO: Handle other types?
      break;
    }
    }//END switch
  }
  else {
    for (int i = 0; state->children[i].encoder_control; ++i) {
      start_block_scaling_jobs(&state->children[i], base_param);
    }
  }

}*/
/*
//TODO: To be depricated
// Scale image by adding a scaling worker job to the job queue
static kvz_picture* deferred_block_scaling(kvz_picture* const pic_in, const scaling_parameter_t *const param, encoder_state_t *state, uint8_t skip_same)
{
  if (pic_in == NULL) {
    return NULL;
  }

  //If no scaling needs to be done, just return pic_in
  if (skip_same && param->src_height == param->trgt_height && param->src_width == param->trgt_width) {
    //scaling_param->skip = 1;
    //scaling_param->pic_out = NULL;
    //kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done);
    state->tqj_ilr_rec_scaling_done = NULL;
    //pic_out = kvz_image_copy_ref(pic_in);
    return kvz_image_copy_ref(pic_in);
  }

  //Prepare img_job_param.
  state->layer->img_job_param.param = param;
  kvz_image_free(state->layer->img_job_param.pic_in); //Free prev pic
  state->layer->img_job_param.pic_in = kvz_image_copy_ref(pic_in);

  kvz_image_free(state->layer->img_job_param.pic_out); //Free prev pic
  state->layer->img_job_param.pic_out = kvz_image_alloc(pic_in->chroma_format,
    param->trgt_width + param->trgt_padding_x,
    param->trgt_height + param->trgt_padding_y);

  //Start jobs
  start_block_scaling_jobs(state, &state->layer->img_job_param);

  return kvz_image_copy_ref(state->layer->img_job_param.pic_out);
}*/

//TODO: To be depricated
//Start the per wpp scaling jobs recursively
//TODO: Add support for downsampling?
/*static void start_block_step_scaling_jobs(encoder_state_t * const state, const kvz_image_scaling_parameter_t * const base_param)
{
  //Go deeper until a leaf state is reached
  if (state->is_leaf) {
    const kvz_image_scaling_parameter_t * const state_param = &state->layer->img_job_param;
    switch (state->type) {
    case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
    {
      //Start a job for each ctu on the wavefront row
      for (int i = 0; i < state->lcu_order_count; ++i) {
        const lcu_order_element_t * const lcu = &state->lcu_order[i];

        //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
        kvz_image_scaling_parameter_t *param_ver = calloc(1, sizeof(kvz_image_scaling_parameter_t));
        param_ver->param = base_param->param;
        param_ver->pic_in = NULL;
        param_ver->pic_out = kvz_image_copy_ref(base_param->pic_out);
        param_ver->src_buffer = state_param->src_buffer;
        param_ver->ver_tmp_buffer = state_param->ver_tmp_buffer;
        param_ver->trgt_buffer = state_param->trgt_buffer;
        param_ver->block_x = state->tile->offset_x + lcu->position_px.x;
        param_ver->block_y = state->tile->offset_y + lcu->position_px.y;
        param_ver->block_width = lcu->size.x;
        param_ver->block_height = lcu->size.y;

        kvz_image_scaling_parameter_t *param_hor = calloc(1, sizeof(kvz_image_scaling_parameter_t));
        memcpy(param_hor, param_ver, sizeof(kvz_image_scaling_parameter_t));
        param_hor->pic_in = kvz_image_copy_ref(base_param->pic_in);
        param_hor->pic_out = NULL;//kvz_image_copy_ref(base_param->pic_out);

                                  //First create job for vertical scaling
        kvz_threadqueue_free_job(&state->layer->image_ver_scaling_jobs[lcu->id]);
        state->layer->image_ver_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param_ver);

        //Add dependencies
        //  Create horizontal scaling jobs that the ver job depends on
        int ver_range[2];
        kvz_blockScalingSrcHeightRange(ver_range, base_param->param, param_ver->block_y, param_ver->block_height);
        int set_job_row = 0; //row of the first job not yet set

                             //Check previous block range to avoid re-creating jobs
        if (lcu->above != NULL) {
          int tmp_block_y = lcu->above->encoder_state->tile->offset_y + lcu->above->position_px.y;
          int tmp_block_height = lcu->above->size.y;
          int tmp_range[2];
          kvz_blockScalingSrcHeightRange(tmp_range, base_param->param, tmp_block_y, tmp_block_height);
          set_job_row = (tmp_range[1] + LCU_WIDTH - 1) / LCU_WIDTH;
        }

        //Map the pixel range to LCU row
        ver_range[0] = ver_range[0] / LCU_WIDTH; //First LCU that is needed
        ver_range[1] = (ver_range[1] + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed

                                                                   //Create hor job for each row that does not have one yet and add dep
                                                                   //int id_offset = (lcu->id % state->lcu_order_count) * state->ILR_state->encoder_control->in.height_in_lcu; //Offset the hor job index by the column number of the current lcu
        int id_offset = lcu->index; //times state->ILR_state->encoder_control->in.height_in_lcu; //Offset the hor job index by the column number of the current lcu
        for (int row = ver_range[0]; row < ver_range[1]; row++) {
          int hor_ind = row * state->encoder_control->in.width_in_lcu + id_offset;

          if (row >= set_job_row) {
            //Set correct block parameters for hor job since it may be on a different lcu row than lcu
            const encoder_state_t * const ilr_state = &state->ILR_state[row];
            param_hor->block_y = ilr_state->tile->offset_y + ilr_state->lcu_order->position_px.y;
            param_hor->block_height = ilr_state->lcu_order->size.y;

            kvz_threadqueue_free_job(&state->layer->image_hor_scaling_jobs[hor_ind]);
            state->layer->image_hor_scaling_jobs[hor_ind] = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param_hor);

            //Add dependency to ILR jobs
            //Calculate vertical range of block scaling
            int range_hor[2]; //Range of blocks needed for scaling
            kvz_blockScalingSrcWidthRange(range_hor, base_param->param, param_hor->block_x, param_hor->block_width);

            //Map the pixel range to LCU pos
            range_hor[0] = range_hor[0] / LCU_WIDTH; //First LCU that is needed
            range_hor[1] = (range_hor[1] + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed

                                                                       //Add ilr state dependencies to hor job
            for (int k = range_hor[0]; k < range_hor[1]; k++) {
              kvz_threadqueue_job_dep_add(state->layer->image_hor_scaling_jobs[hor_ind], ilr_state->tile->wf_jobs[ilr_state->lcu_order[k].id]);
            }

            //Add dependency to left lcu so that copying to src_buffer is not an issue
            if (lcu->left != NULL) {
              kvz_threadqueue_job_dep_add(state->layer->image_hor_scaling_jobs[hor_ind], state->layer->image_hor_scaling_jobs[hor_ind - 1]);
            }

            //Submit hor job
            kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_hor_scaling_jobs[hor_ind]);
          }

          //Add hor step dependency to ver step
          kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], state->layer->image_hor_scaling_jobs[hor_ind]);
        }

        //Dependencies added so submit the ver job
        kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_ver_scaling_jobs[lcu->id]);
      }

      //Set last lcu scaling job as rec done job
      state->tqj_ilr_rec_scaling_done = kvz_threadqueue_copy_ref(state->layer->image_hor_scaling_jobs[state->lcu_order[state->lcu_order_count - 1].id]);

      break;
    }

    case ENCODER_STATE_TYPE_TILE:
    {
      //Make a scaling job for each tile

      //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
      kvz_image_scaling_parameter_t *param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
      param->param = base_param->param;
      param->pic_in = kvz_image_copy_ref(base_param->pic_in);;
      param->pic_out = kvz_image_copy_ref(base_param->pic_out);
      param->src_buffer = state_param->src_buffer;
      param->ver_tmp_buffer = state_param->ver_tmp_buffer;
      param->trgt_buffer = state_param->trgt_buffer;
      param->block_x = state->tile->offset_x;
      param->block_y = state->tile->offset_y;
      param->block_width = state_param->trgt_buffer->y->width; //Trgt buffer should be the size of the block
      param->block_height = state_param->trgt_buffer->y->height;

      //Do hor/ver scaling in the same job
      kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done);
      state->tqj_ilr_rec_scaling_done = kvz_threadqueue_job_create(kvz_tile_step_scaler_worker, (void*)param);

      //Need to add dependency to all ilr tiles that that are within the src range
      int range[4];
      kvz_blockScalingSrcWidthRange(range, param->param, param->block_x, param->block_width);
      kvz_blockScalingSrcHeightRange(range + 2, param->param, param->block_y, param->block_height);

      for (int i = 0; i < state->num_ILR_states; i++) {
        const encoder_state_t *ilr_state = &state->ILR_state[i];
        int ilr_tile_x = ilr_state->tile->offset_x;
        int ilr_tile_y = ilr_state->tile->offset_y;
        int ilr_tile_width = ilr_state->tile->frame->width;
        int ilr_tile_height = ilr_state->tile->frame->height;
        int block_x = range[0];
        int block_y = range[2];
        int block_width = range[1] - range[0] + 1;
        int block_height = range[3] - range[2] + 1;

        if (ilr_state->tqj_recon_done != NULL
          && ABS((block_x << 1) + block_width - (ilr_tile_x << 1) - ilr_tile_width) < block_width + ilr_tile_width
          && ABS((block_y << 1) + block_height - (ilr_tile_y << 1) - ilr_tile_height) < block_height + ilr_tile_height) {
          //Areas intersect so add a dependency
          kvz_threadqueue_job_dep_add(state->tqj_ilr_rec_scaling_done, ilr_state->tqj_recon_done);
        }
      }

      //Submit job
      kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tqj_ilr_rec_scaling_done);

      break;
    }

    default:
    {
      //TODO: Handle other types?
      break;
    }
    }//END switch
  }
  else {
    for (int i = 0; state->children[i].encoder_control; ++i) {
      start_block_step_scaling_jobs(&state->children[i], base_param);
    }
  }

}*/

//TODO: To be depricated
// Scale image by adding scaling worker jobs to the job queue
/*static kvz_picture* deferred_block_step_scaling(kvz_picture* const pic_in, const scaling_parameter_t *const param, encoder_state_t *state, uint8_t skip_same)
{
  if (pic_in == NULL) {
    return NULL;
  }

  //If no scaling needs to be done, just return pic_in
  if (skip_same && param->src_height == param->trgt_height && param->src_width == param->trgt_width) {
    return kvz_image_copy_ref(pic_in);
  }

  //Prepare img_job_param.
  state->layer->img_job_param.param = param;
  kvz_image_free(state->layer->img_job_param.pic_in); //Free prev pic
  state->layer->img_job_param.pic_in = kvz_image_copy_ref(pic_in);

  kvz_image_free(state->layer->img_job_param.pic_out); //Free prev pic
  state->layer->img_job_param.pic_out = kvz_image_alloc(pic_in->chroma_format,
    param->trgt_width + param->trgt_padding_x,
    param->trgt_height + param->trgt_padding_y);

  //Copy other information
  state->layer->img_job_param.pic_out->dts = pic_in->dts;
  state->layer->img_job_param.pic_out->pts = pic_in->pts;
  state->layer->img_job_param.pic_out->interlacing = pic_in->interlacing;

  //Start jobs
  start_block_step_scaling_jobs(state, &state->layer->img_job_param);

  return kvz_image_copy_ref(state->layer->img_job_param.pic_out);
}*/
/*
//TODO: To be depricated
// Scale image by adding a scaling worker job to the job queue
static kvz_picture* deferred_image_scaling(kvz_picture* const pic_in, const scaling_parameter_t *const param, encoder_state_t *state, uint8_t skip_same)
{
  if (pic_in == NULL) {
    return NULL;
  }

  //If no scaling needs to be done, just return pic_in
  if (skip_same && param->src_height == param->trgt_height && param->src_width == param->trgt_width) {
    //scaling_param->skip = 1;
    //scaling_param->pic_out = NULL;
    kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done);
    state->tqj_ilr_rec_scaling_done = NULL;
    //pic_out = kvz_image_copy_ref(pic_in);
    return kvz_image_copy_ref(pic_in);
  }

  //Allocate scaling parameters to give to the worker. Worker should handle freeing.
  kvz_image_scaling_parameter_t *scaling_param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
  kvz_picture* pic_out = NULL;

  pic_out = kvz_image_alloc(pic_in->chroma_format,
    param->trgt_width + param->trgt_padding_x,
    param->trgt_height + param->trgt_padding_y);

  scaling_param->pic_out = kvz_image_copy_ref(pic_out);
  scaling_param->pic_in = kvz_image_copy_ref(pic_in);
  scaling_param->param = param;

  //Make new job and free previous
  kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done);
  state->tqj_ilr_rec_scaling_done = kvz_threadqueue_job_create(kvz_image_scaler_worker, scaling_param);

  //Figure out dependency. ILR recon needs to be completed before scaling can be done.
  add_dep_from_children(state->tqj_ilr_rec_scaling_done, state->ILR_state);

  //Submit job and set it to encoder state
  kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tqj_ilr_rec_scaling_done);

  //Propagate tqj_ilr_rec_scaling_done to child states in order to set it as a dependency
  //TODO: disable for now, until the need for a waitfor after this function is fixed
  //propagate_tqj_ilr_rec_scaling_done_to_children(state);

  return pic_out;
}*/

//TODO: To be depricated
// Scale cu_array by adding a scaling worker job to the job queue
/*static cu_array_t* deferred_cu_array_upsampling(encoder_state_t *state, int32_t mv_scale[2], int32_t pos_scale[2], uint8_t skip_same)
{
  if (skip_same && state->tile->frame->width_in_lcu == state->ILR_state->tile->frame->width_in_lcu &&
    state->tile->frame->height_in_lcu == state->ILR_state->tile->frame->height_in_lcu) {
    kvz_threadqueue_free_job(&state->tqj_ilr_cua_upsampling_done);
    state->tqj_ilr_cua_upsampling_done = NULL;
    return kvz_cu_array_copy_ref(state->ILR_state->tile->frame->cu_array);
  }

  //Allocate the new cua.
  uint32_t n_width = state->tile->frame->width_in_lcu * LCU_WIDTH;
  uint32_t n_height = state->tile->frame->height_in_lcu * LCU_WIDTH;
  cu_array_t *cua = kvz_cu_array_alloc(n_width, n_height);

  //Allocate scaling parameters to give to the worker. Worker should handle freeing.
  kvz_cua_upsampling_parameter_t *param = calloc(1, sizeof(kvz_cua_upsampling_parameter_t));
  param->base_cua = state->ILR_state->tile->frame->cu_array;
  param->cu_pos_scale[0] = pos_scale[0];
  param->cu_pos_scale[1] = pos_scale[1];
  param->mv_scale[0] = mv_scale[0];
  param->mv_scale[1] = mv_scale[1];
  param->nh_in_lcu = state->tile->frame->height_in_lcu;
  param->nw_in_lcu = state->tile->frame->width_in_lcu;
  param->only_init = !state->encoder_control->cfg.tmvp_enable;
  param->out_cua = cua;
  param->lcu_ind = -1;

  //Make new job and free previous
  kvz_threadqueue_free_job(&state->tqj_ilr_cua_upsampling_done);
  state->tqj_ilr_cua_upsampling_done = kvz_threadqueue_job_create(kvz_cu_array_upsampling_worker, param);

  //Figure out dependency. ILR recon needs to be completed before scaling can be done.
  add_dep_from_children(state->tqj_ilr_cua_upsampling_done, state->ILR_state);

  //Submit job and set it to encoder state
  kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tqj_ilr_cua_upsampling_done);

  //Propagate tqj_ilr_cua_upsampling_done to child states in order to set it as a dependency
  //TODO: disable for now, until the need for a waitfor after this function is fixed
  //propagate_tqj_ilr_cua_upsampling_done_to_children(state);

  return cua;
}*/

//TODO: To be depricated
//Start the per wpp scaling jobs recursively
/*static void start_cua_lcu_scaling_jobs(encoder_state_t *state, kvz_cua_upsampling_parameter_t *base_param)
{
  //Go deeper until a leaf state is reached
  if (state->is_leaf) {
    switch (state->type) {
    case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
    {
      //Start a job for each ctu on the wavefront row
      for (int i = 0; i < state->lcu_order_count; ++i) {
        const lcu_order_element_t * const lcu = &state->lcu_order[i];

        //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
        kvz_cua_upsampling_parameter_t *param = calloc(1, sizeof(kvz_cua_upsampling_parameter_t));
        param->base_cua = kvz_cu_array_copy_ref(base_param->base_cua);
        param->out_cua = kvz_cu_array_copy_ref(base_param->out_cua);
        param->nw_in_lcu = base_param->nw_in_lcu;
        param->nh_in_lcu = base_param->nh_in_lcu;
        param->cu_pos_scale[0] = base_param->cu_pos_scale[0];
        param->cu_pos_scale[1] = base_param->cu_pos_scale[1];
        param->mv_scale[0] = base_param->mv_scale[0];
        param->mv_scale[1] = base_param->mv_scale[1];
        param->only_init = base_param->only_init;

        int block_x = state->tile->offset_x + lcu->position_px.x;
        int block_y = state->tile->offset_y + lcu->position_px.y;
        int block_width = lcu->size.x;
        int block_height = lcu->size.y;
        int src_width = base_param->base_cua->width;
        int src_height = base_param->base_cua->height;

        param->lcu_ind = lcu->id;

        kvz_threadqueue_free_job(&state->layer->cua_scaling_jobs[lcu->id]);
        state->layer->cua_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_cu_array_upsampling_worker, (void*)param);

        //Calculate (vertical/horizontal) range of scaling
        int range[4]; //Range of blocks needed for scaling
        kvz_cu_array_upsampling_src_range(range, block_x, block_x + block_width - 1, src_width, param->cu_pos_scale[0]); //Width  
        kvz_cu_array_upsampling_src_range(range + 2, block_y, block_y + block_height - 1, src_height, param->cu_pos_scale[1]); //Height

                                                                                                                           //Map the pixel range to LCU pos
        range[0] = range[0] / LCU_WIDTH; //First LCU that is needed
        range[1] = (range[1] + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
        range[2] = range[2] / LCU_WIDTH;
        range[3] = (range[3] + LCU_WIDTH - 1) / LCU_WIDTH;

        //Add dependencies to ilr states
        for (int j = range[2]; j < range[3]; j++) {
          const encoder_state_t * const ilr_state = &state->ILR_state[j];
          for (int k = range[0]; k < range[1]; k++) {
            kvz_threadqueue_job_dep_add(state->layer->cua_scaling_jobs[lcu->id], ilr_state->tile->wf_jobs[ilr_state->lcu_order[k].id]);
          }
        }

        //Dependencies added so submit the job
        kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->cua_scaling_jobs[lcu->id]);
      }
      break;
    }

    default:
    {
      //TODO: Handle other types?
      break;
    }
    }//END switch
  }
  else {
    for (int i = 0; state->children[i].encoder_control; ++i) {
      start_cua_lcu_scaling_jobs(&state->children[i], base_param);
    }
  }

}*/

//TODO: To be depricated
// Scale cu_array by adding a scaling worker job for each lcu
/*static cu_array_t* deferred_cua_lcu_upsampling(encoder_state_t *state, int32_t mv_scale[2], int32_t pos_scale[2], uint8_t skip_same)
{
  if (skip_same && state->tile->frame->width_in_lcu == state->ILR_state->tile->frame->width_in_lcu &&
    state->tile->frame->height_in_lcu == state->ILR_state->tile->frame->height_in_lcu) {

    return kvz_cu_array_copy_ref(state->ILR_state->tile->frame->cu_array);
  }

  //Allocate the new cua.
  uint32_t n_width = state->tile->frame->width_in_lcu * LCU_WIDTH;
  uint32_t n_height = state->tile->frame->height_in_lcu * LCU_WIDTH;
  cu_array_t *cua = kvz_cu_array_alloc(n_width, n_height);

  //Allocate scaling parameters to give to the worker. Worker should handle freeing.
  kvz_cua_upsampling_parameter_t *param = &state->layer->cua_job_param;
  kvz_cu_array_free(&param->base_cua); //Free prev
  param->base_cua = kvz_cu_array_copy_ref(state->ILR_state->tile->frame->cu_array);
  param->cu_pos_scale[0] = pos_scale[0];
  param->cu_pos_scale[1] = pos_scale[1];
  param->mv_scale[0] = mv_scale[0];
  param->mv_scale[1] = mv_scale[1];
  param->nh_in_lcu = state->tile->frame->height_in_lcu;
  param->nw_in_lcu = state->tile->frame->width_in_lcu;
  param->only_init = !state->encoder_control->cfg.tmvp_enable;
  kvz_cu_array_free(&param->out_cua); //Free prev
  param->out_cua = kvz_cu_array_copy_ref(cua);
  param->lcu_ind = -1; //Set correct ind later

                       //Make new job and free previous
  start_cua_lcu_scaling_jobs(state, param);

  return cua;
}*/

/**
* Handle adding ilr frames to the ref list.
*
* - Add ilr frame to the ref list if ilr is enabled
*//*
//TODO: To be depricated
static void add_ilr_frames(encoder_state_t *state) 
{
  const encoder_control_t * const encoder = state->encoder_control;
  // For SHVC.
  //TODO: Account for adding several ILR frames. Should ilr rec ever bee NULL?
  if (encoder->cfg.ILR_frames > 0 && state->ILR_state != NULL && state->ILR_state->tile->frame->rec != NULL) {
    //Also add base layer to the reference list.
    const encoder_state_t *ILR_state = state->ILR_state;
    kvz_picture *ilr_rec = kvz_image_copy_ref(ILR_state->tile->frame->rec);
    kvz_picture *scaled_pic = NULL;
    if (encoder->cfg.threads > 0) {
      //TODO: fix dependencies etc. so that waitfor does not need to be called here
      scaled_pic = deferred_block_step_scaling(ilr_rec, &encoder->layer.upscaling, state, 1);//deferred_image_scaling(ilr_rec, &encoder->layer.upscaling, state, 1);
                                                                                             //if (state->tqj_ilr_rec_scaling_done != NULL) {
                                                                                             //  kvz_threadqueue_waitfor(state->encoder_control->threadqueue, state->tqj_ilr_rec_scaling_done);
                                                                                             //}
    }
    else {
      scaled_pic = kvz_image_scaling(ilr_rec, &encoder->layer.upscaling, 1);
    }
    if (ilr_rec == NULL || scaled_pic == NULL) {
      return; //TODO: Add error etc?
    }
    //Copy image ref info
    //TODO: Is there a way to determine how many refs irl_rec has? Otherwise just copy everything to be safe
    memcpy(scaled_pic->ref_pocs, ilr_rec->ref_pocs, sizeof(ilr_rec->ref_pocs));//sizeof(int32_t) * ILR_state->frame->ref->used_size);
    memcpy(scaled_pic->picture_info, ilr_rec->picture_info, sizeof(ilr_rec->picture_info));//sizeof(kvz_picture_info_t) * ILR_state->frame->ref->used_size);
    kvz_image_free(ilr_rec);

    //TODO: Account for offsets etc. Need to use something else than original sizes?
    //TODO: Shm doesn't upsample the cua if frame is idr; need to do something else as well when idr? Does it even matter?
    int32_t mv_scale[2] = { GET_SCALE_MV(encoder->layer.upscaling.src_width,encoder->layer.upscaling.trgt_width),
      GET_SCALE_MV(encoder->layer.upscaling.src_height,encoder->layer.upscaling.trgt_height) };
    int32_t pos_scale[2] = { GET_SCALE_POS(encoder->layer.upscaling.src_width,encoder->layer.upscaling.trgt_width),
      GET_SCALE_POS(encoder->layer.upscaling.src_height,encoder->layer.upscaling.trgt_height) };
    cu_array_t* scaled_cu = NULL;

    if (encoder->cfg.threads > 0) {
      //TODO: fix dependencies etc. so that waitfor does not need to be called here
      scaled_cu = deferred_cua_lcu_upsampling(state, mv_scale, pos_scale, 1);//deferred_cu_array_upsampling( state, mv_scale, pos_scale, 1);
                                                                             //scaled_cu = kvz_cu_array_copy_ref(ILR_state->tile->frame->cu_array);
      if (state->tqj_ilr_cua_upsampling_done != NULL) {
        kvz_threadqueue_waitfor(state->encoder_control->threadqueue, state->tqj_ilr_cua_upsampling_done);
      }
    }
    else {
      scaled_cu = kvz_cu_array_upsampling(ILR_state->tile->frame->cu_array,
        state->tile->frame->width_in_lcu,
        state->tile->frame->height_in_lcu,
        mv_scale, pos_scale, 1,
        !state->encoder_control->cfg.tmvp_enable);//(state->frame->pictype >= KVZ_NAL_BLA_W_LP && state->frame->pictype <= KVZ_NAL_CRA_NUT)); //pic type not set yet 
    }

    kvz_image_list_add(state->frame->ref,
      scaled_pic,
      scaled_cu,
      ILR_state->frame->poc,
      ILR_state->frame->ref_LX,
      ILR_state->encoder_control->cfg.gop[ILR_state->frame->gop_offset].tId,
      ILR_state->encoder_control->layer.layer_id,
      1); //Currently only ILR can be a long term references
          //TODO: Add error handling?
    kvz_image_free(scaled_pic);
    kvz_cu_array_free(&scaled_cu);
  }
}*/

// Scale tile specified area
static void block_step_scaling(encoder_state_t * const state )
{
  //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
  kvz_image_scaling_parameter_t * const param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
  kvz_copy_image_scaling_parameters(param, &state->layer->img_job_param);
  param->block_x = state->tile->offset_x;
  param->block_y = state->tile->offset_y;
  param->block_width = state->layer->img_job_param.trgt_buffer->y->width; //Trgt buffer should be the size of the block
  param->block_height = state->layer->img_job_param.trgt_buffer->y->height;

  param->use_tiles = 1;
  kvz_block_step_scaler_worker(param);

  //state->layer->scaling_started = 1;

}

// Start the block step scaling job(s) for the given state
static void start_block_step_scaling_job(encoder_state_t * const state, const lcu_order_element_t * const lcu)
{
  //If scaling has already been started, no need to do it here anymore
  //if( state->layer == NULL || state->layer->scaling_started ){
  //  return;
  //}

  const kvz_image_scaling_parameter_t * const state_param = &state->layer->img_job_param;
  switch (state->type) {
  case ENCODER_STATE_TYPE_WAVEFRONT_ROW: 
  {
#define SINGLE_STAGE_SCALING
#ifndef SINGLE_STAGE_SCALING
    //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
    kvz_image_scaling_parameter_t *param_ver = calloc(1, sizeof(kvz_image_scaling_parameter_t));
    kvz_copy_image_scaling_parameters(param_ver, state_param);
    param_ver->block_x = state->tile->offset_x + lcu->position_px.x;
    param_ver->block_y = state->tile->offset_y + lcu->position_px.y;
    param_ver->block_width = lcu->size.x;
    param_ver->block_height = lcu->size.y;

    kvz_image_free(param_ver->pic_in);
    param_ver->pic_in = NULL;
    
    param_ver->use_tiles = 0;

    //First create job for vertical scaling
    kvz_threadqueue_free_job(&state->layer->image_ver_scaling_jobs[lcu->id]);
    state->layer->image_ver_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param_ver);

    //Add dependencies
    //  Create horizontal scaling jobs that the ver job depends on
    int ver_range[2];
    kvz_blockScalingSrcHeightRange(ver_range, state_param->param, param_ver->block_y, param_ver->block_height);
    int set_job_row = 0; //row of the first job not yet set

    //Need to account for SAO/deblock in the ilr state
    int margin = 0; //4; //Accounts for fracmvest?
    if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF) {
      margin += SAO_DELAY_PX;
    }
    else if (state->ILR_state->encoder_control->cfg.deblock_enable) {
      margin += DEBLOCK_DELAY_PX;
    }
    
    //Check previous block range to avoid re-creating jobs
    if (lcu->above != NULL) {
      int tmp_block_y = lcu->above->encoder_state->tile->offset_y + lcu->above->position_px.y;
      int tmp_block_height = lcu->above->size.y;
      int tmp_range[2];
      kvz_blockScalingSrcHeightRange(tmp_range, state_param->param, tmp_block_y, tmp_block_height);
      set_job_row = ((tmp_range[1] + margin) / LCU_WIDTH) + 1;//(tmp_range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //No need to clip since it works a lower bound and ver_range[1] limits row to the correct range
    }

    //Map the pixel range to LCU row
    //ver_range[0] = ver_range[0] / LCU_WIDTH; //First LCU that is needed
    ver_range[1] = ((ver_range[1] + margin) / LCU_WIDTH) + 1;//(ver_range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
    ver_range[1] = MIN(ver_range[1], state->num_ILR_states);
    
    //Create hor job for each row that does not have one yet and add dep
    int id_offset = lcu->index; //Offset the hor job index by the column number of the current lcu
    //for (int row = ver_range[0]; row < ver_range[1]; row++) {
    const int row = ver_range[1] - 1;
      int hor_ind = row + id_offset * state->num_ILR_states;//row * state->encoder_control->in.width_in_lcu + id_offset;

      if (row >= set_job_row) {
        //Create hor param
        kvz_image_scaling_parameter_t *param_hor = calloc(1, sizeof(kvz_image_scaling_parameter_t));
        kvz_copy_image_scaling_parameters(param_hor, param_ver);
        
        param_hor->pic_in = kvz_image_copy_ref(state_param->pic_in);
        kvz_image_free(param_hor->pic_out);
        param_hor->pic_out = NULL;

        //Set correct block parameters for hor job since it may be on a different lcu row than lcu
        const encoder_state_t * const ilr_state = &state->ILR_state[row];
        param_hor->block_y = ilr_state->tile->offset_y + ilr_state->lcu_order->position_px.y;
        param_hor->block_height = ilr_state->lcu_order->size.y;

        //Account for SAO by only upsampling rows that have been SAOd on this lcu row (exclude bottom most margin pixels if not the last lcu row)
        if( lcu->above != NULL ){
          param_hor->block_y -= margin;
        } else {
          param_hor->block_height -= margin;
        }
        if( lcu->below == NULL ){
          param_hor->block_height += margin;
        }

        kvz_threadqueue_free_job(&state->layer->image_hor_scaling_jobs[hor_ind]);
        state->layer->image_hor_scaling_jobs[hor_ind] = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param_hor);

        //Add dependency to ILR jobs
        //Calculate vertical range of block scaling
        int hor_range[2]; //Range of blocks needed for scaling
        kvz_blockScalingSrcWidthRange(hor_range, state_param->param, param_hor->block_x, param_hor->block_width);

        //Map the pixel range to LCU pos.
        //hor_range[0] = hor_range[0] / LCU_WIDTH; //First LCU that is needed
        hor_range[1] = ((hor_range[1] + margin) / LCU_WIDTH) + 1;//(hor_range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
        hor_range[1] = MIN(hor_range[1], ilr_state->lcu_order_count);
        
        //TODO: Only need to add dependency to last lcu since it already depends on prev lcu?

        //Add ilr state dependencies to hor job
        //for (int k = hor_range[0]; k < hor_range[1]; k++) {
        const int k = hor_range[1] - 1;
        kvz_threadqueue_job_dep_add(state->layer->image_hor_scaling_jobs[hor_ind], ilr_state->tile->wf_jobs[ilr_state->lcu_order[k].id]);
        //}

        //Add dependency to left lcu so that copying to src_buffer is not an issue
        if (lcu->left != NULL) {
          kvz_threadqueue_job_dep_add(state->layer->image_hor_scaling_jobs[hor_ind], state->layer->image_hor_scaling_jobs[hor_ind - state->num_ILR_states /*1*/]);
        }

        //Submit hor job
        kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_hor_scaling_jobs[hor_ind]);
      }

      //Add hor step dependency to ver step
      kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], state->layer->image_hor_scaling_jobs[hor_ind]);
    //}

    //Dependencies added so submit the ver job
    kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_ver_scaling_jobs[lcu->id]);
    
    break;
#else
    //Calculate suitable size for scaling job. Make sure atleast on ctu of the BL is processed
    unsigned job_width = MAX(1, (unsigned)(state_param->param->rnd_trgt_width / state_param->param->scaled_src_width));
    unsigned job_height = MAX(1, (unsigned)(state_param->param->rnd_trgt_height / state_param->param->scaled_src_height)) + 1;

    //Skip creating job if current lcu is already being upsampled
    if (lcu->position.x % job_width != 0 || lcu->position.y % job_height != 0)
    {
      break;
    }

    //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
    kvz_image_scaling_parameter_t *param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
    kvz_copy_image_scaling_parameters(param, state_param);
    param->block_x = state->tile->offset_x + lcu->position_px.x;
    param->block_y = state->tile->offset_y + lcu->position_px.y;
    param->block_width = lcu->size.x;
    param->block_height = lcu->size.y;

    param->use_tiles = 0;

    //Add more lcu to the job based on job_height/width
    const lcu_order_element_t *tmp_lcu = lcu->right;

    for (unsigned i = 1; i < job_width && tmp_lcu != NULL; i++)
    {
      param->block_width += tmp_lcu->size.x;
      tmp_lcu = tmp_lcu->right;
    }

    tmp_lcu = lcu->below;
    for (unsigned i = 1; i < job_height && tmp_lcu != NULL; i++)
    {
      param->block_height += tmp_lcu->size.y;
      tmp_lcu = tmp_lcu->below;
    }

    //Need to account for SAO/deblock in the ilr state
    int margin = 0; //4; //Accounts for fracmvest?
    if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF) {
      margin += SAO_DELAY_PX;
    } else if (state->ILR_state->encoder_control->cfg.deblock_enable) {
      margin += DEBLOCK_DELAY_PX;
    }

    //Account for SAO by only upsampling rows that have been SAOd on this lcu row (exclude bottom most margin pixels if not the last lcu row)
    /*if (lcu->above != NULL) {
      param->block_y -= margin;
    } else {
      param->block_height -= margin;
    }
    if (lcu->below == NULL) {
      param->block_height += margin;
    }

    if (lcu->left != NULL) {
      param->block_x -= margin;
    } else {
      param->block_width -= margin;
    }
    if (lcu->right == NULL) {
      param->block_width += margin;
    }*/

    //First create the job
    kvz_threadqueue_free_job(&state->layer->image_ver_scaling_jobs[lcu->id]);
    state->layer->image_ver_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param);

    //Calculate horizontal range
    int ver_range[2];
    kvz_blockScalingSrcHeightRange(ver_range, state_param->param, param->block_y, param->block_height);
    
    //Map the pixel range to LCU row
    //ver_range[0] = ver_range[0] / LCU_WIDTH; //First LCU that is needed
    ver_range[1] = ((ver_range[1] + margin) / LCU_WIDTH) + 1;//(ver_range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
    ver_range[1] = MIN(ver_range[1], state->num_ILR_states);

    const int row = ver_range[1] - 1;

    //Set correct block parameters for hor job since it may be on a different lcu row than lcu
    const encoder_state_t * const ilr_state = &state->ILR_state[row];

    //Calculate vertical range of block scaling
    int hor_range[2]; //Range of blocks needed for scaling
    kvz_blockScalingSrcWidthRange(hor_range, state_param->param, param->block_x, param->block_width);

    //Map the pixel range to LCU pos.
    //hor_range[0] = hor_range[0] / LCU_WIDTH; //First LCU that is needed
    hor_range[1] = ((hor_range[1] + margin) / LCU_WIDTH) + 1;//(hor_range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
    hor_range[1] = MIN(hor_range[1], ilr_state->lcu_order_count);
    
    //Add dependency to ILR jobs
    //Add ilr state dependencies to hor job
    //for (int k = hor_range[0]; k < hor_range[1]; k++) {
    const int k = hor_range[1] - 1;
    kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], ilr_state->tile->wf_jobs[ilr_state->lcu_order[k].id]);
    //}

    //Add dependency to left and above lcu so that copying to src_buffer is not an issue
    tmp_lcu = lcu->left;
    for (unsigned i = 1; i < job_width && tmp_lcu != NULL; i++)
    {
      tmp_lcu = tmp_lcu->left;
    }
    if (tmp_lcu != NULL) {
      kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], state->layer->image_ver_scaling_jobs[tmp_lcu->id]);
    }

    tmp_lcu = lcu->above;
    for (unsigned i = 1; i < job_height && tmp_lcu != NULL; i++)
    {
      tmp_lcu = tmp_lcu->above;
    }
    if (tmp_lcu != NULL) {
      kvz_threadqueue_job_dep_add(state->layer->image_ver_scaling_jobs[lcu->id], state->layer->image_ver_scaling_jobs[tmp_lcu->id]);
    }

    //Dependencies added so submit the ver job
    kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->image_ver_scaling_jobs[lcu->id]);

    break;

#endif
  }

  case ENCODER_STATE_TYPE_TILE:
  {
    //Make a scaling job for tile

    //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
    kvz_image_scaling_parameter_t * const param = calloc(1, sizeof(kvz_image_scaling_parameter_t));
    kvz_copy_image_scaling_parameters(param, state_param);
    param->block_x = state->tile->offset_x;
    param->block_y = state->tile->offset_y;
    param->block_width = state_param->trgt_buffer->y->width; //Trgt buffer should be the size of the block
    param->block_height = state_param->trgt_buffer->y->height;

    param->use_tiles = 1;

    //Do hor/ver scaling in the same job
    kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done); //Should have been set to NULL anyway
    state->tqj_ilr_rec_scaling_done = kvz_threadqueue_job_create(kvz_block_step_scaler_worker, (void*)param);

    //Need to add dependency to all ilr tiles that that are within the src range
    int range[4];
    kvz_blockScalingSrcWidthRange(range, param->param, param->block_x, param->block_width);
    kvz_blockScalingSrcHeightRange(range + 2, param->param, param->block_y, param->block_height);

    //Need to account for SAO/deblock in the ilr state. TODO: Figure out if it affects tiles
    /*int margin = 4; //Accounts for fracmvest?
    if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF) {
      margin += SAO_DELAY_PX;
    }
    else if (state->ILR_state->encoder_control->cfg.deblock_enable) {
      margin += DEBLOCK_DELAY_PX;
    }*/

    for (int i = 0; i < state->num_ILR_states; i++) {
      const encoder_state_t *ilr_state = &state->ILR_state[i];
      int ilr_tile_x = ilr_state->tile->offset_x;
      int ilr_tile_y = ilr_state->tile->offset_y;
      int ilr_tile_width = ilr_state->tile->frame->width;
      int ilr_tile_height = ilr_state->tile->frame->height;
      int block_x = range[0];
      int block_y = range[2];
      int block_width = range[1] - range[0] + 1;
      int block_height = range[3] - range[2] + 1;

      if (ilr_state->tqj_recon_done != NULL
        && ABS((block_x << 1) + block_width - (ilr_tile_x << 1) - ilr_tile_width) < block_width + ilr_tile_width
        && ABS((block_y << 1) + block_height - (ilr_tile_y << 1) - ilr_tile_height) < block_height + ilr_tile_height) {
        //Areas intersect so add a dependency to ilr tile recon done
        kvz_threadqueue_job_dep_add(state->tqj_ilr_rec_scaling_done, ilr_state->tqj_recon_done);
      }
    }

    //Submit job
    kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tqj_ilr_rec_scaling_done);

    break;
  }

  default:
  {
    //TODO: Handle other types, give error?
    break;
  }
  }//END switch

  //state->layer->scaling_started = 1;
}


//Do necessary preparations for deferred scaling.
//Actual scaling jobs are started later
static kvz_picture* prepare_deferred_block_step_scaling(kvz_picture* const pic_in, const scaling_parameter_t *const param, encoder_state_t * const state, const uint8_t skip_same)
{
  if (pic_in == NULL) {
    return NULL;
  }

  //If no scaling needs to be done, just return pic_in
  if (skip_same && param->src_height == param->trgt_height && param->src_width == param->trgt_width) {
    state->layer->scaling_started = 1;
    return kvz_image_copy_ref(pic_in);
  }

  state->layer->scaling_started = 0;

  //Prepare img_job_param.
  state->layer->img_job_param.param = param;
  kvz_image_free(state->layer->img_job_param.pic_in); //Free prev pic
  state->layer->img_job_param.pic_in = kvz_image_copy_ref(pic_in);

  kvz_image_free(state->layer->img_job_param.pic_out); //Free prev pic
  state->layer->img_job_param.pic_out = kvz_image_alloc(pic_in->chroma_format,
    param->trgt_width + param->trgt_padding_x,
    param->trgt_height + param->trgt_padding_y);

  //Copy other information
  state->layer->img_job_param.pic_out->dts = pic_in->dts;
  state->layer->img_job_param.pic_out->pts = pic_in->pts;
  state->layer->img_job_param.pic_out->interlacing = pic_in->interlacing;

  return kvz_image_copy_ref(state->layer->img_job_param.pic_out);
}

//Loop over a tile
//TODO: make a proper implementation that doesn't just call kvz_cu_array_upsampling_worker
static void tile_cu_array_upsampling_worker(void *opaque_param)
{
  kvz_cua_upsampling_parameter_t *param = opaque_param;

  for (int i = 0; i < param->lcu_order_count; i++) {
    kvz_cua_upsampling_parameter_t *tmp_param = calloc(1, sizeof(kvz_cua_upsampling_parameter_t));
    kvz_copy_cua_upsampling_parameters(tmp_param, param);
    int tile_x = ((encoder_state_config_tile_t*)param->tile)->lcu_offset_x;
    int tile_y = ((encoder_state_config_tile_t*)param->tile)->lcu_offset_y;
    tmp_param->lcu_ind = ((lcu_order_element_t*)param->lcu_order)[i].id + tile_x + tile_y * param->nw_in_lcu; //Calculate the modified ind accounting for tile offset
    tmp_param->lcu_order = NULL;
    tmp_param->tile = NULL;
    tmp_param->lcu_order_count = 0;

    kvz_cu_array_upsampling_worker(tmp_param);
  }

}

//Scale tile specified part of cua
static void cua_lcu_scaling( encoder_state_t * const state )
{
  state->layer->cua_job_param.tile = state->tile;
  state->layer->cua_job_param.lcu_order = state->lcu_order;
  state->layer->cua_job_param.lcu_order_count = state->lcu_order_count;

  tile_cu_array_upsampling_worker(&state->layer->cua_job_param);

  //state->layer->scaling_started = 1;
}

// Start the cu array scaling job for the given state
static void start_cua_lcu_scaling_job(encoder_state_t * const state, const lcu_order_element_t * const lcu)
{
  //If scaling has already been started, no need to do it here anymore
  //if (state->layer == NULL || state->layer->scaling_started) {
  //  return;
  //}

  const kvz_cua_upsampling_parameter_t * const state_param = &state->layer->cua_job_param;
  switch (state->type) {
  case ENCODER_STATE_TYPE_TILE:
  {
    int block_x = state->tile->offset_x + lcu->position_px.x;
    int block_y = state->tile->offset_y + lcu->position_px.y;
    int block_width = lcu->size.x;
    int block_height = lcu->size.y;
    int src_width = state_param->base_cua->width;
    int src_height = state_param->base_cua->height;

    state->layer->cua_job_param.tile = state->tile;
    state->layer->cua_job_param.lcu_order = state->lcu_order;
    state->layer->cua_job_param.lcu_order_count = state->lcu_order_count;

    kvz_threadqueue_free_job(&state->tqj_ilr_cua_upsampling_done);
    state->tqj_ilr_cua_upsampling_done = kvz_threadqueue_job_create(tile_cu_array_upsampling_worker, (void*)state_param);

    //Calculate (vertical/horizontal) range of scaling
    int range[4]; //Range of blocks needed for scaling
    kvz_cu_array_upsampling_src_range(range, block_x, block_x + block_width - 1, src_width, state_param->cu_pos_scale[0]); //Width  
    kvz_cu_array_upsampling_src_range(range + 2, block_y, block_y + block_height - 1, src_height, state_param->cu_pos_scale[1]); //Height

    //Need to account for SAO/deblock in the ilr state. TODO: Figure out if it affects tiles
    /*int margin = 4; //Accounts for fracmvest?
    if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF) {
      margin += SAO_DELAY_PX;
    } else if (state->ILR_state->encoder_control->cfg.deblock_enable) {
      margin += DEBLOCK_DELAY_PX;
    }*/

    for (int i = 0; i < state->num_ILR_states; i++) {
      const encoder_state_t *ilr_state = &state->ILR_state[i];
      int ilr_tile_x = ilr_state->tile->offset_x;
      int ilr_tile_y = ilr_state->tile->offset_y;
      int ilr_tile_width = ilr_state->tile->frame->width;
      int ilr_tile_height = ilr_state->tile->frame->height;
      int block_x = range[0];
      int block_y = range[2];
      int block_width = range[1] - range[0] + 1;
      int block_height = range[3] - range[2] + 1;

      if (ilr_state->tqj_recon_done != NULL
        && ABS((block_x << 1) + block_width - (ilr_tile_x << 1) - ilr_tile_width) < block_width + ilr_tile_width
        && ABS((block_y << 1) + block_height - (ilr_tile_y << 1) - ilr_tile_height) < block_height + ilr_tile_height) {
        //Areas intersect so add a dependency to ilr tile recon done
        kvz_threadqueue_job_dep_add(state->tqj_ilr_cua_upsampling_done, ilr_state->tqj_recon_done);
      }
    }

    //Dependencies added so submit the job
    kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tqj_ilr_cua_upsampling_done);

    break;
  }

  case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
  {
    //Allocate new scaling parameters to pass to the worker and set block info. Worker is in charge of deallocation.
    kvz_cua_upsampling_parameter_t *param = calloc(1, sizeof(kvz_cua_upsampling_parameter_t));
    kvz_copy_cua_upsampling_parameters(param, state_param);

    int block_x = state->tile->offset_x + lcu->position_px.x;
    int block_y = state->tile->offset_y + lcu->position_px.y;
    int block_width = lcu->size.x;
    int block_height = lcu->size.y;
    int src_width = state_param->base_cua->width;
    int src_height = state_param->base_cua->height;
    int tile_x = state->tile->lcu_offset_x;
    int tile_y = state->tile->lcu_offset_y;

    param->lcu_ind = lcu->id + tile_x + tile_y * param->nw_in_lcu; //Calculate the modified ind accounting for tile offset

    kvz_threadqueue_free_job(&state->layer->cua_scaling_jobs[lcu->id]);
    state->layer->cua_scaling_jobs[lcu->id] = kvz_threadqueue_job_create(kvz_cu_array_upsampling_worker, (void*)param);

    //Calculate (vertical/horizontal) range of scaling
    int range[4]; //Range of blocks needed for scaling
    kvz_cu_array_upsampling_src_range(range, block_x, block_x + block_width - 1, src_width, param->cu_pos_scale[0]); //Width  
    kvz_cu_array_upsampling_src_range(range + 2, block_y, block_y + block_height - 1, src_height, param->cu_pos_scale[1]); //Height

    //Need to account for SAO/deblock in the ilr state.
    int margin = 0; //4; //Accounts for fracmvest?
    if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF) {
      margin += SAO_DELAY_PX;
    } else if (state->ILR_state->encoder_control->cfg.deblock_enable) {
      margin += DEBLOCK_DELAY_PX;
    }

    //Map the pixel range to LCU pos
    //range[0] = range[0] / LCU_WIDTH; //First LCU that is needed
    range[1] = ((range[1] + margin) / LCU_WIDTH) + 1;//(range[1] + margin + LCU_WIDTH - 1) / LCU_WIDTH; //Last LCU that is not needed
    //range[2] = range[2] / LCU_WIDTH;
    range[3] = ((range[3] + margin) / LCU_WIDTH) + 1;//(range[3] + margin + LCU_WIDTH - 1) / LCU_WIDTH;

    //TODO: Figure out correct dependency.
    //range[2] = 0;//MAX(0, range[2] - 1); //0;
    range[3] = MIN(range[3] + state->encoder_control->max_inter_ref_lcu.down, state->ILR_state->tile->frame->height_in_lcu);  //state->ILR_state->tile->frame->height_in_lcu;//MIN(range[3]+3, state->ILR_state->tile->frame->height_in_lcu);//state->ILR_state->tile->frame->height_in_lcu;
    //range[0] = 0; //MAX(0, range[0] - 1); //0;
    range[1] = MIN(range[1] + state->encoder_control->max_inter_ref_lcu.right, state->ILR_state->tile->frame->width_in_lcu); //state->ILR_state->tile->frame->width_in_lcu;//MIN(range[1]+3, state->ILR_state->tile->frame->width_in_lcu); //state->ILR_state->tile->frame->width_in_lcu;
    //TODO: Only need to add dependency to last lcu since it already depends on prev lcu?
    //Add dependencies to ilr states
    /*for (int j = range[2]; j < range[3]; j++) {
      const encoder_state_t * const ilr_state = &state->ILR_state[j];
      for (int k = range[0]; k < range[1]; k++) {
        kvz_threadqueue_job_dep_add(state->layer->cua_scaling_jobs[lcu->id], ilr_state->tile->wf_jobs[ilr_state->lcu_order[k].id]);
      }
    }*/
    const encoder_state_t *const ilr_state = &state->ILR_state[range[3] - 1];
    kvz_threadqueue_job_dep_add(state->layer->cua_scaling_jobs[lcu->id], ilr_state->tile->wf_jobs[ilr_state->lcu_order[range[1] - 1].id]);

    //Dependencies added so submit the job
    kvz_threadqueue_submit(state->encoder_control->threadqueue, state->layer->cua_scaling_jobs[lcu->id]);
    
    break;
  }

  default:
  {
    //TODO: Handle other types?
    break;
  }
  }//END switch

  //state->layer->scaling_started = 1;
}


//Do necessary preparations for deferred scaling.
//Actual scaling jobs are started later
static cu_array_t* prepare_deferred_cua_lcu_upsampling(encoder_state_t * const state, const int32_t mv_scale[2], const int32_t pos_scale[2], const uint8_t skip_same)
{
  if (skip_same && pos_scale[0] == POS_SCALE_FAC_1X && pos_scale[1] == POS_SCALE_FAC_1X) {
    state->layer->scaling_started = 1;
    return kvz_cu_array_copy_ref(state->ILR_state->tile->frame->cu_array);
  }

  state->layer->scaling_started = 0;

  //Allocate the new cua.
  uint32_t n_width = state->tile->frame->width_in_lcu * LCU_WIDTH;
  uint32_t n_height = state->tile->frame->height_in_lcu * LCU_WIDTH;
  cu_array_t *cua = kvz_cu_array_alloc(n_width, n_height);

  //Allocate scaling parameters to give to the worker. Worker should handle freeing.
  kvz_cua_upsampling_parameter_t *param = &state->layer->cua_job_param;
  kvz_cu_array_free(&param->base_cua); //Free prev
  param->base_cua = kvz_cu_array_copy_ref(state->ILR_state->tile->frame->cu_array);
  param->cu_pos_scale[0] = pos_scale[0];
  param->cu_pos_scale[1] = pos_scale[1];
  param->mv_scale[0] = mv_scale[0];
  param->mv_scale[1] = mv_scale[1];
  param->nh_in_lcu = state->tile->frame->height_in_lcu;
  param->nw_in_lcu = state->tile->frame->width_in_lcu;
  param->only_init = !state->encoder_control->cfg.tmvp_enable;
  kvz_cu_array_free(&param->out_cua); //Free prev
  param->out_cua = kvz_cu_array_copy_ref(cua);
  param->lcu_ind = -1; //Set correct ind later

  return cua;
}

/**
* Handle adding ilr frames to the ref list.
*
* - Add ilr frame to the ref list if ilr is enabled
* - Prepare scaling jobs if threads are used or do scaling
*/
static void prepare_ilr_frames(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
    //TODO: Account for adding several ILR frames. Should ilr rec ever bee NULL?
  if (encoder->cfg.ILR_frames > 0 && state->ILR_state != NULL && state->ILR_state->tile->frame->rec != NULL) {
    //Add base layer to the reference list.
    const encoder_state_t *ILR_state = state->ILR_state;
    kvz_picture *ilr_rec = kvz_image_copy_ref(ILR_state->tile->frame->rec);
    kvz_picture *scaled_pic = NULL;
    
    //TODO: Account for offsets etc. Need to use something else than original sizes?
    //TODO: Shm doesn't upsample the cua if frame is idr; need to do something else as well when idr? Does it even matter?
    int32_t mv_scale[2] = { GET_SCALE_MV(encoder->layer.upscaling.src_width,encoder->layer.upscaling.trgt_width),
      GET_SCALE_MV(encoder->layer.upscaling.src_height,encoder->layer.upscaling.trgt_height) };
    int32_t pos_scale[2] = { GET_SCALE_POS(encoder->layer.upscaling.src_width,encoder->layer.upscaling.trgt_width),
      GET_SCALE_POS(encoder->layer.upscaling.src_height,encoder->layer.upscaling.trgt_height) };
    cu_array_t* scaled_cu = NULL;    
    
    if (encoder->cfg.threads > 0) {
      scaled_pic = prepare_deferred_block_step_scaling(ilr_rec, &encoder->layer.upscaling, state, 1);
      scaled_cu = prepare_deferred_cua_lcu_upsampling(state, mv_scale, pos_scale, 1);
    }
    else {
      scaled_pic = kvz_image_scaling(ilr_rec, &encoder->layer.upscaling, 1);
      scaled_cu = kvz_cu_array_upsampling(ILR_state->tile->frame->cu_array,
        state->tile->frame->width_in_lcu,
        state->tile->frame->height_in_lcu,
        mv_scale, pos_scale, 1,
        !state->encoder_control->cfg.tmvp_enable);//(state->frame->pictype >= KVZ_NAL_BLA_W_LP && state->frame->pictype <= KVZ_NAL_CRA_NUT)); //pic type not set yet 
      state->layer->scaling_started = 1;
    }

    if (ilr_rec == NULL || scaled_pic == NULL || scaled_cu == NULL) {
      kvz_image_free(ilr_rec);
      kvz_image_free(scaled_pic);
      kvz_cu_array_free(&scaled_cu);
      return; //TODO: Add error etc?
    }

    //Copy image ref info
    //TODO: Is there a way to determine how many refs ilr_rec has? Otherwise just copy everything to be safe
    memcpy(scaled_pic->ref_pocs, ilr_rec->ref_pocs, sizeof(ilr_rec->ref_pocs));//sizeof(int32_t) * ILR_state->frame->ref->used_size);
    memcpy(scaled_pic->picture_info, ilr_rec->picture_info, sizeof(ilr_rec->picture_info));//sizeof(kvz_picture_info_t) * ILR_state->frame->ref->used_size);
    kvz_image_free(ilr_rec);

    kvz_image_list_add(state->frame->ref,
      scaled_pic,
      scaled_cu,
      ILR_state->frame->poc,
      ILR_state->frame->ref_LX,
      ILR_state->encoder_control->cfg.gop[ILR_state->frame->gop_offset].tId,
      ILR_state->encoder_control->layer.layer_id,
      1); //Currently only ILR can be a long term references
          //TODO: Add error handling?
    kvz_image_free(scaled_pic);
    kvz_cu_array_free(&scaled_cu);
  }
}
// ***********************************************

/**
 * \brief Save edge pixels before SAO to buffers.
 *
 * Copies pixels at the edges of the area that will be filtered with SAO to
 * the given buffers. If deblocking is enabled, the pixels must have been
 * deblocked before this.
 *
 * The saved pixels will be needed later when doing SAO for the neighboring
 * areas.
 */
static void encoder_state_recdata_before_sao_to_bufs(
    encoder_state_t * const state,
    const lcu_order_element_t * const lcu,
    yuv_t * const hor_buf,
    yuv_t * const ver_buf)
{
  videoframe_t* const frame = state->tile->frame;

  if (hor_buf && lcu->below) {
    // Copy the bottommost row that will be filtered with SAO to the
    // horizontal buffer.
    vector2d_t pos = {
      .x = lcu->position_px.x,
      .y = lcu->position_px.y + LCU_WIDTH - SAO_DELAY_PX - 1,
    };
    // Copy all pixels that have been deblocked.
    int length = lcu->size.x - DEBLOCK_DELAY_PX;

    if (!lcu->right) {
      // If there is no LCU to the right, the last pixels will be
      // filtered too.
      length += DEBLOCK_DELAY_PX;
    }

    if (lcu->left) {
      // The rightmost pixels of the CTU to the left will also be filtered.
      pos.x -= DEBLOCK_DELAY_PX;
      length += DEBLOCK_DELAY_PX;
    }

    const unsigned from_index = pos.x + pos.y * frame->rec->stride;
    // NOTE: The horizontal buffer is indexed by
    //    x_px + y_lcu * frame->width
    // where x_px is in pixels and y_lcu in number of LCUs.
    const unsigned to_index = pos.x + lcu->position.y * frame->width;

    kvz_pixels_blit(&frame->rec->y[from_index],
                    &hor_buf->y[to_index],
                    length, 1,
                    frame->rec->stride,
                    frame->width);

    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      const unsigned from_index_c = (pos.x / 2) + (pos.y / 2) * frame->rec->stride / 2;
      const unsigned to_index_c = (pos.x / 2) + lcu->position.y * frame->width / 2;

      kvz_pixels_blit(&frame->rec->u[from_index_c],
                      &hor_buf->u[to_index_c],
                      length / 2, 1,
                      frame->rec->stride / 2,
                      frame->width / 2);
      kvz_pixels_blit(&frame->rec->v[from_index_c],
                      &hor_buf->v[to_index_c],
                      length / 2, 1,
                      frame->rec->stride / 2,
                      frame->width / 2);
    }
  }

  if (ver_buf && lcu->right) {
    // Copy the rightmost column that will be filtered with SAO to the
    // vertical buffer.
    vector2d_t pos = {
      .x = lcu->position_px.x + LCU_WIDTH - SAO_DELAY_PX - 1,
      .y = lcu->position_px.y,
    };
    int length = lcu->size.y - DEBLOCK_DELAY_PX;

    if (!lcu->below) {
      // If there is no LCU below, the last pixels will be filtered too.
      length += DEBLOCK_DELAY_PX;
    }

    if (lcu->above) {
      // The bottommost pixels of the CTU above will also be filtered.
      pos.y -= DEBLOCK_DELAY_PX;
      length += DEBLOCK_DELAY_PX;
    }

    const unsigned from_index = pos.x + pos.y * frame->rec->stride;
    // NOTE: The vertical buffer is indexed by
    //    x_lcu * frame->height + y_px
    // where x_lcu is in number of LCUs and y_px in pixels.
    const unsigned to_index = lcu->position.x * frame->height + pos.y;

    kvz_pixels_blit(&frame->rec->y[from_index],
                    &ver_buf->y[to_index],
                    1, length,
                    frame->rec->stride, 1);

    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      const unsigned from_index_c = (pos.x / 2) + (pos.y / 2) * frame->rec->stride / 2;
      const unsigned to_index_c = lcu->position.x * frame->height / 2 + pos.y / 2;

      kvz_pixels_blit(&frame->rec->u[from_index_c],
                      &ver_buf->u[to_index_c],
                      1, length / 2,
                      frame->rec->stride / 2, 1);
      kvz_pixels_blit(&frame->rec->v[from_index_c],
                      &ver_buf->v[to_index_c],
                      1, length / 2,
                      frame->rec->stride / 2, 1);
    }
  }
}

static void encoder_state_recdata_to_bufs(encoder_state_t * const state,
                                          const lcu_order_element_t * const lcu,
                                          yuv_t * const hor_buf,
                                          yuv_t * const ver_buf)
{
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

/**
 * \brief Do SAO reconstuction for all available pixels.
 *
 * Does SAO reconstruction for all pixels that are available after the
 * given LCU has been deblocked. This means the following pixels:
 *  - bottom-right block of SAO_DELAY_PX times SAO_DELAY_PX in the lcu to
 *    the left and up
 *  - the rightmost SAO_DELAY_PX pixels of the LCU to the left (excluding
 *    the bottommost pixel)
 *  - the bottommost SAO_DELAY_PX pixels of the LCU above (excluding the
 *    rightmost pixels)
 *  - all pixels inside the LCU, excluding the rightmost SAO_DELAY_PX and
 *    bottommost SAO_DELAY_PX
 */
static void encoder_sao_reconstruct(const encoder_state_t *const state,
                                    const lcu_order_element_t *const lcu)
{
  videoframe_t *const frame = state->tile->frame;


  // Temporary buffers for SAO input pixels. The buffers cover the pixels
  // inside the LCU (LCU_WIDTH x LCU_WIDTH), SAO_DELAY_PX wide bands to the
  // left and above the LCU, and one pixel border on the left and top
  // sides. We add two extra pixels to the buffers because the AVX2 SAO
  // reconstruction reads up to two extra bytes when using edge SAO in the
  // horizontal direction.
#define SAO_BUF_WIDTH   (1 + SAO_DELAY_PX   + LCU_WIDTH)
#define SAO_BUF_WIDTH_C (1 + SAO_DELAY_PX/2 + LCU_WIDTH_C)
  kvz_pixel sao_buf_y_array[SAO_BUF_WIDTH   * SAO_BUF_WIDTH   + 2];
  kvz_pixel sao_buf_u_array[SAO_BUF_WIDTH_C * SAO_BUF_WIDTH_C + 2];
  kvz_pixel sao_buf_v_array[SAO_BUF_WIDTH_C * SAO_BUF_WIDTH_C + 2];

  // Pointers to the top-left pixel of the LCU in the buffers.
  kvz_pixel *const sao_buf_y = &sao_buf_y_array[(SAO_DELAY_PX + 1) * (SAO_BUF_WIDTH + 1)];
  kvz_pixel *const sao_buf_u = &sao_buf_u_array[(SAO_DELAY_PX/2 + 1) * (SAO_BUF_WIDTH_C + 1)];
  kvz_pixel *const sao_buf_v = &sao_buf_v_array[(SAO_DELAY_PX/2 + 1) * (SAO_BUF_WIDTH_C + 1)];

  const int x_offsets[3] = {
    // If there is an lcu to the left, we need to filter its rightmost
    // pixels.
    lcu->left ? -SAO_DELAY_PX : 0,
    0,
    // If there is an lcu to the right, the rightmost pixels of this LCU
    // are filtered when filtering that LCU. Otherwise we filter them now.
    lcu->size.x - (lcu->right ? SAO_DELAY_PX : 0),
  };

  const int y_offsets[3] = {
    // If there is an lcu above, we need to filter its bottommost pixels.
    lcu->above ? -SAO_DELAY_PX : 0,
    0,
    // If there is an lcu below, the bottommost pixels of this LCU are
    // filtered when filtering that LCU. Otherwise we filter them now.
    lcu->size.y - (lcu->below ? SAO_DELAY_PX : 0),
  };

  // Number of pixels around the block that need to be copied to the
  // buffers.
  const int border_left  = lcu->left  ? 1 : 0;
  const int border_right = lcu->right ? 1 : 0;
  const int border_above = lcu->above ? 1 : 0;
  const int border_below = lcu->below ? 1 : 0;

  // Index of the pixel at the intersection of the top and left borders.
  const int border_index = (x_offsets[0] - border_left) +
                           (y_offsets[0] - border_above) * SAO_BUF_WIDTH;
  const int border_index_c = (x_offsets[0]/2 - border_left) +
                             (y_offsets[0]/2 - border_above) * SAO_BUF_WIDTH_C;
  // Width and height of the whole area to filter.
  const int width  = x_offsets[2] - x_offsets[0];
  const int height = y_offsets[2] - y_offsets[0];

  // Copy bordering pixels from above and left to buffers.
  if (lcu->above) {
    const int from_index = (lcu->position_px.x + x_offsets[0] - border_left) +
                           (lcu->position.y - 1) * frame->width;
    kvz_pixels_blit(&state->tile->hor_buf_before_sao->y[from_index],
                    &sao_buf_y[border_index],
                    width + border_left + border_right,
                    1,
                    frame->width,
                    SAO_BUF_WIDTH);
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      const int from_index_c = (lcu->position_px.x + x_offsets[0])/2 - border_left +
                               (lcu->position.y - 1) * frame->width/2;
      kvz_pixels_blit(&state->tile->hor_buf_before_sao->u[from_index_c],
                      &sao_buf_u[border_index_c],
                      width/2 + border_left + border_right,
                      1,
                      frame->width/2,
                      SAO_BUF_WIDTH_C);
      kvz_pixels_blit(&state->tile->hor_buf_before_sao->v[from_index_c],
                      &sao_buf_v[border_index_c],
                      width/2 + border_left + border_right,
                      1,
                      frame->width/2,
                      SAO_BUF_WIDTH_C);
    }
  }
  if (lcu->left) {
    const int from_index = (lcu->position.x - 1) * frame->height +
                           (lcu->position_px.y + y_offsets[0] - border_above);
    kvz_pixels_blit(&state->tile->ver_buf_before_sao->y[from_index],
                    &sao_buf_y[border_index],
                    1,
                    height + border_above + border_below,
                    1,
                    SAO_BUF_WIDTH);
    if (state->encoder_control->chroma_format != KVZ_CSP_400) {
      const int from_index_c = (lcu->position.x - 1) * frame->height/2 +
                               (lcu->position_px.y + y_offsets[0])/2 - border_above;
      kvz_pixels_blit(&state->tile->ver_buf_before_sao->u[from_index_c],
                      &sao_buf_u[border_index_c],
                      1,
                      height/2 + border_above + border_below,
                      1,
                      SAO_BUF_WIDTH_C);
      kvz_pixels_blit(&state->tile->ver_buf_before_sao->v[from_index_c],
                      &sao_buf_v[border_index_c],
                      1,
                      height/2 + border_above + border_below,
                      1,
                      SAO_BUF_WIDTH_C);
    }
  }
  // Copy pixels that will be filtered and bordering pixels from right and
  // below.
  const int from_index = (lcu->position_px.x + x_offsets[0]) +
                         (lcu->position_px.y + y_offsets[0]) * frame->rec->stride;
  const int to_index = x_offsets[0] + y_offsets[0] * SAO_BUF_WIDTH;
  kvz_pixels_blit(&frame->rec->y[from_index],
                  &sao_buf_y[to_index],
                  width + border_right,
                  height + border_below,
                  frame->rec->stride,
                  SAO_BUF_WIDTH);
  if (state->encoder_control->chroma_format != KVZ_CSP_400) {
    const int from_index_c = (lcu->position_px.x + x_offsets[0])/2 +
                             (lcu->position_px.y + y_offsets[0])/2 * frame->rec->stride/2;
    const int to_index_c = x_offsets[0]/2 + y_offsets[0]/2 * SAO_BUF_WIDTH_C;
    kvz_pixels_blit(&frame->rec->u[from_index_c],
                    &sao_buf_u[to_index_c],
                    width/2 + border_right,
                    height/2 + border_below,
                    frame->rec->stride/2,
                    SAO_BUF_WIDTH_C);
    kvz_pixels_blit(&frame->rec->v[from_index_c],
                    &sao_buf_v[to_index_c],
                    width/2 + border_right,
                    height/2 + border_below,
                    frame->rec->stride/2,
                    SAO_BUF_WIDTH_C);
  }

  // We filter the pixels in four parts:
  //  1. Pixels that belong to the LCU above and to the left
  //  2. Pixels that belong to the LCU above
  //  3. Pixels that belong to the LCU to the left
  //  4. Pixels that belong to the current LCU
  for (int y_offset_index = 0; y_offset_index < 2; y_offset_index++) {
    for (int x_offset_index = 0; x_offset_index < 2; x_offset_index++) {
      const int x = x_offsets[x_offset_index];
      const int y = y_offsets[y_offset_index];
      const int width = x_offsets[x_offset_index + 1] - x;
      const int height = y_offsets[y_offset_index + 1] - y;

      if (width == 0 || height == 0) continue;

      const int lcu_x = (lcu->position_px.x + x) >> LOG2_LCU_WIDTH;
      const int lcu_y = (lcu->position_px.y + y) >> LOG2_LCU_WIDTH;
      const int lcu_index = lcu_x + lcu_y * frame->width_in_lcu;
      const sao_info_t *sao_luma   = &frame->sao_luma[lcu_index];
      const sao_info_t *sao_chroma = &frame->sao_chroma[lcu_index];

      kvz_sao_reconstruct(state,
                          &sao_buf_y[x + y * SAO_BUF_WIDTH],
                          SAO_BUF_WIDTH,
                          lcu->position_px.x + x,
                          lcu->position_px.y + y,
                          width,
                          height,
                          sao_luma,
                          COLOR_Y);

      if (state->encoder_control->chroma_format != KVZ_CSP_400) {
        // Coordinates in chroma pixels.
        int x_c = x >> 1;
        int y_c = y >> 1;

        kvz_sao_reconstruct(state,
                            &sao_buf_u[x_c + y_c * SAO_BUF_WIDTH_C],
                            SAO_BUF_WIDTH_C,
                            lcu->position_px.x / 2 + x_c,
                            lcu->position_px.y / 2 + y_c,
                            width / 2,
                            height / 2,
                            sao_chroma,
                            COLOR_U);
        kvz_sao_reconstruct(state,
                            &sao_buf_v[x_c + y_c * SAO_BUF_WIDTH_C],
                            SAO_BUF_WIDTH_C,
                            lcu->position_px.x / 2 + x_c,
                            lcu->position_px.y / 2 + y_c,
                            width / 2,
                            height / 2,
                            sao_chroma,
                            COLOR_V);
      }
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
 * The QPs are used in deblocking and QP prediction.
 *
 * The QP delta for a quantization group is coded when the first CU with
 * coded block flag set is encountered. Hence, for the purposes of
 * deblocking and QP prediction, all CUs in before the first one that has
 * cbf set use the QP predictor and all CUs after that use (QP predictor
 * + QP delta).
 *
 * \param state           encoder state
 * \param x               x-coordinate of the left edge of the root CU
 * \param y               y-coordinate of the top edge of the root CU
 * \param depth           depth in the CU quadtree
 * \param last_qp         QP of the last CU in the last quantization group
 * \param prev_qp         -1 if QP delta has not been coded in current QG,
 *                        otherwise the QP of the current QG
 */
static void set_cu_qps(encoder_state_t *state, int x, int y, int depth, int *last_qp, int *prev_qp)
{

  // Stop recursion if the CU is completely outside the frame.
  if (x >= state->tile->frame->width || y >= state->tile->frame->height) return;

  cu_info_t *cu = kvz_cu_array_at(state->tile->frame->cu_array, x, y);
  const int cu_width = LCU_WIDTH >> depth;

  if (depth <= state->encoder_control->max_qp_delta_depth) {
    *prev_qp = -1;
  }

  if (cu->depth > depth) {
    // Recursively process sub-CUs.
    const int d = cu_width >> 1;
    set_cu_qps(state, x,     y,     depth + 1, last_qp, prev_qp);
    set_cu_qps(state, x + d, y,     depth + 1, last_qp, prev_qp);
    set_cu_qps(state, x,     y + d, depth + 1, last_qp, prev_qp);
    set_cu_qps(state, x + d, y + d, depth + 1, last_qp, prev_qp);

  } else {
    bool cbf_found = *prev_qp >= 0;

    if (cu->tr_depth > depth) {
      // The CU is split into smaller transform units. Check whether coded
      // block flag is set for any of the TUs.
      const int tu_width = LCU_WIDTH >> cu->tr_depth;
      for (int y_scu = y; !cbf_found && y_scu < y + cu_width; y_scu += tu_width) {
        for (int x_scu = x; !cbf_found && x_scu < x + cu_width; x_scu += tu_width) {
          cu_info_t *tu = kvz_cu_array_at(state->tile->frame->cu_array, x_scu, y_scu);
          if (cbf_is_set_any(tu->cbf, cu->depth)) {
            cbf_found = true;
          }
        }
      }
    } else if (cbf_is_set_any(cu->cbf, cu->depth)) {
      cbf_found = true;
    }

    int8_t qp;
    if (cbf_found) {
      *prev_qp = qp = cu->qp;
    } else {
      qp = kvz_get_cu_ref_qp(state, x, y, *last_qp);
    }

    // Set the correct QP for all state->tile->frame->cu_array elements in
    // the area covered by the CU.
    for (int y_scu = y; y_scu < y + cu_width; y_scu += SCU_WIDTH) {
      for (int x_scu = x; x_scu < x + cu_width; x_scu += SCU_WIDTH) {
        kvz_cu_array_at(state->tile->frame->cu_array, x_scu, y_scu)->qp = qp;
      }
    }

    if (is_last_cu_in_qg(state, x, y, depth)) {
      *last_qp = cu->qp;
    }
  }
}

//Debug stuff for printing thread info
#if 1 && defined(linux)
#define _GNU_SOURCES
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#define PRINT_TID_LCU_JOB_INFO(x,y,w,h,ind,poc,lid) fprintf(stderr, "TID: %ld, pos: (%d,%d), size: (%d,%d), lcu_ind: %d, poc: %d, lid: %d\n", syscall(SYS_gettid), x, y, w, h, ind, poc, lid)
#else
#define PRINT_TID_LCU_JOB_INFO(x,y,w,h,ind,poc,lid)
#endif

static void encoder_state_worker_encode_lcu(void * opaque)
{
  const lcu_order_element_t * const lcu = opaque;
  encoder_state_t *state = lcu->encoder_state;
  const encoder_control_t * const encoder = state->encoder_control;
  videoframe_t* const frame = state->tile->frame;

  PRINT_TID_LCU_JOB_INFO(lcu->position_px.x, lcu->position_px.y, lcu->size.x, lcu->size.y, lcu->id, frame->poc, encoder->layer.layer_id);

  kvz_set_lcu_lambda_and_qp(state, lcu->position);

  lcu_coeff_t coeff;
  state->coeff = &coeff;

  //This part doesn't write to bitstream, it's only search, deblock and sao
  kvz_search_lcu(state, lcu->position_px.x, lcu->position_px.y, state->tile->hor_buf_search, state->tile->ver_buf_search);

  encoder_state_recdata_to_bufs(state, lcu, state->tile->hor_buf_search, state->tile->ver_buf_search);

  if (encoder->max_qp_delta_depth >= 0) {
    int last_qp = state->last_qp;
    int prev_qp = -1;
    set_cu_qps(state, lcu->position_px.x, lcu->position_px.y, 0, &last_qp, &prev_qp);
  }

  if (encoder->cfg.deblock_enable) {
    kvz_filter_deblock_lcu(state, lcu->position_px.x, lcu->position_px.y);
  }

  if (encoder->cfg.sao_type) {
    // Save the post-deblocking but pre-SAO pixels of the LCU to a buffer
    // so that they can be used in SAO reconstruction later.
    encoder_state_recdata_before_sao_to_bufs(state,
                                             lcu,
                                             state->tile->hor_buf_before_sao,
                                             state->tile->ver_buf_before_sao);
    kvz_sao_search_lcu(state, lcu->position.x, lcu->position.y);
    encoder_sao_reconstruct(state, lcu);
  }

  //Now write data to bitstream (required to have a correct CABAC state)
  const uint64_t existing_bits = kvz_bitstream_tell(&state->stream);

  //Encode SAO
  if (encoder->cfg.sao_type) {
    encode_sao(state, lcu->position.x, lcu->position.y, &frame->sao_luma[lcu->position.y * frame->width_in_lcu + lcu->position.x], &frame->sao_chroma[lcu->position.y * frame->width_in_lcu + lcu->position.x]);
  }

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
}

static void encoder_state_encode_leaf(encoder_state_t * const state)
{
  assert(state->is_leaf);
  assert(state->lcu_order_count > 0);

  const encoder_control_t *ctrl = state->encoder_control;
  const kvz_config *cfg = &ctrl->cfg;

  state->last_qp = state->frame->QP;

  if (cfg->crypto_features) {
    state->crypto_hdl = kvz_crypto_create(cfg);
    state->crypto_prev_pos = 0;
  }

  // Select whether to encode the frame/tile in current thread or to define
  // wavefront jobs for other threads to handle.
  bool wavefront = state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW;
  bool use_parallel_encoding = (wavefront && state->parent->children[1].encoder_control);
  if (!use_parallel_encoding) {
    //*********************************************
    //For scalable extension.
    // Need to check if scaling jobs have been started. If not, need to do scaling first. 
    if( state->layer != NULL && !state->layer->scaling_started ){
      block_step_scaling(state);
      cua_lcu_scaling(state);
      state->layer->scaling_started = 1; //Propably no need to set it here, but shouldn't hurt either. Only set if one wavefront row or using tiles in which case layer struct should be separate.
    }
    //*********************************************

    // Encode every LCU in order and perform SAO reconstruction after every
    // frame is encoded. Deblocking and SAO search is done during LCU encoding.

    for (int i = 0; i < state->lcu_order_count; ++i) {
      encoder_state_worker_encode_lcu(&state->lcu_order[i]);
    }
  } else {
    // Add each LCU in the wavefront row as it's own job to the queue.

    // Select which frame dependancies should be set to.
    const encoder_state_t * ref_state = NULL;

    if (state->frame->slicetype == KVZ_SLICE_I) {
      // I-frames have no references.
      ref_state = NULL;
    } else if (cfg->gop_lowdelay &&
               cfg->gop_len > 0 &&
               state->previous_encoder_state != state)
    {
      // For LP-gop, depend on the state of the first reference.
      int ref_neg = cfg->gop[state->frame->gop_offset].ref_neg[0];
      if (ref_neg > cfg->owf) {
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

      kvz_threadqueue_free_job(&state->tile->wf_jobs[lcu->id]);
      state->tile->wf_jobs[lcu->id] = kvz_threadqueue_job_create(encoder_state_worker_encode_lcu, (void*)lcu);
      threadqueue_job_t **job = &state->tile->wf_jobs[lcu->id];

      // If job object was returned, add dependancies and allow it to run.
      if (job[0]) {
        // Add inter frame dependancies when ecoding more than one frame at
        // once. The added dependancy is for the first LCU of each wavefront
        // row to depend on the reconstruction status of the row below in the
        // previous frame.
        if (ref_state != NULL &&
            state->previous_encoder_state->tqj_recon_done &&
            state->frame->slicetype != KVZ_SLICE_I)
        {
          // We need to wait until the CTUs whose pixels we refer to are
          // done before we can start this CTU.
          const lcu_order_element_t *dep_lcu = lcu;
          for (int i = 0; dep_lcu->below && i < ctrl->max_inter_ref_lcu.down; i++) {
            dep_lcu = dep_lcu->below;
          }
          for (int i = 0; dep_lcu->right && i < ctrl->max_inter_ref_lcu.right; i++) {
            dep_lcu = dep_lcu->right;
          }
          kvz_threadqueue_job_dep_add(job[0], ref_state->tile->wf_jobs[dep_lcu->id]);
        }

        // Add local WPP dependancy to the LCU on the left.
        if (lcu->left) {
          kvz_threadqueue_job_dep_add(job[0], job[-1]);
        }
        // Add local WPP dependancy to the LCU on the top right.
        if (lcu->above) {
          if (lcu->above->right) {
            kvz_threadqueue_job_dep_add(job[0], job[-state->tile->frame->width_in_lcu + 1]);
          } else {
            kvz_threadqueue_job_dep_add(job[0], job[-state->tile->frame->width_in_lcu]);
          }
        }

        //*********************************************
        //For scalable extension.
        //Add dependency to ilr recon upscaling

        if (state->ILR_state != NULL) {

          //Set up scaling jobs
          if (state->layer != NULL && !state->layer->scaling_started) {
            start_block_step_scaling_job(state, lcu);
            start_cua_lcu_scaling_job(state, lcu);
            //Don't set scaling started to true here since it prevents scaling in other lcu and wavefront rows
          }

          //Add a direct dependecy from ilr states wf_job to this (only for SNR)
          if (state->encoder_control->cfg.width == state->ILR_state->encoder_control->cfg.width &&
            state->encoder_control->cfg.width == state->ILR_state->encoder_control->cfg.width) {
            //Account for deblock/SAO
            const lcu_order_element_t *ilr_lcu = lcu;
            if (state->ILR_state->encoder_control->cfg.sao_type != KVZ_SAO_OFF || state->ILR_state->encoder_control->cfg.deblock_enable) {
              ilr_lcu = (ilr_lcu->below != NULL) ? lcu->below : ilr_lcu;
              ilr_lcu = (ilr_lcu->right != NULL) ? lcu->right : ilr_lcu;
            }
            if (state->ILR_state->tile->wf_jobs[ilr_lcu->id] != NULL)
            {
              kvz_threadqueue_job_dep_add(job[0], state->ILR_state->tile->wf_jobs[ilr_lcu->id]);
            }
          } else if(state->layer != NULL) {
            if (state->layer->image_ver_scaling_jobs[lcu->id] != NULL) {
              kvz_threadqueue_job_dep_add(job[0], state->layer->image_ver_scaling_jobs[lcu->id]);
            }
            if (state->layer->cua_scaling_jobs[lcu->id] != NULL) {
              kvz_threadqueue_job_dep_add(job[0], state->layer->cua_scaling_jobs[lcu->id]);
            }
          }
        } 

        
        ////should be enough to add it to the first only?
        //if (i == 0) {
        //  encoder_state_t* parent = NULL;
        //  for (parent = state->parent; parent->parent != NULL; parent = parent->parent);

        //  /*if (parent->tqj_ilr_rec_scaling_done != NULL) {
        //    kvz_threadqueue_job_dep_add(job[0], parent->tqj_ilr_rec_scaling_done);
        //  }*/
        //  //Add dependency to ilr cua upsampling
        //  if (parent->tqj_ilr_cua_upsampling_done != NULL) {
        //    kvz_threadqueue_job_dep_add(job[0], parent->tqj_ilr_cua_upsampling_done);
        //  }

        //}
      
        //*********************************************

        kvz_threadqueue_submit(state->encoder_control->threadqueue, state->tile->wf_jobs[lcu->id]);

        // The wavefront row is done when the last LCU in the row is done.
        if (i + 1 == state->lcu_order_count) {
          assert(!state->tqj_recon_done);
          state->tqj_recon_done =
            kvz_threadqueue_copy_ref(state->tile->wf_jobs[lcu->id]);

          //*********************************************
          //For scalable extension.
          //Set scaling job dones for good measure
          if (state->layer != NULL) {
            assert(!state->tqj_ilr_rec_scaling_done);
            assert(!state->tqj_ilr_cua_upsampling_done);
            if (state->layer->image_ver_scaling_jobs[lcu->id] != NULL) {
              state->tqj_ilr_rec_scaling_done = kvz_threadqueue_copy_ref(state->layer->image_ver_scaling_jobs[lcu->id]);
            }
            if (state->layer->cua_scaling_jobs[lcu->id] != NULL) {
              state->tqj_ilr_cua_upsampling_done = kvz_threadqueue_copy_ref(state->layer->cua_scaling_jobs[lcu->id]);
            }
          }
          //*********************************************
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

static int encoder_state_tree_is_a_chain(const encoder_state_t * const state) {
  if (!state->children[0].encoder_control) return 1;
  if (state->children[1].encoder_control) return 0;
  return encoder_state_tree_is_a_chain(&state->children[0]);
}

static void encoder_state_encode(encoder_state_t * const main_state) {
  //If we have children, encode at child level
  if (main_state->children[0].encoder_control) {
    //If we have only one child, than it cannot be the last split in tree
    int node_is_the_last_split_in_tree = (main_state->children[1].encoder_control != 0);

    for (int i = 0; main_state->children[i].encoder_control; ++i) {
      encoder_state_t *sub_state = &(main_state->children[i]);

      if (sub_state->tile != main_state->tile) {
        const int offset_x = sub_state->tile->offset_x;
        const int offset_y = sub_state->tile->offset_y;
        const int width = MIN(sub_state->tile->frame->width_in_lcu * LCU_WIDTH, main_state->tile->frame->width - offset_x);
        const int height = MIN(sub_state->tile->frame->height_in_lcu * LCU_WIDTH, main_state->tile->frame->height - offset_y);

        kvz_image_free(sub_state->tile->frame->source);
        sub_state->tile->frame->source = NULL;

        kvz_image_free(sub_state->tile->frame->rec);
        sub_state->tile->frame->rec = NULL;

        kvz_cu_array_free(&sub_state->tile->frame->cu_array);

        sub_state->tile->frame->source = kvz_image_make_subimage(
            main_state->tile->frame->source,
            offset_x,
            offset_y,
            width,
            height
        );
        sub_state->tile->frame->rec = kvz_image_make_subimage(
            main_state->tile->frame->rec,
            offset_x,
            offset_y,
            width,
            height
        );
        sub_state->tile->frame->cu_array = kvz_cu_subarray(
            main_state->tile->frame->cu_array,
            offset_x,
            offset_y,
            sub_state->tile->frame->width_in_lcu * LCU_WIDTH,
            sub_state->tile->frame->height_in_lcu * LCU_WIDTH
        );
      }

      //*********************************************
      //For scalable extension.
      //Propagate layer scaling parameters to children
      //TODO: Could use subimage/array for out/in(?) images/cua?
      if (main_state->layer != NULL && sub_state->layer != main_state->layer) {
        if (!main_state->layer->scaling_started) {
          //kvz_copy_image_scaling_parameters(&sub_state->layer->img_job_param, &main_state->layer->img_job_param);
          kvz_image_free(sub_state->layer->img_job_param.pic_in);
          sub_state->layer->img_job_param.pic_in = kvz_image_copy_ref(main_state->layer->img_job_param.pic_in);
          kvz_image_free(sub_state->layer->img_job_param.pic_out);
          sub_state->layer->img_job_param.pic_out = kvz_image_copy_ref(main_state->layer->img_job_param.pic_out);
          sub_state->layer->img_job_param.param = main_state->layer->img_job_param.param;
          kvz_copy_cua_upsampling_parameters(&sub_state->layer->cua_job_param, &main_state->layer->cua_job_param);
        }
        sub_state->layer->scaling_started = main_state->layer->scaling_started;
      }
      //*********************************************

      //To be the last split, we require that every child is a chain
      node_is_the_last_split_in_tree =
        node_is_the_last_split_in_tree &&
        encoder_state_tree_is_a_chain(&main_state->children[i]);
    }
    //If it's the latest split point
    if (node_is_the_last_split_in_tree) {
      for (int i = 0; main_state->children[i].encoder_control; ++i) {
        //If we don't have wavefronts, parallelize encoding of children.
        if (main_state->children[i].type != ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
          kvz_threadqueue_free_job(&main_state->children[i].tqj_recon_done);
          main_state->children[i].tqj_recon_done =
            kvz_threadqueue_job_create(encoder_state_worker_encode_children, &main_state->children[i]);
          if (main_state->children[i].previous_encoder_state != &main_state->children[i] &&
              main_state->children[i].previous_encoder_state->tqj_recon_done &&
              !main_state->children[i].frame->is_irap)
          {
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

          //*********************************************
          //For scalable extension.
          //Set up scaling jobs
          if (main_state->children[i].layer != NULL && !main_state->children[i].layer->scaling_started) {
            start_block_step_scaling_job(&main_state->children[i], NULL);
            start_cua_lcu_scaling_job(&main_state->children[i], NULL);
            main_state->children[i].layer->scaling_started = 1; //Signal scaling started. Each tile should have its own layer struct so does not interfere with other tiles
          }

          //Add dependency to ilr recon upscaling and cua upsampling
          if (main_state->children[i].tqj_ilr_rec_scaling_done != NULL) {
            kvz_threadqueue_job_dep_add(main_state->children[i].tqj_recon_done, main_state->children[i].tqj_ilr_rec_scaling_done);
          }
          if (main_state->children[i].tqj_ilr_cua_upsampling_done != NULL) {
            kvz_threadqueue_job_dep_add(main_state->children[i].tqj_recon_done, main_state->children[i].tqj_ilr_cua_upsampling_done);
          }
          //*********************************************

          kvz_threadqueue_submit(main_state->encoder_control->threadqueue, main_state->children[i].tqj_recon_done);
        } else {
          //Wavefront rows have parallelism at LCU level, so we should not launch multiple threads here!
          //FIXME: add an assert: we can only have wavefront children
          encoder_state_worker_encode_children(&(main_state->children[i]));
        }
      }
    } else {
      for (int i = 0; main_state->children[i].encoder_control; ++i) {
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


static void encoder_ref_insertion_sort(const encoder_state_t *const state,
                                       uint8_t reflist[16],
                                       uint8_t length,
                                       bool reverse)
{

  for (uint8_t i = 1; i < length; ++i) {
    const uint8_t cur_idx = reflist[i];
    const int32_t cur_poc = state->frame->ref->pocs[cur_idx];
    int8_t j = i;
    while ((j > 0 && !reverse && cur_poc > state->frame->ref->pocs[reflist[j - 1]]) ||
           (j > 0 &&  reverse && cur_poc < state->frame->ref->pocs[reflist[j - 1]]))
    {
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
  const kvz_config *cfg = &state->encoder_control->cfg;

  FILL_ARRAY(state->frame->ref_LX_size, 0, 2);

  int num_negative = 0;
  int num_positive = 0;

  // Add positive references to L1 list
  for (int i = 0; i < state->frame->ref->used_size; i++) {
    if (state->frame->ref->pocs[i] > state->frame->poc) {
      state->frame->ref_LX[1][state->frame->ref_LX_size[1]] = i;
      state->frame->ref_LX_size[1] += 1;
      num_positive++;
    }
    
  }

  // Add negative references to L1 list when bipred is enabled and GOP is
  // either disabled or does not use picture reordering.
  bool l1_negative_refs =
    (cfg->bipred && (cfg->gop_len == 0 || cfg->gop_lowdelay));

  // Add negative references to L0 and L1 lists.
  for (int i = 0; i < state->frame->ref->used_size; i++) {
  // Modified for SHVC. TODO: Does <= really help?
  // ***********************************************
    if (state->frame->ref->pocs[i] <= state->frame->poc) {
      state->frame->ref_LX[0][state->frame->ref_LX_size[0]] = i;
      state->frame->ref_LX_size[0] += 1;
      if (l1_negative_refs) {
        state->frame->ref_LX[1][state->frame->ref_LX_size[1]] = i;
        state->frame->ref_LX_size[1] += 1;
      }
      num_negative++;
    }
  // ***********************************************
  }

  // Fill the rest with -1.
  for (int i = state->frame->ref_LX_size[0]; i < 16; i++) {
    state->frame->ref_LX[0][i] = 0xff;
  }
  for (int i = state->frame->ref_LX_size[1]; i < 16; i++) {
    state->frame->ref_LX[1][i] = 0xff;
  }

  // Sort reference lists.
  encoder_ref_insertion_sort(state, state->frame->ref_LX[0], num_negative, false);
  encoder_ref_insertion_sort(state, state->frame->ref_LX[1], num_positive, true);
  if (l1_negative_refs) {
    encoder_ref_insertion_sort(state, state->frame->ref_LX[1] + num_positive, num_negative, false);
  }
}

/**
 * \brief Remove any references that should no longer be used.
 */
static void encoder_state_remove_refs(encoder_state_t *state) {
  const encoder_control_t * const encoder = state->encoder_control;

  //*********************************************
  //For scalable extension.
  kvz_image_list_rem_ILR(state->frame->ref,
                         state->frame->poc,
                         encoder->cfg.gop[state->frame->gop_offset].tId,
                         state->encoder_control->layer.layer_id );
  //*********************************************

  int neg_refs = encoder->cfg.gop[state->frame->gop_offset].ref_neg_count;
  int pos_refs = encoder->cfg.gop[state->frame->gop_offset].ref_pos_count;

  unsigned target_ref_num;
  if (encoder->cfg.gop_len) {
    target_ref_num = neg_refs + pos_refs;
  } else {
    target_ref_num = encoder->cfg.ref_frames;
  }
  
  if (state->frame->pictype == KVZ_NAL_IDR_W_RADL ||
      state->frame->pictype == KVZ_NAL_IDR_N_LP)
  {
    target_ref_num = 0;
  }

  //*********************************************
  //For scalable extension.
  //Add space for irl to the list
  target_ref_num += encoder->cfg.ILR_frames;
  //*********************************************

  if (encoder->cfg.gop_len && target_ref_num > 0) {
    // With GOP in use, go through all the existing reference pictures and
    // remove any picture that is not referenced by the current picture.

    for (int ref = state->frame->ref->used_size - 1; ref >= 0; --ref) {
      //*********************************************
      //For scalable extension.
      //If ref is ILR no need to remove. (invalid ILR refs should already be removed)
      if( state->frame->ref->image_info[ref].layer_id < state->encoder_control->layer.layer_id){
        continue;
      }
      //*********************************************

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

      if (ref_poc < state->frame->irap_poc &&
          state->frame->irap_poc < state->frame->poc)
      {
        // Trailing frames cannot refer to leading frames.
        is_referenced = false;
      }

      if (encoder->cfg.intra_period > 0 &&
          ref_poc < state->frame->irap_poc - encoder->cfg.intra_period)
      {
        // No frame can refer past the two preceding IRAP frames.
        is_referenced = false;
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
  //*********************************************
  //For scalable extension.
  kvz_threadqueue_free_job(&state->tqj_ilr_rec_scaling_done);
  kvz_threadqueue_free_job(&state->tqj_ilr_cua_upsampling_done);  
  //*********************************************

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
  
  assert(!state->tile->frame->cu_array);
  state->tile->frame->cu_array = kvz_cu_array_alloc(
      state->tile->frame->width,
      state->tile->frame->height
  );

  // Set POC.
  if (state->frame->num == 0) {
    state->frame->poc = 0;
  } else if (cfg->gop_len && !cfg->gop_lowdelay) {
    // Calculate POC according to the global frame counter and GOP structure
    int32_t poc = state->frame->num - 1;
    int32_t poc_offset = cfg->gop[state->frame->gop_offset].poc_offset;
    state->frame->poc = poc - poc % cfg->gop_len + poc_offset;
    kvz_videoframe_set_poc(state->tile->frame, state->frame->poc);
  } else if (cfg->intra_period > 0) {
    state->frame->poc = state->frame->num % cfg->intra_period;
  } else {
    state->frame->poc = state->frame->num;
  }

  // Check whether the frame is a keyframe or not.
  if (state->frame->num == 0) {
    state->frame->is_irap = true;
  } else {
    state->frame->is_irap =
      cfg->intra_period > 0 &&
      (state->frame->poc % cfg->intra_period) == 0;
  }
 
  if (state->frame->is_irap) {
    state->frame->irap_poc = state->frame->poc;
  }

  // Set pictype.
  if (state->frame->is_irap) {
    if (state->frame->num == 0 ||
        cfg->intra_period == 1 ||
        cfg->gop_len == 0 ||
        cfg->gop_lowdelay)
    {
      state->frame->pictype = KVZ_NAL_IDR_W_RADL;
    } else {
      state->frame->pictype = KVZ_NAL_CRA_NUT;
    }
  } else if (state->frame->poc < state->frame->irap_poc) {
    state->frame->pictype = KVZ_NAL_RASL_R;
  } else {
    state->frame->pictype = KVZ_NAL_TRAIL_R;
  }

  encoder_state_remove_refs(state);
  kvz_encoder_create_ref_lists(state);

  //*********************************************
  //For scalable extension. TODO: Enable encoding el frames with intra based on intra period?
  
  // Set slicetype.
  if (state->frame->is_irap && state->encoder_control->layer.layer_id == 0) {
    state->frame->slicetype = KVZ_SLICE_I;
  } else if (state->frame->ref_LX_size[1] > 0) {
    state->frame->slicetype = KVZ_SLICE_B;
  } else {
    state->frame->slicetype = KVZ_SLICE_P;
  }
  
  //*********************************************

  if (cfg->target_bitrate > 0 && state->frame->num > cfg->owf) {
    normalize_lcu_weights(state);
  }
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
  // ***********************************************
  // Modified for SHVC.
  // Do scalability preparation here so that the ilr rec is set when using gop
  scalability_prepare(state);
  prepare_ilr_frames(state);
  // ***********************************************

  encoder_state_init_new_frame(state, frame);
  // ***********************************************
  // Modified for SHVC.
  encoder_state_set_rps(state);
  // ***********************************************
  encoder_state_encode(state);

  threadqueue_job_t *job =
    kvz_threadqueue_job_create(kvz_encoder_state_worker_write_bitstream, state);

  _encode_one_frame_add_bitstream_deps(state, job);
  if (state->previous_encoder_state != state && state->previous_encoder_state->tqj_bitstream_written) {
    //We need to depend on previous bitstream generation
    kvz_threadqueue_job_dep_add(job, state->previous_encoder_state->tqj_bitstream_written);
  }
  kvz_threadqueue_submit(state->encoder_control->threadqueue, job);
  assert(!state->tqj_bitstream_written);
  state->tqj_bitstream_written = job;

  state->frame->done = 0;
}


// ***********************************************
// Modified for SHVC.

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
    state->frame->poc = 0;
    state->frame->irap_poc = 0;
    assert(!state->tile->frame->source);
    assert(!state->tile->frame->rec);
    assert(!state->tile->frame->cu_array);
    //scalability_prepare(state);
    //add_ilr_frames(state);//Need to do this here so that the first el frame can be inter

    state->frame->prepared = 1;
    return;
  }
  
  // NOTE: prev_state is equal to state when OWF is zero
  encoder_state_t *prev_state = state->previous_encoder_state;

  if (state->previous_encoder_state != state) {
    kvz_cu_array_free(&state->tile->frame->cu_array);
    unsigned width  = state->tile->frame->width_in_lcu  * LCU_WIDTH;
    unsigned height = state->tile->frame->height_in_lcu * LCU_WIDTH;
    state->tile->frame->cu_array = kvz_cu_array_alloc(width, height);

    kvz_image_list_copy_contents(state->frame->ref, prev_state->frame->ref);
    kvz_encoder_create_ref_lists(state);
  }

  //TODO: remove
  // For SHVC.
  //if (encoder->layer.layer_id > 0 && prev_state->ILR_state != NULL) {
  //    kvz_image_list_rem_ILR(state->frame->ref, prev_state->frame->poc); //Remove old ILR pics from the ref list so they don't interfere.
  //}
  //scalability_prepare(state);

  if (!encoder->cfg.gop_len ||
      !prev_state->frame->poc ||
      encoder->cfg.gop[prev_state->frame->gop_offset].is_ref) {

    // Store current list of POCs for use in TMVP derivation
    memcpy(prev_state->tile->frame->rec->ref_pocs, state->frame->ref->pocs, sizeof(int32_t)*state->frame->ref->used_size);
    // For SHVC.
    // Also store image info
    memcpy(prev_state->tile->frame->rec->picture_info, state->frame->ref->image_info, sizeof(kvz_picture_info_t) * state->frame->ref->used_size);
    // For SHVC.
    // Add previous reconstructed picture as a reference
    kvz_image_list_add(state->frame->ref,
                   prev_state->tile->frame->rec,
                   prev_state->tile->frame->cu_array,
                   prev_state->frame->poc,
                   prev_state->frame->ref_LX,
                   prev_state->encoder_control->cfg.gop[prev_state->frame->gop_offset].tId,
                   prev_state->encoder_control->layer.layer_id,
                   0); //Currently only ILR can be a long term references
    kvz_cu_array_free(&state->tile->frame->cu_array);
    unsigned height = state->tile->frame->height_in_lcu * LCU_WIDTH;
    unsigned width  = state->tile->frame->width_in_lcu  * LCU_WIDTH;
    state->tile->frame->cu_array = kvz_cu_array_alloc(width, height);
  }

  //add_ilr_frames(state);

  // Remove source and reconstructed picture.
  kvz_image_free(state->tile->frame->source);
  state->tile->frame->source = NULL;

  kvz_image_free(state->tile->frame->rec);
  state->tile->frame->rec = NULL;

  kvz_cu_array_free(&state->tile->frame->cu_array);

  // Update POC and frame count.
  state->frame->num = prev_state->frame->num + 1;
  state->frame->poc = prev_state->frame->poc + 1;
  state->frame->irap_poc = prev_state->frame->irap_poc;

  state->frame->prepared = 1;
}

// ***********************************************

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

int kvz_get_cu_ref_qp(const encoder_state_t *state, int x, int y, int last_qp)
{
  const encoder_control_t *ctrl = state->encoder_control;
  const cu_array_t *cua = state->tile->frame->cu_array;
  // Quantization group width
  const int qg_width = LCU_WIDTH >> MIN(ctrl->max_qp_delta_depth, kvz_cu_array_at_const(cua, x, y)->depth);

  // Coordinates of the top-left corner of the quantization group
  const int x_qg = x & ~(qg_width - 1);
  const int y_qg = y & ~(qg_width - 1);

  int qp_pred_a = last_qp;
  if (x_qg % LCU_WIDTH > 0) {
    qp_pred_a = kvz_cu_array_at_const(cua, x_qg - 1, y_qg)->qp;
  }

  int qp_pred_b = last_qp;
  if (y_qg % LCU_WIDTH > 0) {
    qp_pred_b = kvz_cu_array_at_const(cua, x_qg, y_qg - 1)->qp;
  }

  return ((qp_pred_a + qp_pred_b + 1) >> 1);
}
