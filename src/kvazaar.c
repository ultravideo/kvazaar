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

#include "kvazaar.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitstream.h"
#include "cfg.h"
#include "checkpoint.h"
#include "encoder.h"
#include "encoder_state-bitstream.h"
#include "encoder_state-ctors_dtors.h"
#include "encoderstate.h"
#include "global.h"
#include "image.h"
#include "input_frame_buffer.h"
#include "kvazaar_internal.h"
#include "strategyselector.h"
#include "threadqueue.h"
#include "videoframe.h"

// ***********************************************
  // Modified for SHVC
#include "scaler/scaler.h"

// ***********************************************


static void kvazaar_close(kvz_encoder *encoder)
{
  // ***********************************************
  // Modified for SHVC. TODO: Account for more complex ref structures?
  kvz_encoder *next = NULL;
  if (encoder) {
    next = encoder->next_enc;
    encoder->next_enc = NULL;
    encoder->prev_enc = NULL; //Prev encoder should be closed in a previous call

    // The threadqueue must be stopped before freeing states.
    if (encoder->control) {
      kvz_threadqueue_stop(encoder->control->threadqueue);
    }

    if (encoder->states) {
      // Flush input frame buffer.
      kvz_picture *pic = NULL;
      while ((pic = kvz_encoder_feed_frame(&encoder->input_buffer,
                                           &encoder->states[0],
                                           NULL)) != NULL) {
        kvz_image_free(pic);
        pic = NULL;
      }

      for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
        kvz_encoder_state_finalize(&encoder->states[i]);
      }
    }
    FREE_POINTER(encoder->states);

    //Threadqueue is shared so free it only in the last encoder and set to NULL for others
    if(next!=NULL) {
      encoder_control_t *ctrl = (encoder_control_t*)encoder->control;
      ctrl->threadqueue = NULL;
    }

    // Discard const from the pointer.
    kvz_encoder_control_free((void*) encoder->control);
    encoder->control = NULL;
  }
  else {
    return;
  }
  FREE_POINTER(encoder);
 
  kvazaar_close(next);
  // ***********************************************
}

// ***********************************************
// Modified for SHVC.

static void print_encoderstate_hierarchy( encoder_state_t* state, int indent )
{
  //Formatting definitions
  static const char *ind_str = "  "; //String used for indentation
  static const char *hor_con = "|"; //Horisontal connector
  static const char *ver_con = "-"; //Vertical connector

  //Print cur state
  for (int i = indent; i > 0; i--) {
    fputs(ind_str, stderr);
    //fputs(hor_con, stderr);
    if( i == 1 ){
      fputs(hor_con, stderr);
      fputs(ver_con, stderr);
    }
  }

  switch( state->type ){
  
  case ENCODER_STATE_TYPE_INVALID:
    fputs("ENCODER_STATE_TYPE_INVLID\n", stderr);
    break;

  case ENCODER_STATE_TYPE_MAIN:
    fputs("ENCODER_STATE_TYPE_MAIN\n", stderr);
    break;

  case ENCODER_STATE_TYPE_SLICE:
    fputs("ENCODER_STATE_TYPE_SLICE\n", stderr);
    break;

  case ENCODER_STATE_TYPE_TILE:
    fputs("ENCODER_STATE_TYPE_TILE\n", stderr);
    break;

  case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
    fputs("ENCODER_STATE_TYPE_WAVEFRONT_ROW\n", stderr);
    break;

  default:
    //Nothing to do
    break;
  }

  //Print hierarchy tree in a depth first manner
  for (int i = 0; state->children[i].encoder_control; i++) {
    print_encoderstate_hierarchy(&state->children[i], indent + 1);
  }

}
// ***********************************************

static kvz_encoder * kvazaar_open(const kvz_config *cfg)
{
  kvz_encoder *encoder = NULL; //The base layer encoder

  //Initialize strategies
  // TODO: Make strategies non-global
  if (!kvz_strategyselector_init(cfg->cpuid, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "Failed to initialize strategies.\n");
    goto kvazaar_open_failure;
  }

  // ***********************************************
  // Modified for SHVC. TODO: Account for more complex ref structures?
  kvz_encoder *cur_enc = NULL;
  kvz_encoder *prev_enc = NULL;
  const encoder_control_t *ctrl = kvz_encoder_control_init(cfg); //Initializes the control structures for different layers/encoders
  if (ctrl == NULL) {
    goto kvazaar_open_failure;
  }

  //TODO: Just use a while loop?
  for (; ctrl != NULL ; ctrl = ctrl->next_enc_ctrl) {
    
    cur_enc = calloc(1, sizeof(kvz_encoder));
    if (!cur_enc) {
      goto kvazaar_open_failure;
    }
    if( encoder == NULL ) encoder = cur_enc;

    //Connect the sub encoders
    cur_enc->next_enc = NULL;
    cur_enc->prev_enc = prev_enc;
    if (prev_enc != NULL) prev_enc->next_enc = cur_enc;

    
    cur_enc->control = ctrl;

    cur_enc->num_encoder_states = encoder->control->cfg.owf + 1;
    cur_enc->cur_state_num = 0;
    cur_enc->out_state_num = 0;
    cur_enc->frames_started = 0;
    cur_enc->frames_done = 0;

    kvz_init_input_frame_buffer(&cur_enc->input_buffer);

    cur_enc->states = calloc(cur_enc->num_encoder_states, sizeof(encoder_state_t));
    if (!cur_enc->states) {
      goto kvazaar_open_failure;
    }

    for (unsigned i = 0; i < cur_enc->num_encoder_states; ++i) {
      cur_enc->states[i].encoder_control = cur_enc->control;

      if (!kvz_encoder_state_init(&cur_enc->states[i], NULL)) {
        goto kvazaar_open_failure;
      }

      cur_enc->states[i].frame->QP = (int8_t)cfg->qp;
    }

    for (int i = 0; i < cur_enc->num_encoder_states; ++i) {
      if (i == 0) {
        cur_enc->states[i].previous_encoder_state = &cur_enc->states[cur_enc->num_encoder_states - 1];
      } else {
        cur_enc->states[i].previous_encoder_state = &cur_enc->states[(i - 1) % cur_enc->num_encoder_states];
      }
      kvz_encoder_state_match_children_of_previous_frame(&cur_enc->states[i]);
    }

    //Set ILR states for the current encoder's states. TODO: Account for a more complex ref structure
    //TODO: error checking
    if (prev_enc != NULL) {
      for (int i = 0; i < cur_enc->num_encoder_states; ++i) {
        //Should give the state that encodes the ILR frame
        cur_enc->states[i].ILR_state = &prev_enc->states[i];
        cur_enc->states[i].num_ILR_states = 1;
        kvz_encoder_state_match_ILR_states_of_children(&cur_enc->states[i]);
      }
    }

    cur_enc->states[cur_enc->cur_state_num].frame->num = -1;
    
    //Print encoder state hierarchy for debugging purposes.
    if (cfg->shared != NULL && cfg->shared->print_es_hierarchy) {
      fprintf(stderr, "Layer %d encoder state hierarchy:\n", cfg->layer);
      print_encoderstate_hierarchy(&cur_enc->states[cur_enc->cur_state_num], 0);
    }

    //Prepare for the next loop
    prev_enc = cur_enc;
    cfg = cfg->next_cfg;

  }


  // ***********************************************
  return encoder;

kvazaar_open_failure:
  kvazaar_close(encoder);
  return NULL;
}


static void set_frame_info(kvz_frame_info *const info, const encoder_state_t *const state)
{
  info->poc = state->frame->poc,
  info->qp = state->frame->QP;
  info->nal_unit_type = state->frame->pictype;
  info->slice_type = state->frame->slicetype;

  memset(info->ref_list[0], 0, 16 * sizeof(int));
  memset(info->ref_list[1], 0, 16 * sizeof(int));

  // ***********************************************
  // Modified for SHVC.
  
  for (size_t i = 0; i < state->frame->ref_LX_size[0]; i++) {
    if (!state->local_rps->is_used[state->frame->ref_LX[0][i]]) continue;

    info->ref_list[0][i] = state->frame->ref->pocs[state->frame->ref_LX[0][i]];
    info->ref_list_len[0]++;
  }

  for (size_t i = 0; i < state->frame->ref_LX_size[1]; i++) {
    if(!state->local_rps->is_used[state->frame->ref_LX[1][i]]) continue;

    info->ref_list[1][i] = state->frame->ref->pocs[state->frame->ref_LX[1][i]];
    info->ref_list_len[1]++;
  }

  info->ref_list_len[0] = state->local_rps->num_ref_idx_LX_active[0]; //state->frame->ref_LX_size[0];
  info->ref_list_len[1] = state->local_rps->num_ref_idx_LX_active[1]; //state->frame->ref_LX_size[1];

  info->lid = state->encoder_control->layer.layer_id;
  info->tid = state->encoder_control->cfg.gop[state->frame->gop_offset].tId;
  // ***********************************************
}


static int kvazaar_headers(kvz_encoder *enc,
                           kvz_data_chunk **data_out,
                           uint32_t *len_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) *len_out = 0;

  bitstream_t stream;
  kvz_bitstream_init(&stream);

  kvz_encoder_state_write_parameter_sets(&stream, &enc->states[enc->cur_state_num]);

  // Get stream length before taking chunks since that clears the stream.
  if (len_out) *len_out = kvz_bitstream_tell(&stream) / 8;
  if (data_out) *data_out = kvz_bitstream_take_chunks(&stream);

  kvz_bitstream_finalize(&stream);
  return 1;
}


/**
* \brief Separate a single field from a frame.
*
* \param frame_in           input frame to extract field from
* \param source_scan_type   scan type of input material (0: progressive, 1:top field first, 2:bottom field first)
* \param field parity   
* \param field_out
*
* \return              1 on success, 0 on failure
*/
static int yuv_io_extract_field(const kvz_picture *frame_in, unsigned source_scan_type, unsigned field_parity, kvz_picture *field_out)
{
  if ((source_scan_type != 1) && (source_scan_type != 2)) return 0;
  if ((field_parity != 0)     && (field_parity != 1))     return 0;

  unsigned offset = 0;
  if (source_scan_type == 1) offset = field_parity ? 1 : 0;
  else if (source_scan_type == 2) offset = field_parity ? 0 : 1;  

  //Luma
  for (int i = 0; i < field_out->height; ++i){
    kvz_pixel *row_in  = frame_in->y + MIN(frame_in->height - 1, 2 * i + offset) * frame_in->stride;
    kvz_pixel *row_out = field_out->y + i * field_out->stride;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width);
  }

  //Chroma
  for (int i = 0; i < field_out->height / 2; ++i){
    kvz_pixel *row_in = frame_in->u + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    kvz_pixel *row_out = field_out->u + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width / 2);
  }

  for (int i = 0; i < field_out->height / 2; ++i){
    kvz_pixel *row_in = frame_in->v + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    kvz_pixel *row_out = field_out->v + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width / 2);
  }

  return 1;
}


static int kvazaar_encode(kvz_encoder *enc,
                          kvz_picture *pic_in,
                          kvz_data_chunk **data_out,
                          uint32_t *len_out,
                          kvz_picture **pic_out,
                          kvz_picture **src_out,
                          kvz_frame_info *info_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) *len_out = 0;
  if (pic_out) *pic_out = NULL;
  if (src_out) *src_out = NULL;

  encoder_state_t *state = &enc->states[enc->cur_state_num];

  if (!state->frame->prepared) {
    
    kvz_encoder_prepare(state);

  }

  if (pic_in != NULL) {
    // FIXME: The frame number printed here is wrong when GOP is enabled.
    CHECKPOINT_MARK("read source frame: %d", state->frame->num + enc->control->cfg.seek);
  }

  kvz_picture* frame = kvz_encoder_feed_frame(&enc->input_buffer, state, pic_in);
  if (frame) {
    assert(state->frame->num == enc->frames_started);
    // Start encoding.
    kvz_encode_one_frame(state, frame);
    enc->frames_started += 1;
  }

  // If we have finished encoding as many frames as we have started, we are done.
  if (enc->frames_done == enc->frames_started) {
    return 1;
  }

  if (!state->frame->done) {
    // We started encoding a frame; move to the next encoder state.
    enc->cur_state_num = (enc->cur_state_num + 1) % (enc->num_encoder_states);
  }

  encoder_state_t *output_state = &enc->states[enc->out_state_num];
  if (!output_state->frame->done &&
      (pic_in == NULL || enc->cur_state_num == enc->out_state_num)) {

    kvz_threadqueue_waitfor(enc->control->threadqueue, output_state->tqj_bitstream_written);
    // The job pointer must be set to NULL here since it won't be usable after
    // the next frame is done.
    kvz_threadqueue_free_job(&output_state->tqj_bitstream_written);

    // Get stream length before taking chunks since that clears the stream.
    if (len_out) *len_out = kvz_bitstream_tell(&output_state->stream) / 8;
    if (data_out) *data_out = kvz_bitstream_take_chunks(&output_state->stream);
    if (pic_out) *pic_out = kvz_image_copy_ref(output_state->tile->frame->rec);
    if (src_out) *src_out = kvz_image_copy_ref(output_state->tile->frame->source);
    if (info_out) set_frame_info(info_out, output_state);

    output_state->frame->done = 1;
    output_state->frame->prepared = 0;
    enc->frames_done += 1;

    enc->out_state_num = (enc->out_state_num + 1) % (enc->num_encoder_states);
  }

  return 1;
}


//An unit delay to the encoding process
/*static void encode_delay(const int layer, kvz_picture **pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
{
  assert(layer < MAX_LAYERS);

  //Store values in static variables
  static kvz_picture* d_pic_in[MAX_LAYERS] = { NULL };
  static kvz_data_chunk* d_data_out[MAX_LAYERS] = { NULL };
  static uint32_t d_len_out[MAX_LAYERS] = { 0 };
  static kvz_picture* d_pic_out[MAX_LAYERS] = { NULL };
  static kvz_picture* d_src_out[MAX_LAYERS] = { NULL };
  static kvz_frame_info d_info_out[MAX_LAYERS] = { { 0 } };

  //Store new value and return the old one
  void *tmp = (void *)*pic_in;
  *pic_in = d_pic_in[layer];
  d_pic_in[layer] = (kvz_picture *)tmp;

  tmp = (void *)*data_out;
  *data_out = d_data_out[layer];
  d_data_out[layer] = (kvz_data_chunk *)tmp;

  tmp = (void *)*pic_out;
  *pic_out = d_pic_out[layer];
  d_pic_out[layer] = (kvz_picture *)tmp;

  tmp = (void *)*src_out;
  *src_out = d_src_out[layer];
  d_src_out[layer] = (kvz_picture *)tmp;

  uint32_t tmp_len = *len_out;
  *len_out = d_len_out[layer];
  d_len_out[layer] = tmp_len;

  kvz_frame_info tmp_info = *info_out;
  *info_out = d_info_out[layer];
  d_info_out[layer] = tmp_info;
}*/

//TODO: make a note of this: Asume that info_out is an array with an element for each layer
//TODO: Allow scaling "step-wise" instead of allways from the original, for a potentially reduced complexity?
//TODO: Account for pic_in containing several input images for different layers
//Use this function to aggregate the results etc. but otherwise just call kvazaar_encode with the correct encoder
/*static int kvazaar_scalable_encode(kvz_encoder* enc, kvz_picture* pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) memset(len_out, 0, sizeof(uint32_t)*enc->control->layer.max_layers);
  if (pic_out) *pic_out = NULL;
  if (src_out) *src_out = NULL;

  //Pic_in should contain the input images chained using base_image.
  //Move them to a list for easier access
  kvz_picture **pics_in = calloc(enc->control->layer.max_layers, sizeof(kvz_picture*));
  if (pics_in == NULL) {
    fprintf(stderr, "Memory error: Could not allocate picture array.\n");
    return 0;
  }
  pics_in[0] = pic_in;

  if (pic_in != NULL) {
    for (int i = 1; pic_in->base_image != pic_in; i++) {
      pics_in[i] = pic_in->base_image;
      pic_in = pic_in->base_image;
    }
  }
  
  kvz_encoder *cur_enc = enc;

  int el_tmvp_enabled = false;// enc->next_enc->control->cfg.tmvp_enable && enc->control->cfg.owf > 1; //TODO: does not work when owf == 1; find a work-around

  //Pre prepare statest to prevent data-races between ILR states (when copying stuff needed for tmvp)
  while (cur_enc != NULL) {

    encoder_state_t *state = &cur_enc->states[cur_enc->cur_state_num];

    if (!state->frame->prepared) {
      kvz_encoder_prepare(state);
    }

    cur_enc = cur_enc->next_enc;
  }

  cur_enc = enc;

  //TODO: Use a while loop instead?
  //for( unsigned i = 0; i < enc->control->layer.max_layers; i++) {
  for( unsigned i = 0; cur_enc != NULL; i++) {  
    
    //Use these to pass stuff to the actual encoder function and aggregate the results into the actual parameters
    kvz_picture *cur_pic_in = kvz_image_scaling(pics_in[cur_enc->control->layer.input_layer], &cur_enc->control->layer.downscaling, 1);

    kvz_data_chunk* cur_data_out = NULL;
    uint32_t cur_len_out = 0;
    kvz_picture *cur_pic_out = NULL;
    kvz_picture *cur_src_out = NULL;


    if (el_tmvp_enabled && i > 0) {
      //TODO: Account for multiple EL layers needing to be delayed (when EL references another EL)
      encode_delay(i, &cur_pic_in, &cur_data_out, &cur_len_out, &cur_pic_out, &cur_src_out, &(info_out[i]));
    }

    if(!kvazaar_encode(cur_enc, cur_pic_in, &cur_data_out, &cur_len_out, &cur_pic_out, &cur_src_out, &(info_out[i]))) {
      kvz_image_free(cur_pic_in);
      free(pics_in);
      return 0;
    }

    if (el_tmvp_enabled && i == 0) {
      encode_delay(i, &cur_pic_in, &cur_data_out, &cur_len_out, &cur_pic_out, &cur_src_out, &(info_out[i]));
    }

    kvz_image_free(cur_pic_in);
    cur_pic_in = NULL;

    //Aggregate new stuff
    if (data_out) {
      if (*data_out == NULL) {
        *data_out = cur_data_out;
      } else {
        while ((*data_out)->next != NULL) {
          data_out = &(*data_out)->next;
        }
        (*data_out)->next = cur_data_out;
      }
      cur_data_out = NULL;
    }
    if (len_out) len_out[i] += cur_len_out;
    if (pic_out) {
      if (*pic_out == NULL) {
        *pic_out = cur_pic_out;
      } else {
        if (cur_pic_out != *pic_out) {
          (*pic_out)->base_image = kvz_image_copy_ref(cur_pic_out);
          pic_out = &(*pic_out)->base_image;
        }
        kvz_image_free(cur_pic_out);
      }
      cur_pic_out = NULL;
    }
    if (src_out) {
      if (*src_out == NULL) {
        *src_out = cur_src_out;
      } else {
        if (cur_src_out != *src_out) {
          (*src_out)->base_image = kvz_image_copy_ref(cur_src_out);
          src_out = &(*src_out)->base_image;
        }
        kvz_image_free(cur_src_out);
      }
      cur_src_out = NULL;
    }

    //Update other values for the next layer
    cur_enc = cur_enc->next_enc;
  }

  free(pics_in);
  
  return 1;
}*/

//TODO: make a note of this: Asume that info_out is an array with an element for each layer
//TODO: Allow scaling "step-wise" instead of allways from the original, for a potentially reduced complexity?
//Custom encoding loop for scalable encoding
static int kvazaar_scalable_encode(kvz_encoder *enc,
  kvz_picture *pic_in,
  kvz_data_chunk **data_out,
  uint32_t *len_out,
  kvz_picture **pic_out,
  kvz_picture **src_out,
  kvz_frame_info *info_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) memset(len_out, 0, sizeof(uint32_t)*enc->control->layer.max_layers);
  if (pic_out) *pic_out = NULL;
  if (src_out) *src_out = NULL;

  //Pic_in should contain the input images chained using base_image.
  //Move them to a list for easier access
  kvz_picture *pics_in[MAX_LAYERS] = { NULL };
  pics_in[0] = pic_in;
  if (pic_in != NULL) {
    for (int i = 1; pic_in->base_image != pic_in; i++) {
      pics_in[i] = pic_in->base_image;
      pic_in = pic_in->base_image;
    }
  }

  //Make a list of encoders
  kvz_encoder *enc_list[MAX_LAYERS] = { NULL };
  int num_enc = 0;
  while (enc != NULL) {
    enc_list[num_enc] = enc;
    num_enc++;
    enc = enc->next_enc;
  }

  //For keeping track of states
  int frame_initialized[MAX_LAYERS] = { 0 };
  int pic_in_is_null[MAX_LAYERS] = { 0 };

  //Prepare current states
  for (int i = 0; i < num_enc; i++)
  {
    encoder_state_t *state = &enc_list[i]->states[enc_list[i]->cur_state_num];

    if (!state->frame->prepared) {
      kvz_encoder_prepare(state);
    }

    //Use these to store intermediate values of each encoder and aggregate the results into the actual output parameters
    kvz_picture *cur_pic_in = kvz_image_scaling(pics_in[enc_list[i]->control->layer.input_layer], &enc_list[i]->control->layer.downscaling, 1);
    pic_in_is_null[i] = cur_pic_in == NULL;

    if (cur_pic_in != NULL) {
      // FIXME: The frame number printed here is wrong when GOP is enabled.
      CHECKPOINT_MARK("read source frame: %d", state->frame->num + enc_list[i]->control->cfg.seek);
    }

    kvz_picture* frame = kvz_encoder_feed_frame(&enc_list[i]->input_buffer, state, cur_pic_in);
    if (frame) {
      assert(state->frame->num == enc_list[i]->frames_started);

      kvz_scalability_prepare(state);

      kvz_init_one_frame(state, frame);
      frame_initialized[i] = true;
    }

    kvz_image_free(cur_pic_in);
  }

  //Start encoding frames
  for (int i = 0; i < num_enc; i++) {

    encoder_state_t *state = &enc_list[i]->states[enc_list[i]->cur_state_num];

    if (frame_initialized[i]) {
      kvz_start_encode_one_frame(state);
      enc_list[i]->frames_started += 1;
    }
  }

  //Get output and move to next state etc.
  for (int i = 0; i < num_enc; i++) {

    // If we have finished encoding as many frames as we have started, we are done.
    if (enc_list[i]->frames_done == enc_list[i]->frames_started) {
      continue;
    }

    encoder_state_t *state = &enc_list[i]->states[enc_list[i]->cur_state_num];

    if (!state->frame->done) {
      // We started encoding a frame; move to the next encoder state.
      enc_list[i]->cur_state_num = (enc_list[i]->cur_state_num + 1) % (enc_list[i]->num_encoder_states);
    }

    encoder_state_t *output_state = &enc_list[i]->states[enc_list[i]->out_state_num];
    if (!output_state->frame->done &&
      (pic_in_is_null[i] || enc_list[i]->cur_state_num == enc_list[i]->out_state_num)) {

      kvz_threadqueue_waitfor(enc_list[i]->control->threadqueue, output_state->tqj_bitstream_written);
      // The job pointer must be set to NULL here since it won't be usable after
      // the next frame is done.
      kvz_threadqueue_free_job(&output_state->tqj_bitstream_written);

      //Aggregate new stuff
      // Get stream length before taking chunks since that clears the stream.
      if (len_out) len_out[i] += kvz_bitstream_tell(&output_state->stream) / 8;
      if (data_out) {
        if (*data_out == NULL) {
          *data_out = kvz_bitstream_take_chunks(&output_state->stream);
        } else {
          while ((*data_out)->next != NULL) {
            data_out = &(*data_out)->next;
          }
          (*data_out)->next = kvz_bitstream_take_chunks(&output_state->stream);
        }
      }
      if (pic_out) {
        if (*pic_out == NULL) {
          *pic_out = kvz_image_copy_ref(output_state->tile->frame->rec);
        } else if (output_state->tile->frame->rec != *pic_out) {
          (*pic_out)->base_image = kvz_image_copy_ref(output_state->tile->frame->rec);
          pic_out = &(*pic_out)->base_image;
        }
      }
      if (src_out) {
        if (*src_out == NULL) {
          *src_out = kvz_image_copy_ref(output_state->tile->frame->source);
        } else if (output_state->tile->frame->source != *src_out) {
          (*src_out)->base_image = kvz_image_copy_ref(output_state->tile->frame->source);
          src_out = &(*src_out)->base_image;
        }
      }
      if (info_out) set_frame_info(&info_out[i], output_state);

      output_state->frame->done = 1;
      output_state->frame->prepared = 0;
      enc_list[i]->frames_done += 1;

      enc_list[i]->out_state_num = (enc_list[i]->out_state_num + 1) % (enc_list[i]->num_encoder_states);
    }
  }

  return 1;
}


//TODO: Handle interlaced
static int kvazaar_field_encoding_adapter(kvz_encoder *enc,
                                          kvz_picture *pic_in,
                                          kvz_data_chunk **data_out,
                                          uint32_t *len_out,
                                          kvz_picture **pic_out,
                                          kvz_picture **src_out,
                                          kvz_frame_info *info_out)
{
  if (enc->control->cfg.source_scan_type == KVZ_INTERLACING_NONE) {
    // For progressive, simply call the normal encoding function.
    //If several layers are used, call the apropriate function
    //If base layer and input layer differ in size, use scalable
    uint8_t bl_scaling = enc->control->layer.downscaling.src_width != enc->control->layer.downscaling.trgt_width ||
                         enc->control->layer.downscaling.src_height != enc->control->layer.downscaling.trgt_height;
    if(enc->control->layer.max_layers > 1 || bl_scaling) return kvazaar_scalable_encode(enc, pic_in, data_out, len_out, pic_out, src_out, info_out);
    return kvazaar_encode(enc, pic_in, data_out, len_out, pic_out, src_out, info_out);
  }

// ***********************************************

  // For interlaced, make two fields out of the input frame and call encode on them separately.
  encoder_state_t *state = &enc->states[enc->cur_state_num];
  kvz_picture *first_field = NULL, *second_field = NULL;
  struct {
    kvz_data_chunk* data_out;
    uint32_t len_out;
  } first = { 0, 0 }, second = { 0, 0 };

  if (pic_in != NULL) {
    first_field = kvz_image_alloc(state->encoder_control->chroma_format, state->encoder_control->in.width, state->encoder_control->in.height);
    if (first_field == NULL) {
      goto kvazaar_field_encoding_adapter_failure;
    }
    second_field = kvz_image_alloc(state->encoder_control->chroma_format, state->encoder_control->in.width, state->encoder_control->in.height);
    if (second_field == NULL) {
      goto kvazaar_field_encoding_adapter_failure;
    }

    yuv_io_extract_field(pic_in, pic_in->interlacing, 0, first_field);
    yuv_io_extract_field(pic_in, pic_in->interlacing, 1, second_field);
    
    first_field->pts = pic_in->pts;
    first_field->dts = pic_in->dts;
    first_field->interlacing = pic_in->interlacing;

    // Should the second field have higher pts and dts? It shouldn't affect anything.
    second_field->pts = pic_in->pts;
    second_field->dts = pic_in->dts;
    second_field->interlacing = pic_in->interlacing;
  }

  if (!kvazaar_encode(enc, first_field, &first.data_out, &first.len_out, pic_out, NULL, info_out)) {
    goto kvazaar_field_encoding_adapter_failure;
  }
  if (!kvazaar_encode(enc, second_field, &second.data_out, &second.len_out, NULL, NULL, NULL)) {
    goto kvazaar_field_encoding_adapter_failure;
  }

  kvz_image_free(first_field);
  kvz_image_free(second_field);

  // Concatenate bitstreams.
  if (len_out != NULL) {
    *len_out = first.len_out + second.len_out;
  }
  if (data_out != NULL) {
    *data_out = first.data_out;
    if (first.data_out != NULL) {
      kvz_data_chunk *chunk = first.data_out;
      while (chunk->next != NULL) {
        chunk = chunk->next;
      }
      chunk->next = second.data_out;
    }
  }

  if (src_out != NULL) {
    // TODO: deinterlace the fields to one picture.
  }

  return 1;

kvazaar_field_encoding_adapter_failure:
  kvz_image_free(first_field);
  kvz_image_free(second_field);
  kvz_bitstream_free_chunks(first.data_out);
  kvz_bitstream_free_chunks(second.data_out);
  return 0;
}


static const kvz_api kvz_8bit_api = {
  .config_alloc = kvz_config_alloc,
  .config_init = kvz_config_init,
  .config_destroy = kvz_config_destroy,
  .config_parse = kvz_config_parse,

  .picture_alloc = kvz_image_alloc_420,
  .picture_free = kvz_image_free,

  .chunk_free = kvz_bitstream_free_chunks,

  .encoder_open = kvazaar_open,
  .encoder_close = kvazaar_close,
  .encoder_headers = kvazaar_headers,
  .encoder_encode = kvazaar_field_encoding_adapter,

  .picture_alloc_csp = kvz_image_alloc,
};


const kvz_api * kvz_api_get(int bit_depth)
{
  return &kvz_8bit_api;
}
