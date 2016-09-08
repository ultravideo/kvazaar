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
  if (encoder) {
    if (encoder->states) {
      for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
        kvz_encoder_state_finalize(&encoder->states[i]);
      }
    }
    FREE_POINTER(encoder->states);

    // ***********************************************
    // Modified for SHVC
    int layers = encoder->control->cfg->max_layers;
    // ***********************************************

    kvz_encoder_control_free(encoder->control);
    encoder->control = NULL;

    // ***********************************************
    // Modified for SHVC
    for (int layer_id_minus1 = 0; layer_id_minus1 < layers-1; layer_id_minus1++) {
      if (encoder->el_states[layer_id_minus1]) {
        for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
          kvz_encoder_state_finalize(&encoder->el_states[layer_id_minus1][i]);
        }
      }
      FREE_POINTER(encoder->el_states[layer_id_minus1]);

      //kvz_config* el_cfg = (kvz_config*)encoder->el_control[layer_id_minus1]->cfg;
      kvz_encoder_control_free(encoder->el_control[layer_id_minus1]);
      encoder->el_control = NULL;
      //kvz_config_destroy(el_cfg); //TODO: Figure out a better way 
    }
    FREE_POINTER(encoder->el_control);
    FREE_POINTER(encoder->el_states);
    FREE_POINTER(encoder->cur_el_state_num);
    FREE_POINTER(encoder->out_el_state_num);
    FREE_POINTER(encoder->el_input_buffer);
    FREE_POINTER(encoder->el_frames_started);
    FREE_POINTER(encoder->el_frames_done);
    FREE_POINTER(encoder->upscaling);
    FREE_POINTER(encoder->downscaling);
    
  // ***********************************************
  }
  FREE_POINTER(encoder);

  

}



static kvz_encoder * kvazaar_open(const kvz_config *cfg)
{
  kvz_encoder *encoder = NULL;

  //Initialize strategies
  // TODO: Make strategies non-global
  if (!kvz_strategyselector_init(cfg->cpuid, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "Failed to initialize strategies.\n");
    goto kvazaar_open_failure;
  }

  kvz_init_exp_golomb();

  encoder = calloc(1, sizeof(kvz_encoder));
  if (!encoder) {
    goto kvazaar_open_failure;
  }

  // FIXME: const qualifier disgarded. I don't want to change kvazaar_open
  // but I really need to change cfg.
  encoder->control = kvz_encoder_control_init((kvz_config*)cfg);
  if (!encoder->control) {
    goto kvazaar_open_failure;
  }

  encoder->num_encoder_states = encoder->control->owf + 1;
  encoder->cur_state_num = 0;
  encoder->out_state_num = 0;
  encoder->frames_started = 0;
  encoder->frames_done = 0;

  // ***********************************************
  // Modified for SHVC
  //TODO: Make a better implementaino. el_cfg is needed to pass the layer id, figure out a better way
  //TODO: Add error checking
  //Allocate the needed arrays etc.
  int el_layers = cfg->max_layers-1;

  if( el_layers > 0 ) {  
    encoder->el_control = MALLOC(encoder_control_t*, el_layers);
    encoder->el_states = MALLOC(encoder_state_t*, el_layers);
    encoder->cur_el_state_num = MALLOC(unsigned, el_layers);
    encoder->out_el_state_num = MALLOC(unsigned, el_layers);
    encoder->el_input_buffer = MALLOC(input_frame_buffer_t, el_layers);
    encoder->el_frames_started = MALLOC(unsigned, el_layers);
    encoder->el_frames_done = MALLOC(unsigned, el_layers);
    encoder->downscaling = MALLOC(scaling_parameter_t, el_layers+1);
    encoder->upscaling = MALLOC(scaling_parameter_t, el_layers+1);

    if (!encoder->el_control || !encoder->el_states || !encoder->cur_el_state_num ||
      !encoder->out_el_state_num || !encoder->el_input_buffer || !encoder->el_frames_started ||
      !encoder->el_frames_done || !encoder->downscaling || !encoder->upscaling ) {
      goto kvazaar_open_failure;
    }

    //Set scaling param for base layer
    //Need to use the padded size
    uint8_t padding_x = (CU_MIN_SIZE_PIXELS - cfg->in_width % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    uint8_t padding_y = (CU_MIN_SIZE_PIXELS - cfg->in_height % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    encoder->downscaling[0] = newScalingParameters(cfg->in_width + padding_x,cfg->in_height + padding_y,encoder->control->in.width,encoder->control->in.height,CHROMA_420); //TODO: get proper width/height for each layer from cfg etc.
    encoder->upscaling[0] = newScalingParameters(encoder->control->in.width,encoder->control->in.height,encoder->control->in.width,encoder->control->in.height,CHROMA_420);
  }
  else {
    encoder->el_control = NULL;
    encoder->el_states = NULL;
    encoder->cur_el_state_num = NULL;
    encoder->out_el_state_num = NULL;
    encoder->el_input_buffer = NULL;
    encoder->el_frames_started = NULL;
    encoder->el_frames_done = NULL;
    encoder->downscaling = NULL;
    encoder->upscaling = NULL;
  }

  
  for (int layer_id_minus1 = 0; layer_id_minus1 < el_layers; layer_id_minus1++) {
    //TODO: Set correct size etc. based on cfg
    //kvz_config* el_cfg = kvz_config_alloc();
    //kvz_config_init(el_cfg);
    //*el_cfg = *cfg;
    //el_cfg->qp = 5;
    //el_cfg->width = el_cfg->el_width;
    //el_cfg->height = el_cfg->el_height;
    //el_cfg->layer = layer_id_minus1;
    
    encoder->el_control[layer_id_minus1] = kvz_encoder_control_init(cfg->el_cfg[layer_id_minus1]);
    
    encoder->cur_el_state_num[layer_id_minus1] = 0;
    encoder->out_el_state_num[layer_id_minus1] = 0;
    encoder->el_frames_started[layer_id_minus1] = 0;
    encoder->el_frames_done[layer_id_minus1] = 0;

    kvz_init_input_frame_buffer(&encoder->el_input_buffer[layer_id_minus1]);

    encoder->el_states[layer_id_minus1] = calloc(encoder->num_encoder_states, sizeof(encoder_state_t));
    if (!encoder->el_states[layer_id_minus1]) {
      goto kvazaar_open_failure;
    }

    for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
      encoder->el_states[layer_id_minus1][i].encoder_control = encoder->el_control[layer_id_minus1];

      if (!kvz_encoder_state_init(&encoder->el_states[layer_id_minus1][i], NULL)) {
        goto kvazaar_open_failure;
      }

      encoder->el_states[layer_id_minus1][i].global->QP = (int8_t)cfg->el_cfg[layer_id_minus1]->qp;
    }

    for (int i = 0; i < encoder->num_encoder_states; ++i) {

      encoder->el_states[layer_id_minus1][i].previous_encoder_state = &encoder->el_states[layer_id_minus1][abs((i - 1) % encoder->num_encoder_states)];

      kvz_encoder_state_match_children_of_previous_frame(&encoder->el_states[layer_id_minus1][i]);
    }

    encoder->el_states[layer_id_minus1][encoder->cur_el_state_num[layer_id_minus1]].global->frame = -1;

    //Prepare scaling parameters so that up/downscaling[layer_id] gives the correct parameters for up/downscaling from orig/prev_layer to layer_id
    //Need to use the padded size
    uint8_t padding_x = (CU_MIN_SIZE_PIXELS - cfg->in_width % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    uint8_t padding_y = (CU_MIN_SIZE_PIXELS - cfg->in_height % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    encoder->downscaling[layer_id_minus1 + 1] = newScalingParameters(cfg->in_width + padding_x,
                                                                     cfg->in_height + padding_y,
                                                                     encoder->el_control[layer_id_minus1]->in.width,
                                                                     encoder->el_control[layer_id_minus1]->in.height,
                                                                     CHROMA_420); //TODO: get proper width/height for each layer from cfg etc.
    encoder->upscaling[layer_id_minus1 + 1] = newScalingParameters(encoder->upscaling[layer_id_minus1].trgt_width,
                                                                   encoder->upscaling[layer_id_minus1].trgt_height,
                                                                   encoder->el_control[layer_id_minus1]->in.width,
                                                                   encoder->el_control[layer_id_minus1]->in.height,
                                                                   CHROMA_420); //TODO: Account for irrecular reference structures?
  }
  // ***********************************************
  


  kvz_init_input_frame_buffer(&encoder->input_buffer);

  encoder->states = calloc(encoder->num_encoder_states, sizeof(encoder_state_t));
  if (!encoder->states) {
    goto kvazaar_open_failure;
  }

  for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
    encoder->states[i].encoder_control = encoder->control;

    if (!kvz_encoder_state_init(&encoder->states[i], NULL)) {
      goto kvazaar_open_failure;
    }

    encoder->states[i].frame->QP = (int8_t)cfg->qp;
  }

  for (int i = 0; i < encoder->num_encoder_states; ++i) {
    if (i == 0) {
      encoder->states[i].previous_encoder_state = &encoder->states[encoder->num_encoder_states - 1];
    } else {
      encoder->states[i].previous_encoder_state = &encoder->states[(i - 1) % encoder->num_encoder_states];
    }
    kvz_encoder_state_match_children_of_previous_frame(&encoder->states[i]);
  }

  encoder->states[encoder->cur_state_num].frame->num = -1;

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
  kvz_encoder_get_ref_lists(state, info->ref_list_len, info->ref_list);
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

  if (!state->prepared) {
    kvz_encoder_prepare(state);
  }

  if (pic_in != NULL) {
    // FIXME: The frame number printed here is wrong when GOP is enabled.
    CHECKPOINT_MARK("read source frame: %d", state->frame->num + enc->control->cfg->seek);
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

  if (!state->frame_done) {
    // We started encoding a frame; move to the next encoder state.
    enc->cur_state_num = (enc->cur_state_num + 1) % (enc->num_encoder_states);
  }

  encoder_state_t *output_state = &enc->states[enc->out_state_num];
  if (!output_state->frame_done &&
      (pic_in == NULL || enc->cur_state_num == enc->out_state_num)) {

    kvz_threadqueue_waitfor(enc->control->threadqueue, output_state->tqj_bitstream_written);
    // The job pointer must be set to NULL here since it won't be usable after
    // the next frame is done.
    output_state->tqj_bitstream_written = NULL;

    // Get stream length before taking chunks since that clears the stream.
    if (len_out) *len_out = kvz_bitstream_tell(&output_state->stream) / 8;
    if (data_out) *data_out = kvz_bitstream_take_chunks(&output_state->stream);
    if (pic_out) *pic_out = kvz_image_copy_ref(output_state->tile->frame->rec);
    if (src_out) *src_out = kvz_image_copy_ref(output_state->tile->frame->source);
    if (info_out) set_frame_info(info_out, output_state);

    output_state->frame_done = 1;
    output_state->prepared = 0;
    enc->frames_done += 1;

    enc->out_state_num = (enc->out_state_num + 1) % (enc->num_encoder_states);
  }

  return 1;
}
// ***********************************************
  // Modified for SHVC
//static void set_el_frame_info(kvz_frame_info *const info, const encoder_state_t *const state, int layer_id_minus1)
//{
//
//
//  info->el_qp[layer_id_minus1] = state->global->QP;
//  info->el_nal_unit_type[layer_id_minus1] = state->global->pictype;
//  info->el_slice_type[layer_id_minus1] = state->global->slicetype;
//  kvz_encoder_get_ref_lists(state, info->el_ref_list_len[layer_id_minus1], info->el_ref_list[layer_id_minus1]);
//}

//TODO: Reuse buffers? Or not, who cares. Use a scaler struct to hold all relevant info for different layers?
//Create a new kvz picture based on pic_in with size given by width and height
kvz_picture* kvazaar_scaling(const kvz_picture* const pic_in, scaling_parameter_t* param)
{
  //Create the buffers that are passed to the scaling function
  //TODO: Consider the case when kvz_pixel is not uint8
  assert(pic_in->width==pic_in->stride); //Should be equal or the data transfer may fail.
  assert(pic_in->width==param->src_width); //in pic size should match the param size
  assert(pic_in->height==param->src_height);
  
  yuv_buffer_t* src_pic = newYuvBuffer_uint8(pic_in->y, pic_in->u, pic_in->v, param->src_width, param->src_height, param->chroma, 0);
  yuv_buffer_t* trgt_pic = yuvScaling(src_pic, param, NULL );
  
  if( trgt_pic == NULL ) {
    deallocateYuvBuffer(src_pic);
    return NULL;
  }

  //Create a new kvz picture from the buffer
  kvz_picture* pic_out = kvz_image_alloc(param->trgt_width, param->trgt_height);
  if( pic_out == NULL) {
    deallocateYuvBuffer(src_pic);
    deallocateYuvBuffer(trgt_pic);
    return NULL;
  }
  pic_out->dts = pic_in->dts;
  pic_out->pts = pic_in->pts;

  //Copy data to kvz picture
  int luma_size = param->trgt_width*param->trgt_height;
  int chroma_size = luma_size/(param->chroma == CHROMA_420 ? 4 : param->chroma == CHROMA_444 ? 1 : 0);
  int full_size = luma_size + chroma_size*2;
  for(int i = 0; i < full_size; i++) {
    pic_out->fulldata[i] = i < luma_size ? trgt_pic->y->data[i] : (i < luma_size+chroma_size ? trgt_pic->u->data[i-luma_size] : trgt_pic->v->data[i-luma_size-chroma_size]);
  }

  //Do deallocation
  deallocateYuvBuffer(src_pic);
  deallocateYuvBuffer(trgt_pic);

  return pic_out;
}

//TODO: make a note of this: Asume that info_out is an array with an element for each layer
//TODO: Allow scaling "step-wise" instead of allways from the original, for a potentially reduced complexity
int kvazaar_scalable_encode(kvz_encoder* enc, kvz_picture* pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
{
  //DO scaling here
  //Pic_in for the layer being currently encoded
  kvz_picture* l_pic_in = pic_in == NULL ? NULL : kvazaar_scaling(pic_in, &enc->downscaling[0]);//pic_in->width/2, pic_in->height/2); 

  //Encode Bl first
  if (!kvazaar_encode(enc, l_pic_in, data_out, len_out, pic_out, src_out, info_out)) {
    return 0;
  }

  kvz_image_free(l_pic_in);

  //TODO: checks ?

  //Check if 
  //if (data_out) *data_out = NULL;
  //if (len_out) *len_out = 0;
  //if (pic_out) *pic_out = NULL;
  //if (src_out) *src_out = NULL;

  //Store the pic and data pointers to the most resently encoded el layer to allow chaining them
  kvz_data_chunk* last_l_chunk = *data_out;
  kvz_picture* last_l_pic_out = *pic_out;
  kvz_picture* last_l_src_out = *src_out;

  //Calculate data for Els
  for (int layer_id_minus1 = 0; layer_id_minus1 < enc->control->cfg->max_layers-1; layer_id_minus1++) {

    unsigned* cur_el_state_num = &enc->cur_el_state_num[layer_id_minus1];
    unsigned* out_el_state_num = &enc->cur_el_state_num[layer_id_minus1];
    unsigned* el_frames_started = &enc->el_frames_started[layer_id_minus1];
    unsigned* el_frames_done = &enc->el_frames_done[layer_id_minus1];
    encoder_state_t *state = &enc->el_states[layer_id_minus1][*cur_el_state_num];

    kvz_picture* last_l_pic_in = pic_in; //This is the src frame that should be used when downscaling. TODO: Use higher layers pic to speed up downscaling?
    l_pic_in = last_l_pic_in == NULL ? NULL : kvazaar_scaling(last_l_pic_in, &enc->downscaling[layer_id_minus1+1]);

    if (!state->prepared) {

      kvz_encoder_next_frame(state);

      //Also add base layer to the reference list. Need to do it here so that the ILR is after negative delta pocs
      encoder_state_t *bl_state = &enc->states[*cur_el_state_num]; //Should return the bl state with the same poc as state.
      assert(state->global->poc == bl_state->global->poc);
      //TODO: Add upscaling, Handle memory leak of kvz_cu_array_?
      //Skip on first frame?
      if (state->global->frame > 0) {
        kvz_image_list_add_back(state->global->ref,
                           kvazaar_scaling(bl_state->tile->frame->rec, &enc->upscaling[layer_id_minus1 + 1]),
                           bl_state->tile->frame->cu_array, //kvz_cu_array_alloc(enc->upscaling[layer_id_minus1 + 1].trgt_width, enc->upscaling[layer_id_minus1 + 1].trgt_height),
                           bl_state->global->poc);//bl_state->tile->frame->cu_array, bl_state->global->poc );//
      }
    }

    if (l_pic_in != NULL) {
      // FIXME: The frame number printed here is wrong when GOP is enabled.
      CHECKPOINT_MARK("read source frame: %d", state->global->frame + enc->el_control[layer_id_minus1]->cfg->seek);
    }

    if (kvz_encoder_feed_frame(&enc->el_input_buffer[layer_id_minus1], state, l_pic_in)) {
      assert(state->global->frame == *el_frames_started);
      // Start encoding.
      kvz_encode_one_frame(state);
      *el_frames_started += 1;
    }

    // If we have finished encoding as many frames as we have started, we are done.
    if (*el_frames_done == *el_frames_started) {
      return 1;
    }

    if (!state->frame_done) {
      // We started encoding a frame; move to the next encoder state.
      *cur_el_state_num = (*cur_el_state_num + 1) % (enc->num_encoder_states);
    }

    encoder_state_t *output_state = &enc->el_states[layer_id_minus1][*out_el_state_num];
    if (!output_state->frame_done &&
      (l_pic_in == NULL || *cur_el_state_num == *out_el_state_num)) {

      kvz_threadqueue_waitfor(enc->el_control[layer_id_minus1]->threadqueue, output_state->tqj_bitstream_written);
      // The job pointer must be set to NULL here since it won't be usable after
      // the next frame is done.
      output_state->tqj_bitstream_written = NULL;

      // Get stream length before taking chunks since that clears the stream.
      if (len_out) *len_out += kvz_bitstream_tell(&output_state->stream) / 8;
      if (last_l_chunk) {
        //Concatenate data to the end of last layer's data chunk
        while (last_l_chunk->next != NULL) {
          last_l_chunk = last_l_chunk->next;
        }

        last_l_chunk->next = kvz_bitstream_take_chunks(&output_state->stream);
      }
      if (last_l_pic_out) {
        //Set cur layer pic_out as the base pic for last_l out pic
        last_l_pic_out->base_image = kvz_image_copy_ref(output_state->tile->frame->rec);
        last_l_pic_out = last_l_pic_out->base_image;
      }
      if (last_l_src_out) {
        //Set cur layer src_out as the base pic for last_l out pic
        last_l_src_out->base_image = kvz_image_copy_ref(output_state->tile->frame->source);
        last_l_src_out = last_l_src_out->base_image;
      }
      if (&info_out[layer_id_minus1+1]) {
        set_frame_info(&info_out[layer_id_minus1+1], output_state);
      }

      output_state->frame_done = 1;
      output_state->prepared = 0;
      *el_frames_done += 1;

      *out_el_state_num = (*out_el_state_num + 1) % (enc->num_encoder_states);
    }

    kvz_image_free(l_pic_in);
  }
  return 1;
}


//TODO: Remove dublicated code? Handle interlaced
static int kvazaar_field_encoding_adapter(kvz_encoder *enc,
                                          kvz_picture *pic_in,
                                          kvz_data_chunk **data_out,
                                          uint32_t *len_out,
                                          kvz_picture **pic_out,
                                          kvz_picture **src_out,
                                          kvz_frame_info *info_out)
{
  if (enc->control->cfg->source_scan_type == KVZ_INTERLACING_NONE) {
    // For progressive, simply call the normal encoding function.
    //If several layers are used, call the apropriate function
    if(enc->control->cfg->max_layers > 1) return kvazaar_scalable_encode(enc, pic_in, data_out, len_out, pic_out, src_out, info_out);
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
