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

    kvz_encoder_control_free(encoder->control);
    encoder->control = NULL;

    // ***********************************************
    // Modified for SHVC
    if (encoder->el_states) {
      for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
        kvz_encoder_state_finalize(&encoder->el_states[i]);
      }
    }
    FREE_POINTER(encoder->el_states);

    const kvz_config* el_cfg = encoder->el_control->cfg;
    kvz_encoder_control_free(encoder->el_control);
    encoder->el_control = NULL;
    kvz_config_destroy((kvz_config*)el_cfg); //TODO: Figure out a better way
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
  kvz_config* el_cfg = kvz_config_alloc();
  *el_cfg = *cfg;
  el_cfg->qp = 1;
  el_cfg->width = el_cfg->el_width;
  el_cfg->height = el_cfg->el_height;
  encoder->el_control = kvz_encoder_control_init(el_cfg);

  encoder->cur_el_state_num = 0;
  encoder->out_el_state_num = 0;
  encoder->el_frames_started = 0;
  encoder->el_frames_done = 0;

  kvz_init_input_frame_buffer(&encoder->el_input_buffer);

  encoder->el_states = calloc(encoder->num_encoder_states, sizeof(encoder_state_t));
  if (!encoder->el_states) {
    goto kvazaar_open_failure;
  }

  for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
    encoder->el_states[i].encoder_control = encoder->el_control;

    if (!kvz_encoder_state_init(&encoder->el_states[i], NULL)) {
      goto kvazaar_open_failure;
    }

    encoder->el_states[i].global->QP = (int8_t)el_cfg->qp;
  }

  for (int i = 0; i < encoder->num_encoder_states; ++i) {
    
    encoder->el_states[i].previous_encoder_state = &encoder->el_states[abs((i - 1) % encoder->num_encoder_states)];
    
    kvz_encoder_state_match_children_of_previous_frame(&encoder->el_states[i]);
  }

  encoder->el_states[encoder->cur_el_state_num].global->frame = -1;
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
static void set_el_frame_info(kvz_frame_info *const info, const encoder_state_t *const state)
{
  info->el_qp = state->global->QP;
  info->el_nal_unit_type = state->global->pictype;
  info->el_slice_type = state->global->slicetype;
  kvz_encoder_get_ref_lists(state, info->el_ref_list_len, info->el_ref_list);
}

//TODO: Reuse buffers? Or not, who cares. Use a scaler struct to hold all relevant info for different layers?
//Create a new kvz picture based on pic_in with size given by width and height
kvz_picture* kvazaar_scaling(const kvz_picture* const pic_in, int width, int height)
{
  //Create the buffers that are passed to the scaling function
  //TODO: Consider the case when kvz_pixel is not uint8
  assert(pic_in->width==pic_in->stride); //Should be equal or the data transfer may fail.
  
  yuv_buffer_t* src_pic = newYuvBuffer_uint8(pic_in->y, pic_in->u, pic_in->v, pic_in->width, pic_in->height, CHROMA_420, 0); //TODO: Get chroma format from cfg?
  scaling_parameter_t param = newScalingParameters(pic_in->width, pic_in->height, width, height, CHROMA_420);
  yuv_buffer_t* trgt_pic = yuvScaling(src_pic, &param, NULL );
  
  //Create a new kvz picture from the buffer
  kvz_picture* pic_out = kvz_image_alloc(width, height);
  pic_out->dts = pic_in->dts;
  pic_out->pts = pic_in->pts;

  //Copy data to kvz picture
  int luma_size = width*height;
  int chroma_size = luma_size/4;
  int full_size = luma_size + chroma_size*2;
  for(int i = 0; i < full_size; i++) {
    pic_out->fulldata[i] = i < luma_size ? trgt_pic->y->data[i] : (i < luma_size+chroma_size ? trgt_pic->u->data[i-luma_size] : trgt_pic->v->data[i-luma_size-chroma_size]);
  }

  //Do deallocation
  deallocateYuvBuffer(src_pic);
  deallocateYuvBuffer(trgt_pic);

  return pic_out;
}

int kvazaar_scalable_encode(kvz_encoder* enc, kvz_picture* pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
{
  //DO scaling here
  kvz_picture* bl_pic_in = kvazaar_scaling(pic_in, enc->control->in.width, enc->control->in.height);//pic_in->width/2, pic_in->height/2); 

  //Encode Bl first
  if (!kvazaar_encode(enc, bl_pic_in, data_out, len_out, pic_out, src_out, info_out)) {
    return 0;
  }

  //TODO: checks ?

  //Check if 
  //if (data_out) *data_out = NULL;
  //if (len_out) *len_out = 0;
  //if (pic_out) *pic_out = NULL;
  //if (src_out) *src_out = NULL;

  //Calculate data for El

  encoder_state_t *state = &enc->el_states[enc->cur_el_state_num];
  

  if (!state->prepared) {
    kvz_encoder_next_frame(state);

    //Also add base layer to the reference list
    encoder_state_t *bl_state = &enc->states[enc->cur_el_state_num]; //Should return the bl state with the same poc as state.
    assert( state->global->poc == bl_state->global->poc );
    //TODO: Add upscaling, Handle memory leak of kvz_cu_array_?
    kvz_image_list_add( state->global->ref, kvazaar_scaling(bl_state->tile->frame->rec,pic_in->width,pic_in->height), kvz_cu_array_alloc(pic_in->width,pic_in->height), bl_state->global->poc);//bl_state->tile->frame->cu_array, bl_state->global->poc );
  }

  if (pic_in != NULL) {
    // FIXME: The frame number printed here is wrong when GOP is enabled.
    CHECKPOINT_MARK("read source frame: %d", state->global->frame + enc->el_control->cfg->seek);
  }

  if (kvz_encoder_feed_frame(&enc->el_input_buffer, state, pic_in)) {
    assert(state->global->frame == enc->el_frames_started);
    // Start encoding.
    kvz_encode_one_frame(state);
    enc->el_frames_started += 1;
  }

  // If we have finished encoding as many frames as we have started, we are done.
  if (enc->el_frames_done == enc->el_frames_started) {
    return 1;
  }

  if (!state->frame_done) {
    // We started encoding a frame; move to the next encoder state.
    enc->cur_el_state_num = (enc->cur_el_state_num + 1) % (enc->num_encoder_states);
  }

  encoder_state_t *output_state = &enc->el_states[enc->out_el_state_num];
  if (!output_state->frame_done &&
      (pic_in == NULL || enc->cur_el_state_num == enc->out_el_state_num)) {

    kvz_threadqueue_waitfor(enc->el_control->threadqueue, output_state->tqj_bitstream_written);
    // The job pointer must be set to NULL here since it won't be usable after
    // the next frame is done.
    output_state->tqj_bitstream_written = NULL;

    // Get stream length before taking chunks since that clears the stream.
    if (len_out) *len_out += kvz_bitstream_tell(&output_state->stream) / 8;
    if (data_out) {
      //Concatenate data to the end of data_out from bl
      kvz_data_chunk* bl_end_chunk = *data_out;
      while(bl_end_chunk->next != NULL) {
        bl_end_chunk = bl_end_chunk->next;
      }

      bl_end_chunk->next = kvz_bitstream_take_chunks(&output_state->stream);
    }
    if (pic_out) {
      //Set el pic_out as the base pic for bl out pic
      (*pic_out)->base_image = kvz_image_copy_ref(output_state->tile->frame->rec);
    }
    if (src_out) {
      //Set el src_out as the base pic for bl out pic
      (*src_out)->base_image = kvz_image_copy_ref(output_state->tile->frame->source);
    }
    if (info_out) {
      set_el_frame_info(info_out, output_state);
    }

    output_state->frame_done = 1;
    output_state->prepared = 0;
    enc->el_frames_done += 1;

    enc->out_el_state_num = (enc->out_el_state_num + 1) % (enc->num_encoder_states);
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
