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
#include <crtdbg.h>

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

    if (encoder->states) {
      for (unsigned i = 0; i < encoder->num_encoder_states; ++i) {
        kvz_encoder_state_finalize(&encoder->states[i]);
      }
    }
    FREE_POINTER(encoder->states);

    kvz_encoder_control_free(encoder->control);
    encoder->control = NULL;
  }
  FREE_POINTER(encoder);
 
  kvazaar_close(next);
  // ***********************************************
}



static kvz_encoder * kvazaar_open(const kvz_config *cfg)
{
  kvz_encoder *encoder = NULL; //The base layer encoder

  //Initialize strategies
  // TODO: Make strategies non-global
  if (!kvz_strategyselector_init(cfg->cpuid, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "Failed to initialize strategies.\n");
    goto kvazaar_open_failure;
  }

  kvz_init_exp_golomb();

  // ***********************************************
  // Modified for SHVC. TODO: Account for more complex ref structures?
  kvz_encoder *cur_enc = NULL;
  kvz_encoder *prev_enc = NULL;
  //TODO: Just use a while loop?
  for (unsigned j = 0; j < *cfg->max_layers; j++) {
    cur_enc = calloc(1, sizeof(kvz_encoder));
    if (!cur_enc) {
      goto kvazaar_open_failure;
    }
    if( j == 0 ) encoder = cur_enc;

    // FIXME: const qualifier disgarded. I don't want to change kvazaar_open
    // but I really need to change cfg.
    cur_enc->control = kvz_encoder_control_init((kvz_config*)cfg);
    if (!cur_enc->control) {
      goto kvazaar_open_failure;
    }

    cur_enc->num_encoder_states = cur_enc->control->owf + 1;
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
      cur_enc->states[i].encoder_control = encoder->control;

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

    cur_enc->states[cur_enc->cur_state_num].frame->num = -1;

    //Set scaling parameters
    //Prepare scaling parameters so that up/downscaling gives the correct parameters for up/downscaling from prev_layer/orig to current layer
    enum kvz_chroma_format csp = KVZ_FORMAT2CSP(cfg->input_format);
    cur_enc->downscaling = newScalingParameters(cfg->in_width,
                                                cfg->in_height,
                                                cur_enc->control->in.real_width,
                                                cur_enc->control->in.real_height,
                                                csp);
    if( prev_enc ){
      cur_enc->upscaling = newScalingParameters(prev_enc->upscaling.trgt_width,
                                                prev_enc->upscaling.trgt_height,
                                                cur_enc->control->in.real_width,
                                                cur_enc->control->in.real_height,
                                                csp);
    }
    else {
      cur_enc->upscaling = newScalingParameters(cur_enc->control->in.real_width,
                                                cur_enc->control->in.real_height,
                                                cur_enc->control->in.real_width,
                                                cur_enc->control->in.real_height,
                                                csp);
    }
    //Need to set the source (target?) to the padded size (because reasons) to conform with SHM. TODO: Trgt needs to be padded as well?
    //Scaling parameters need to be calculated for the true sizes.
    cur_enc->upscaling.src_padding_x = (CU_MIN_SIZE_PIXELS - cur_enc->upscaling.src_width % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    cur_enc->upscaling.src_padding_y = (CU_MIN_SIZE_PIXELS - cur_enc->upscaling.src_height % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;

    //Connect the sub encoders
    cur_enc->next_enc = NULL;
    cur_enc->prev_enc = prev_enc;
    if( prev_enc ) prev_enc->next_enc = cur_enc;
    
    //Prepare for the next loop
    prev_enc = cur_enc;
    cfg  = cfg->next_cfg;
  }

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
//  info->el_qp[layer_id_minus1] = state->frame->QP;
//  info->el_nal_unit_type[layer_id_minus1] = state->frame->pictype;
//  info->el_slice_type[layer_id_minus1] = state->frame->slicetype;
//  kvz_encoder_get_ref_lists(state, info->el_ref_list_len[layer_id_minus1], info->el_ref_list[layer_id_minus1]);
//}

//TODO: Move somewhere more appropriate.
void remove_ILR_pics( encoder_state_t* const state)
{
  //Loop over refs and remove IL refs from the list
  for( unsigned i = 0; i < state->frame->ref->used_size; i++) {
    //TODO: Figure out a better way? Eg. extra info.
    //If a ref_poc matches the prev frames poc, the ref should have been an ILR in the prev frame.
    if( state->frame->ref->pocs[i] == state->previous_encoder_state->frame->poc ) {
      kvz_image_list_rem( state->frame->ref, i);
    }
  }
}

//TODO: Reuse buffers? Or not, who cares. Use a scaler struct to hold all relevant info for different layers?
//TODO: remove memory db stuff
//Create a new kvz picture based on pic_in with size given by width and height
kvz_picture* kvazaar_scaling(const kvz_picture* const pic_in, scaling_parameter_t* param)
{
  //Create the buffers that are passed to the scaling function
  //TODO: Consider the case when kvz_pixel is not uint8
  //assert(pic_in->width==pic_in->stride); //Should be equal or the data transfer may fail.
  //assert(pic_in->width==param->src_width); //in pic size should match the param size
  //assert(pic_in->height==param->src_height);
  
  _ASSERTE( _CrtCheckMemory() );
  if( pic_in == NULL) {
    return NULL;
  }


  yuv_buffer_t* src_pic = newYuvBuffer_padded_uint8(pic_in->y, pic_in->u, pic_in->v, param->src_width+param->src_padding_x, param->src_height+param->src_padding_y, pic_in->stride, param->chroma, 0);
  //yuv_buffer_t* src_pic = newYuvBuffer_uint8(pic_in->y, pic_in->u, pic_in->v, pic_in->width, pic_in->height, param->chroma, 0);
  
  
  yuv_buffer_t* trgt_pic = yuvScaling(src_pic, param, NULL );
  
  _ASSERTE( _CrtCheckMemory() );

  if( trgt_pic == NULL ) {
    deallocateYuvBuffer(src_pic);
    return NULL;
  }
  
  //TODO: Add proper padding
  uint8_t padding_x = (CU_MIN_SIZE_PIXELS - param->trgt_width % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
  uint8_t padding_y = (CU_MIN_SIZE_PIXELS - param->trgt_height % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;

  //Create a new kvz picture from the buffer
  kvz_picture* pic_out = kvz_image_alloc(pic_in->chroma_format,param->trgt_width+padding_x, param->trgt_height+padding_y);
  if( pic_out == NULL) {
    deallocateYuvBuffer(src_pic);
    deallocateYuvBuffer(trgt_pic);
    return NULL;
  }
  pic_out->dts = pic_in->dts;
  pic_out->pts = pic_in->pts;

  //Copy data to kvz picture
  /*int luma_size = param->trgt_width*param->trgt_height;
  int chroma_size = luma_size/(param->chroma == CHROMA_420 ? 4 : param->chroma == CHROMA_444 ? 1 : 0);
  int full_size = luma_size + chroma_size*2;
  for(int i = 0; i < full_size; i++) {
    pic_out->fulldata[i] = i < luma_size ? trgt_pic->y->data[i] : (i < luma_size+chroma_size ? trgt_pic->u->data[i-luma_size] : trgt_pic->v->data[i-luma_size-chroma_size]);
  }*/

  int chroma_shift = param->chroma == CHROMA_444 ? 0 : 1;
  pic_data_t* comp_list[] = {trgt_pic->y->data, trgt_pic->u->data, trgt_pic->v->data};
  int stride_list[] = {trgt_pic->y->width,trgt_pic->u->width,trgt_pic->v->width};
  int height_list[] = {trgt_pic->y->height,trgt_pic->u->height,trgt_pic->v->height};
  int padd_x[] = {padding_x,padding_x>>chroma_shift,padding_x>>chroma_shift};
  int padd_y[] = {padding_y,padding_y>>chroma_shift,padding_y>>chroma_shift};
  assert(sizeof(kvz_pixel)==sizeof(char)); //Image copy (memset) only works if the pixels are the same size as char 
  
  _ASSERTE( _CrtCheckMemory() );

  //Loop over components
  for (int comp = 0, i = 0; comp < 3; comp++) {
    int comp_size = height_list[comp]*stride_list[comp];
    int pic_out_stride = pic_out->stride >> ( comp < 1 ? 0 : chroma_shift );
    for (int src_ind = 0; src_ind < comp_size; i++, src_ind++) {
      //TODO: go over src image correctly
      //TODO: Make a better loop
      //Copy value normally
      pic_out->fulldata[i] = comp_list[comp][src_ind];
        
      //_ASSERTE( _CrtCheckMemory() );
      if ( padding_x != 0 && (src_ind % stride_list[comp] == stride_list[comp]-1) ) { //Padd end of row by copying last pixel
        memset(pic_out->fulldata + i + 1, pic_out->fulldata[i], padd_x[comp]);
        //_ASSERTE( _CrtCheckMemory() );
        i += padd_x[comp];
      }
    }
    if (padd_y[comp] != 0 ) { //Padd image with lines copied from the prev row
        for (int j = 0; j < padd_y[comp]; j++) {
          memcpy(pic_out->fulldata + i, pic_out->fulldata + i - pic_out_stride, pic_out_stride);
          //_ASSERTE( _CrtCheckMemory() );
          i += pic_out_stride;
        }
    }
  }

  _ASSERTE( _CrtCheckMemory() );

  //Do deallocation
  deallocateYuvBuffer(src_pic);
  _ASSERTE( _CrtCheckMemory() );
  deallocateYuvBuffer(trgt_pic);
  _ASSERTE( _CrtCheckMemory() );
  return pic_out;
}

//TODO: make a note of this: Asume that info_out is an array with an element for each layer
//TODO: Allow scaling "step-wise" instead of allways from the original, for a potentially reduced complexity?
//TODO: Account for pic_in containing several input images for different layers
//Use this function to aggregate the results etc. but otherwise just call kvazaar_encode with the correct encoder
int kvazaar_scalable_encode(kvz_encoder* enc, kvz_picture* pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
{
  if (data_out) *data_out = NULL;
  if (len_out) *len_out = 0;
  if (pic_out) *pic_out = NULL;
  if (src_out) *src_out = NULL;

  //Use these to pass stuff to the actual encoder function and aggregate the results into the actual parameters
  kvz_encoder *cur_enc = enc;
  kvz_picture *cur_pic_in; 
  kvz_data_chunk* cur_data_out = NULL;
  uint32_t cur_len_out = 0;
  kvz_picture *cur_pic_out = NULL;
  kvz_picture *cur_src_out = NULL;

  //TODO: Use a while loop instead?
  for( unsigned i = 0; i < *enc->control->cfg->max_layers; i++) {
    
    cur_pic_in = kvazaar_scaling(pic_in, &cur_enc->downscaling);

    if(!kvazaar_encode(cur_enc, cur_pic_in, &cur_data_out, &cur_len_out, &cur_pic_out, &cur_src_out, &(info_out[i]))) {
      kvz_image_free(cur_pic_in);
      return 0;
    }

    kvz_image_free(cur_pic_in);
    cur_pic_in = NULL;

    //Aggregate new stuff
    if(data_out) {
      while((*data_out)->next != NULL) {
        data_out = &(*data_out)->next;
      }
      (*data_out)->next = cur_data_out;
      cur_data_out = NULL;
    }
    if(len_out) *len_out += cur_len_out;
    if(pic_out) {
      (*pic_out)->base_image = kvz_image_copy_ref(cur_pic_out);
      pic_out = &(*pic_out)->base_image;
      kvz_image_free(cur_pic_out);
      cur_pic_out = NULL;
    }
    if(src_out) {
      (*src_out)->base_image = kvz_image_copy_ref(cur_src_out);
      src_out = &(*src_out)->base_image;
      kvz_image_free(cur_src_out);
      cur_src_out = NULL;
    }

    //Update other values for the next layer
    cur_enc = cur_enc->next_enc;
  }
  
  return 1;
}

//TODO: Remove
//TODO: make a note of this: Asume that info_out is an array with an element for each layer
//TODO: Allow scaling "step-wise" instead of allways from the original, for a potentially reduced complexity?
//TODO: Merge with kvazaar_encode?
//int _kvazaar_scalable_encode(kvz_encoder* enc, kvz_picture* pic_in, kvz_data_chunk** data_out, uint32_t* len_out, kvz_picture** pic_out, kvz_picture** src_out, kvz_frame_info* info_out)
//{
//  //DO scaling here
//  //Pic_in for the layer being currently encoded
//  kvz_picture* l_pic_in = pic_in == NULL ? NULL : kvazaar_scaling(pic_in, &enc->downscaling[0]);//pic_in->width/2, pic_in->height/2); 
//
//  //Encode Bl first
//  if (!kvazaar_encode(enc, l_pic_in, data_out, len_out, pic_out, src_out, info_out)) {
//    return 0;
//  }
//
//  kvz_image_free(l_pic_in);
//
//  //TODO: checks ?
//
//  //Check if 
//  //if (data_out) *data_out = NULL;
//  //if (len_out) *len_out = 0;
//  //if (pic_out) *pic_out = NULL;
//  //if (src_out) *src_out = NULL;
//
//  //Store the pic and data pointers to the most resently encoded el layer to allow chaining them
//  kvz_data_chunk* last_l_chunk = *data_out;
//  kvz_picture* last_l_pic_out = *pic_out;
//  kvz_picture* last_l_src_out = *src_out;
//
//  //Calculate data for Els
//  for (int layer_id_minus1 = 0; layer_id_minus1 < *enc->control->cfg->max_layers-1; layer_id_minus1++) {
//
//    unsigned* cur_el_state_num = &enc->cur_el_state_num[layer_id_minus1];
//    unsigned* out_el_state_num = &enc->cur_el_state_num[layer_id_minus1];
//    unsigned* el_frames_started = &enc->el_frames_started[layer_id_minus1];
//    unsigned* el_frames_done = &enc->el_frames_done[layer_id_minus1];
//    encoder_state_t *state = &enc->el_states[layer_id_minus1][*cur_el_state_num];
//
//    kvz_picture* last_l_pic_in = pic_in; //This is the src frame that should be used when downscaling. TODO: Use higher layers pic to speed up downscaling?
//    l_pic_in = last_l_pic_in == NULL ? NULL : kvazaar_scaling(last_l_pic_in, &enc->downscaling[layer_id_minus1+1]);
//
//    if (!state->prepared) {
//      
//      //TODO: Find a better way.
//      //deallocate dummy cu_array
//      cu_array_t* cua = state->frame->ref->cu_arrays[0]; //ILR pic should be first
//      int used_size = state->frame->ref->used_size;
//      remove_ILR_pics(state); //Remove old ILR pics from the ref list so they don't interfere. TODO: Move somewhere else?
//      if( used_size > 0 ) kvz_cu_array_free(cua);
//      
//        kvz_encoder_prepare(state);
//
//      //TODO: Move somewhere else. slicetype still refers to the prev slice?
//      //TODO: Allow first EL layer to be a P-slice
//      if (state->frame->num > 0) {//(state->frame->slicetype != KVZ_SLICE_I) {
//        //Also add base layer to the reference list.
//        encoder_state_t *bl_state = &enc->states[*cur_el_state_num]; //Should return the bl state with the same poc as state.
//        //assert(state->frame->poc == bl_state->frame->poc);
//        //TODO: Add upscaling, Handle memory leak of kvz_cu_array_?
//        //Skip on first frame? Skip if inter frame. 
//        if (bl_state->tile->frame->rec != NULL) {
//          kvz_image_list_add/*_back*/(state->frame->ref,
//            kvazaar_scaling(bl_state->tile->frame->rec, &enc->upscaling[layer_id_minus1 + 1]),
//            /*bl_state->tile->frame->cu_array,*/ kvz_cu_array_alloc(enc->upscaling[layer_id_minus1 + 1].trgt_width, enc->upscaling[layer_id_minus1 + 1].trgt_height),
//            bl_state->frame->poc);//bl_state->tile->frame->cu_array, bl_state->frame->poc );//
//        }
//      }
//    }
//
//    if (l_pic_in != NULL) {
//      // FIXME: The frame number printed here is wrong when GOP is enabled.
//      CHECKPOINT_MARK("read source frame: %d", state->frame->frame + enc->el_control[layer_id_minus1]->cfg->seek);
//    }
//
//    kvz_picture* frame = kvz_encoder_feed_frame(&enc->el_input_buffer[layer_id_minus1], state, l_pic_in);
//    if (frame) {
//      assert(state->frame->num == *el_frames_started);
//      // Start encoding.
//      kvz_encode_one_frame(state, frame);
//
//      *el_frames_started += 1;
//    }
//
//    // If we have finished encoding as many frames as we have started, we are done.
//    if (*el_frames_done == *el_frames_started) {
//      return 1;
//    }
//
//    if (!state->frame_done) {
//      // We started encoding a frame; move to the next encoder state.
//      *cur_el_state_num = (*cur_el_state_num + 1) % (enc->num_encoder_states);
//    }
//
//    encoder_state_t *output_state = &enc->el_states[layer_id_minus1][*out_el_state_num];
//    if (!output_state->frame_done &&
//      (l_pic_in == NULL || *cur_el_state_num == *out_el_state_num)) {
//
//      kvz_threadqueue_waitfor(enc->el_control[layer_id_minus1]->threadqueue, output_state->tqj_bitstream_written);
//      // The job pointer must be set to NULL here since it won't be usable after
//      // the next frame is done.
//      output_state->tqj_bitstream_written = NULL;
//
//      // Get stream length before taking chunks since that clears the stream.
//      if (len_out) *len_out += kvz_bitstream_tell(&output_state->stream) / 8;
//      if (last_l_chunk) {
//        //Concatenate data to the end of last layer's data chunk
//        while (last_l_chunk->next != NULL) {
//          last_l_chunk = last_l_chunk->next;
//        }
//
//        last_l_chunk->next = kvz_bitstream_take_chunks(&output_state->stream);
//      }
//      if (last_l_pic_out) {
//        //Set cur layer pic_out as the base pic for last_l out pic
//        last_l_pic_out->base_image = kvz_image_copy_ref(output_state->tile->frame->rec);
//        last_l_pic_out = last_l_pic_out->base_image;
//      }
//      if (last_l_src_out) {
//        //Set cur layer src_out as the base pic for last_l out pic
//        last_l_src_out->base_image = kvz_image_copy_ref(output_state->tile->frame->source);
//        last_l_src_out = last_l_src_out->base_image;
//      }
//      if (&info_out[layer_id_minus1+1]) {
//        set_frame_info(&info_out[layer_id_minus1+1], output_state);
//      }
//
//      output_state->frame_done = 1;
//      output_state->prepared = 0;
//      *el_frames_done += 1;
//
//      *out_el_state_num = (*out_el_state_num + 1) % (enc->num_encoder_states);
//    }
//
//    kvz_image_free(l_pic_in);
//  }
//  return 1;
//}


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
    if(*enc->control->cfg->max_layers > 1) return kvazaar_scalable_encode(enc, pic_in, data_out, len_out, pic_out, src_out, info_out);
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
