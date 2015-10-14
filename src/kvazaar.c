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

#include "kvazaar_internal.h"

#include <stdlib.h>

#include "config.h"
#include "encoder.h"
#include "strategyselector.h"
#include "encoderstate.h"
#include "checkpoint.h"
#include "bitstream.h"
#include "input_frame_buffer.h"


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

  encoder->control = kvz_encoder_control_init(cfg);
  if (!encoder->control) {
    goto kvazaar_open_failure;
  }

  encoder->num_encoder_states = encoder->control->owf + 1;
  encoder->cur_state_num = 0;
  encoder->out_state_num = 0;
  encoder->frames_started = 0;
  encoder->frames_done = 0;

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

    encoder->states[i].global->QP = (int8_t)cfg->qp;
  }

  for (int i = 0; i < encoder->num_encoder_states; ++i) {
    if (i == 0) {
      encoder->states[i].previous_encoder_state = &encoder->states[encoder->num_encoder_states - 1];
    } else {
      encoder->states[i].previous_encoder_state = &encoder->states[(i - 1) % encoder->num_encoder_states];
    }
    kvz_encoder_state_match_children_of_previous_frame(&encoder->states[i]);
  }

  encoder->states[encoder->cur_state_num].global->frame = -1;

  return encoder;

kvazaar_open_failure:
  kvazaar_close(encoder);
  return NULL;
}


static void set_frame_info(kvz_frame_info *const info, const encoder_state_t *const state)
{
  info->poc = state->global->poc,
  info->qp = state->global->QP;
  info->nal_unit_type = state->global->pictype;
  info->slice_type = state->global->slicetype;
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
    kvz_encoder_next_frame(state);
  }

  if (pic_in != NULL) {
    // FIXME: The frame number printed here is wrong when GOP is enabled.
    CHECKPOINT_MARK("read source frame: %d", state->global->frame + enc->control->cfg->seek);
  }

  if (kvz_encoder_feed_frame(&enc->input_buffer, state, pic_in)) {
    assert(state->global->frame == enc->frames_started);
    // Start encoding.
    kvz_encode_one_frame(state);
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


static const kvz_api kvz_8bit_api = {
  .config_alloc = kvz_config_alloc,
  .config_init = kvz_config_init,
  .config_destroy = kvz_config_destroy,
  .config_parse = kvz_config_parse,

  .picture_alloc = kvz_image_alloc,
  .picture_free = kvz_image_free,

  .chunk_free = kvz_bitstream_free_chunks,

  .encoder_open = kvazaar_open,
  .encoder_close = kvazaar_close,
  .encoder_headers = kvazaar_headers,
  .encoder_encode = kvazaar_encode,
};


const kvz_api * kvz_api_get(int bit_depth)
{
  return &kvz_8bit_api;
}
