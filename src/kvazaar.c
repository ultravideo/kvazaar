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

#include <stdlib.h>

#include "config.h"
#include "encoder.h"
#include "strategyselector.h"
#include "encoderstate.h"
#include "checkpoint.h"


static void kvazaar_close(kvz_encoder *encoder)
{
  if (encoder) {
    if (encoder->control) {
      encoder_control_finalize(encoder->control);
    }
    FREE_POINTER(encoder->control);
    FREE_POINTER(encoder->states);
  }
  FREE_POINTER(encoder);
}


static kvz_encoder * kvazaar_open(config_t *cfg)
{
  kvz_encoder *encoder = NULL;

  //Initialize strategies
  // TODO: Make strategies non-global
  if (!strategyselector_init(cfg->cpuid)) {
    fprintf(stderr, "Failed to initialize strategies.\n");
    goto kvazaar_open_failure;
  }

  //Allocate and init exp golomb table
  if (!init_exp_golomb(4096 * 8)) {
    fprintf(stderr, "Failed to allocate the exp golomb code table, shutting down!\n");
    goto kvazaar_open_failure;
  }

  encoder = MALLOC(kvz_encoder, 1);
  if (!encoder) {
    goto kvazaar_open_failure;
  }

  encoder->control = MALLOC(encoder_control_t, 1);
  if (!encoder->control || !encoder_control_init(encoder->control, cfg)) {
    goto kvazaar_open_failure;
  }

  encoder->num_encoder_states = cfg->owf + 1;
  encoder->cur_state_num = 0;
  encoder->frames_started = 0;
  encoder->states = MALLOC(encoder_state_t, encoder->num_encoder_states);
  if (!encoder->states) {
    goto kvazaar_open_failure;
  }

  for (unsigned i = 0; i <= cfg->owf; ++i) {
    encoder->states[i].encoder_control = encoder->control;

    if (!encoder_state_init(&encoder->states[i], NULL)) {
      goto kvazaar_open_failure;
    }

    encoder->states[i].global->QP = (int8_t)cfg->qp;
  }

  for (int i = 0; i <= cfg->owf; ++i) {
    if (i == 0) {
      encoder->states[i].previous_encoder_state = &encoder->states[encoder->num_encoder_states - 1];
    } else {
      encoder->states[i].previous_encoder_state = &encoder->states[(i - 1) % encoder->num_encoder_states];
    }
    encoder_state_match_children_of_previous_frame(&encoder->states[i]);
  }

  encoder->states[encoder->cur_state_num].global->frame = -1;

  return encoder;

kvazaar_open_failure:
  kvazaar_close(encoder);
  return NULL;
}


static int kvazaar_encode(kvz_encoder *enc, kvz_picture *img_in, kvz_picture **img_out, kvz_payload **payload)
{
  // If img_in is NULL, just return the next unfinished frame.
  if (img_in != NULL) {
    encoder_state_t *state = &enc->states[enc->cur_state_num];

    enc->frames_started += 1;
    encoder_next_frame(state, img_in);

    CHECKPOINT_MARK("read source frame: %d", state->global->frame + enc->control->cfg->seek);

    // The actual coding happens here, after this function we have a coded frame
    encode_one_frame(state);
  }

  enc->cur_state_num = (enc->cur_state_num + 1) % (enc->num_encoder_states);
  encoder_state_t *state = &enc->states[enc->cur_state_num];

  if (enc->frames_started >= enc->num_encoder_states && !state->stats_done) {
    threadqueue_waitfor(enc->control->threadqueue, state->tqj_bitstream_written);
    *img_out = image_make_subimage(state->tile->frame->rec, 0, 0, state->tile->frame->width, state->tile->frame->height);
  }

  return 0;
}

kvz_api kvz_8bit_api = {
  .config_alloc = config_alloc,
  .config_init = config_init,
  .config_read = config_read,
  .config_destroy = config_destroy,

  .picture_create = NULL,
  .picture_destroy = NULL,

  .encoder_open = kvazaar_open,
  .encoder_close = kvazaar_close,
  .encoder_encode = kvazaar_encode,
};


const kvz_api * kvz_api_get(int bit_depth)
{
  return &kvz_8bit_api;
}
