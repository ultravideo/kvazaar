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

#include "input_frame_buffer.h"
#include "encoderstate.h"
#include <assert.h>

void kvz_init_input_frame_buffer(input_frame_buffer_t *input_buffer)
{
  FILL(input_buffer->pic_buffer, 0);
  input_buffer->pictures_available = 0;
  input_buffer->write_idx = 0;
  input_buffer->read_idx = 0;
  input_buffer->gop_offset = 0;
}

/**
 * \brief Pass an input frame to the encoder state.
 *
 * Sets the source image of the encoder state if there is a suitable image
 * available.
 *
 * The caller must not modify img_in after calling this function.
 *
 * \param buf     an input frame buffer
 * \param state   a main encoder state
 * \param img_in  input frame or NULL
 * \return        1 if the source image was set, 0 if not
 */
int kvz_encoder_feed_frame(input_frame_buffer_t *buf,
                           encoder_state_t *const state,
                           kvz_picture *const img_in)
{
  const encoder_control_t* const encoder = state->encoder_control;
  const kvz_config* const cfg = encoder->cfg;

  const int gop_buf_size = 2 * cfg->gop_len;

  assert(state->global->frame >= 0);
  assert(state->tile->frame->source == NULL);

  if (cfg->gop_len == 0 || state->global->frame == 0) {
    if (img_in == NULL) return 0;
    state->tile->frame->source = kvz_image_copy_ref(img_in);
    state->tile->frame->rec->pts = img_in->pts;
    state->global->gop_offset = 0;
    return 1;
  }
  // GOP enabled and not the first frame

  if (img_in != NULL) {
    // Save the input image in the buffer.
    assert(buf->pictures_available < gop_buf_size);
    assert(buf->pic_buffer[buf->write_idx] == NULL);
    buf->pic_buffer[buf->write_idx] = kvz_image_copy_ref(img_in);

    buf->pictures_available++;
    buf->write_idx++;
    if (buf->write_idx >= gop_buf_size) {
      buf->write_idx = 0;
    }
  }

  if (buf->pictures_available < cfg->gop_len) {
    if (img_in != NULL || buf->pictures_available == 0) {
      // Either start of the sequence with no full GOP available yet, or the
      // end of the sequence with all pics encoded.
      return 0;
    }
    // End of the sequence and a full GOP is not available.
    // Skip pictures until an available one is found.
    for (; buf->gop_offset < cfg->gop_len &&
           cfg->gop[buf->gop_offset].poc_offset - 1 >= buf->pictures_available;
           buf->gop_offset++) {}

    if (buf->gop_offset >= cfg->gop_len) {
      // All available pictures used.
      buf->gop_offset = 0;
      buf->pictures_available = 0;
      return 0;
    }
  }

  // Move image from buffer to state.
  int buffer_index = buf->read_idx + cfg->gop[buf->gop_offset].poc_offset - 1;
  assert(buf->pic_buffer[buffer_index] != NULL);
  assert(state->tile->frame->source == NULL);
  state->tile->frame->source = buf->pic_buffer[buffer_index];
  state->tile->frame->rec->pts = buf->pic_buffer[buffer_index]->pts;
  buf->pic_buffer[buffer_index] = NULL;

  state->global->gop_offset = buf->gop_offset;

  buf->gop_offset++;
  if (buf->gop_offset >= cfg->gop_len) {
    buf->gop_offset = 0;
    buf->pictures_available = MAX(0, buf->pictures_available - cfg->gop_len);
    buf->read_idx = (buf->read_idx + cfg->gop_len) % gop_buf_size;
  }

  return 1;
}
