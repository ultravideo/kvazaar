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

/**
 * \brief Pass an input frame to the encoder state.
 *
 * Sets the source image of the encoder state if there is a suitable image
 * available.
 *
 * The caller must not modify img_in after calling this function.
 *
 * \param state   a main encoder state
 * \param img_in  input frame or NULL
 * \return        1 if the source image was set, 0 if not
 */
int kvz_encoder_feed_frame(encoder_state_t *const state, kvz_picture *const img_in)
{
  const encoder_control_t* const encoder = state->encoder_control;
  const kvz_config* const cfg = encoder->cfg;

  // TODO: Get rid of static variables.
  static kvz_picture *gop_buffer[2 * KVZ_MAX_GOP_LENGTH] = { NULL };
  static int gop_buf_write_idx = 0;
  static int gop_buf_read_idx = 0;
  static int gop_pictures_available = 0;
  static int gop_offset = 0;

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
    assert(gop_pictures_available < gop_buf_size);
    assert(gop_buffer[gop_buf_write_idx] == NULL);
    gop_buffer[gop_buf_write_idx] = kvz_image_copy_ref(img_in);

    ++gop_pictures_available;
    if (++gop_buf_write_idx >= gop_buf_size) {
      gop_buf_write_idx = 0;
    }
  }

  if (gop_pictures_available < cfg->gop_len) {
    if (img_in != NULL || gop_pictures_available == 0) {
      // Either start of the sequence with no full GOP available yet, or the
      // end of the sequence with all pics encoded.
      return 0;
    }
    // End of the sequence and a full GOP is not available.
    // Skip pictures until an available one is found.
    for (; gop_offset < cfg->gop_len &&
           cfg->gop[gop_offset].poc_offset - 1 >= gop_pictures_available;
           ++gop_offset);

    if (gop_offset >= cfg->gop_len) {
      // All available pictures used.
      gop_offset = 0;
      gop_pictures_available = 0;
      return 0;
    }
  }

  // Move image from buffer to state.
  int buffer_index = gop_buf_read_idx + cfg->gop[gop_offset].poc_offset - 1;
  assert(gop_buffer[buffer_index] != NULL);
  assert(state->tile->frame->source == NULL);
  state->tile->frame->source = gop_buffer[buffer_index];
  state->tile->frame->rec->pts = gop_buffer[buffer_index]->pts;
  gop_buffer[buffer_index] = NULL;

  state->global->gop_offset = gop_offset;

  if (++gop_offset >= cfg->gop_len) {
    gop_offset = 0;
    gop_pictures_available = MAX(0, gop_pictures_available - cfg->gop_len);
    gop_buf_read_idx = (gop_buf_read_idx + cfg->gop_len) % gop_buf_size;
  }

  return 1;
}
