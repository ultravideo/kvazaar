#ifndef INPUT_FRAME_BUFFER_H_
#define INPUT_FRAME_BUFFER_H_
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

#include "global.h"

// Forward declaration.
struct encoder_state_t;

typedef struct input_frame_buffer_t {
  /** \brief An array for stroring the input frames. */
  struct kvz_picture *pic_buffer[2 * KVZ_MAX_GOP_LENGTH];

  /** \brief Number of pictures in the buffer. */
  int pictures_available;

  /** \brief Index where the next input frame is put to. */
  int write_idx;

  /** \brief Index of the first frame of the current GOP. */
  int read_idx;

  /** \brief Number of the next frame in the current GOP. */
  int gop_offset;

} input_frame_buffer_t;

void kvz_init_input_frame_buffer(input_frame_buffer_t *input_buffer);

int kvz_encoder_feed_frame(input_frame_buffer_t *buf,
                           struct encoder_state_t *const state,
                           struct kvz_picture *const img_in);

#endif // INPUT_FRAME_BUFFER_H_
