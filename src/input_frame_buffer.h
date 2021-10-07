#ifndef INPUT_FRAME_BUFFER_H_
#define INPUT_FRAME_BUFFER_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Control
 * \file
 * Buffering of input for reordering.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


// Forward declaration.
struct encoder_state_t;

typedef struct input_frame_buffer_t {
  /** \brief An array for stroring the input frames. */
  struct kvz_picture *pic_buffer[3 * KVZ_MAX_GOP_LENGTH];

  /** \brief An array for stroring the timestamps. */
  int64_t pts_buffer[3 * KVZ_MAX_GOP_LENGTH];

  /** \brief Number of pictures input. */
  uint64_t num_in;

  /** \brief Number of pictures output. */
  uint64_t num_out;

  /** \brief Value to subtract from the DTS values of the first frames.
   *
   * This will be set to the difference of the PTS values of the first and
   * (cfg->gop_len)th frames, unless the sequence has less that cfg->gop_len
   * frames.
   */
  int64_t delay;

  /** \brief Number of GOP pictures skipped.
   *
   * This is used when the last GOP of the sequence is not full.
   */
  int gop_skipped;

} input_frame_buffer_t;

void kvz_init_input_frame_buffer(input_frame_buffer_t *input_buffer);

kvz_picture* kvz_encoder_feed_frame(input_frame_buffer_t *buf,
                                    struct encoder_state_t *const state,
                                    struct kvz_picture *const img_in,
                                    int first_done);

#endif // INPUT_FRAME_BUFFER_H_
