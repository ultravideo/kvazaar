#ifndef KVAZAAR_INTERNAL_H_
#define KVAZAAR_INTERNAL_H_
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
 * \brief Definitions for opaque structs in kvazaar.h
 */

#include "global.h" // IWYU pragma: keep

#include "kvazaar.h"
#include "input_frame_buffer.h"


// Forward declarations.
struct encoder_state_t;
struct encoder_control_t;

struct kvz_encoder {
  const struct encoder_control_t* control;
  struct encoder_state_t* states;
  unsigned num_encoder_states;

  /**
   * \brief Number of the current encoder state.
   *
   * The current state is the one that will be used for encoding the frame
   * that is started next.
   */
  unsigned cur_state_num;

  /**
   * \brief Number of the next encoder state to be finished.
   */
  unsigned out_state_num;

  /**
   * \brief Buffer for input frames.
   */
  input_frame_buffer_t input_buffer;

  unsigned frames_started;
  unsigned frames_done;
};

#endif // KVAZAAR_INTERNAL_H_
