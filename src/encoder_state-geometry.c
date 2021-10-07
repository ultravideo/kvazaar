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

#include "encoder_state-geometry.h"

#include "encoder.h"
#include "encoderstate.h"
#include "videoframe.h"


int kvz_lcu_at_slice_start(const encoder_control_t * const encoder, int lcu_addr_in_ts) {
  int i;
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == 0) return 1;
  for (i = 0; i < encoder->slice_count; ++i) {
    if (encoder->slice_addresses_in_ts[i] == lcu_addr_in_ts) return 1;
  }
  return 0;
}

int kvz_lcu_at_slice_end(const encoder_control_t * const encoder, int lcu_addr_in_ts) {
  int i;
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == encoder->in.height_in_lcu * encoder->in.width_in_lcu - 1) return 1;
  for (i = 0; i < encoder->slice_count; ++i) {
    if (encoder->slice_addresses_in_ts[i] == lcu_addr_in_ts + 1) return 1;
  }
  return 0;
}

int kvz_lcu_at_tile_start(const encoder_control_t * const encoder, int lcu_addr_in_ts) {
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == 0) return 1;
  if (encoder->tiles_tile_id[lcu_addr_in_ts - 1] != encoder->tiles_tile_id[lcu_addr_in_ts]) {
    return 1;
  }
  return 0;
}

int kvz_lcu_at_tile_end(const encoder_control_t * const encoder, int lcu_addr_in_ts) {
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == encoder->in.height_in_lcu * encoder->in.width_in_lcu - 1) return 1;
  if (encoder->tiles_tile_id[lcu_addr_in_ts + 1] != encoder->tiles_tile_id[lcu_addr_in_ts]) {
    return 1;
  }
  return 0;
}

//Return 1 if the LCU is at the first row of a structure (tile or slice)
int kvz_lcu_in_first_row(const encoder_state_t * const state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];

  if (lcu_addr_in_rs / state->encoder_control->in.width_in_lcu == state->tile->lcu_offset_y) {
    return 1;
  }

  if (lcu_addr_in_rs / state->encoder_control->in.width_in_lcu == state->slice->start_in_rs / state->encoder_control->in.width_in_lcu) {
    return 1;
  }

  //One row above is before the start of the slice => it's also a boundary
  if (lcu_addr_in_rs - state->encoder_control->in.width_in_lcu < state->slice->start_in_rs) {
    return 1;
  }

  return 0;
}

//Return 1 if the LCU is at the first row of a structure (tile or slice)
int kvz_lcu_in_last_row(const encoder_state_t * const state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];

  if (lcu_addr_in_rs / state->encoder_control->in.width_in_lcu == state->tile->lcu_offset_y + state->tile->frame->height_in_lcu - 1) {
    return 1;
  }

  if (lcu_addr_in_rs / state->encoder_control->in.width_in_lcu == state->slice->end_in_rs / state->encoder_control->in.width_in_lcu) {
    return 1;
  }

  //One row below is before the end of the slice => it's also a boundary
  if (lcu_addr_in_rs + state->encoder_control->in.width_in_lcu > state->slice->end_in_rs) {
    return 1;
  }

  return 0;
}


//Return 1 if the LCU is at the first column of a structure (tile or slice)
int kvz_lcu_in_first_column(const encoder_state_t * const state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];

  //First column of tile?
  if (lcu_addr_in_rs % state->encoder_control->in.width_in_lcu == state->tile->lcu_offset_x) {
    return 1;
  }

  //Slice start may not be aligned with the tile, so we need to allow this
  if (lcu_addr_in_rs == state->slice->start_in_rs) {
    return 1;
  }

  return 0;
}

//Return 1 if the LCU is at the last column of a structure (tile or slice)
int kvz_lcu_in_last_column(const encoder_state_t * const state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];

  //First column of tile?
  if (lcu_addr_in_rs % state->encoder_control->in.width_in_lcu == state->tile->lcu_offset_x + state->tile->frame->width_in_lcu - 1) {
    return 1;
  }

  //Slice start may not be aligned with the tile, so we need to allow this
  if (lcu_addr_in_rs == state->slice->end_in_rs) {
    return 1;
  }

  return 0;
}
