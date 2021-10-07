#ifndef ENCODER_STATE_GEOMETRY_H_
#define ENCODER_STATE_GEOMETRY_H_
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
 * Helper functions for tiles and slices.
 */

#include "global.h" // IWYU pragma: keep


// Forward declare because including the header would lead  to a cyclic
// dependency.
struct encoder_control_t;
struct encoder_state_t;


int kvz_lcu_at_slice_start(const struct  encoder_control_t * const encoder, int lcu_addr_in_ts);
int kvz_lcu_at_slice_end(const struct  encoder_control_t * const encoder, int lcu_addr_in_ts);
int kvz_lcu_at_tile_start(const struct  encoder_control_t * const encoder, int lcu_addr_in_ts);
int kvz_lcu_at_tile_end(const struct  encoder_control_t * const encoder, int lcu_addr_in_ts);
int kvz_lcu_in_first_row(const struct encoder_state_t * const encoder_state, int lcu_addr_in_ts);
int kvz_lcu_in_last_row(const struct encoder_state_t * const encoder_state, int lcu_addr_in_ts);
int kvz_lcu_in_first_column(const struct encoder_state_t * const encoder_state, int lcu_addr_in_ts);
int kvz_lcu_in_last_column(const struct encoder_state_t * const encoder_state, int lcu_addr_in_ts);


#endif // ENCODER_STATE_GEOMETRY_H_
