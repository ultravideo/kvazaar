#ifndef CONTEXT_H_
#define CONTEXT_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 * \brief Context derivation for CABAC.
 */

#include "global.h"

#include "encoder.h"
#include "encoderstate.h"
#include "cabac.h"


// Functions
void ctx_init(cabac_ctx* ctx, uint32_t qp, uint32_t init_value);
void init_contexts(encoder_state *encoder_state, int8_t QP, int8_t slice);

void context_copy(encoder_state * target_encoder_state, const encoder_state * source_encoder_state);
int32_t context_calc_pattern_sig_ctx( const uint32_t *sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width);

uint32_t context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,uint32_t pos_x, uint32_t pos_y,int32_t width);


int32_t context_get_sig_ctx_inc(int32_t pattern_sig_ctx,uint32_t scan_idx,int32_t pos_x,
                                int32_t pos_y,int32_t block_type, int8_t texture_type);

#define CNU 154

#endif
