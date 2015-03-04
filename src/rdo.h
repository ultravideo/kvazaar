#ifndef RDO_H_
#define RDO_H_
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
 * \brief Handling Rate-Distortion Optimization related functionality
 */

#include "global.h"

#include "encoder.h"
#include "encoderstate.h"


extern const uint32_t g_go_rice_range[5];
extern const uint32_t g_go_rice_prefix_len[5];

int intra_rdo_cost_compare(uint32_t *rdo_costs,int8_t rdo_modes_to_check, uint32_t cost);

void  rdoq(encoder_state_t *encoder_state, coeff_t *coef, coeff_t *dest_coeff, int32_t width,
           int32_t height, int8_t type, int8_t scan_mode, int8_t block_type, int8_t tr_depth);

uint32_t rdo_cost_intra(encoder_state_t *encoder, pixel* pred, pixel* orig_block, int width, int8_t mode, int tr_depth);

int32_t get_coeff_cost(const encoder_state_t *encoder_state, coeff_t *coeff, int32_t width, int32_t type, int8_t scan_mode);

int32_t get_ic_rate(encoder_state_t *encoder_state, uint32_t abs_level, uint16_t ctx_num_one,uint16_t ctx_num_abs,
                     uint16_t abs_go_rice, uint32_t c1_idx, uint32_t c2_idx, int8_t type);
uint32_t get_coded_level ( encoder_state_t * encoder_state, double* coded_cost, double* coded_cost0, double* coded_cost_sig,
                           int32_t level_double, uint32_t max_abs_level,
                           uint16_t ctx_num_sig, uint16_t ctx_num_one, uint16_t ctx_num_abs,
                           uint16_t abs_go_rice,
                           uint32_t c1_idx, uint32_t c2_idx,
                           int32_t q_bits,double temp, int8_t last, int8_t type);

extern const float f_entropy_bits[128];
#define CTX_ENTROPY_FBITS(ctx,val) f_entropy_bits[(ctx)->uc_state ^ (val)]

#endif
