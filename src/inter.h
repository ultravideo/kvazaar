#ifndef INTER_H_
#define INTER_H_
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
 * \brief Handling Coding Units (CU's) for inter frames.
 */

#include "global.h"

#include "picture.h"
#include "encoder.h"


void inter_set_block(picture* pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, cu_info *cur_cu);
void inter_recon_lcu(const encoder_control *encoder, picture* ref,int32_t xpos, int32_t ypos,int32_t width, const int16_t mv_param[2], lcu_t *lcu);

void inter_get_spatial_merge_candidates(int32_t x, int32_t y, int8_t depth, cu_info **b0, cu_info **b1,
                                        cu_info **b2,cu_info **a0,cu_info **a1, lcu_t *lcu);
void inter_get_mv_cand(const encoder_state *encoder_state, int32_t x, int32_t y, int8_t depth, int16_t mv_cand[2][2], cu_info* cur_cu, lcu_t *lcu);
uint8_t inter_get_merge_cand(int32_t x, int32_t y, int8_t depth, int16_t mv_cand[MRG_MAX_NUM_CANDS][3], lcu_t *lcu);
#endif
