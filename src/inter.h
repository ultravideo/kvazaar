#ifndef INTER_H_
#define INTER_H_
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
 * \brief Handling Coding Units (CU's) for inter frames.
 */

#include "global.h"

#include "image.h"
#include "encoder.h"
#include "encoderstate.h"

typedef struct {
  uint8_t dir;
  uint8_t ref[2];
  int16_t mv[2][2];

} inter_merge_cand_t;


//void kvz_inter_set_block(image* im,uint32_t x_cu, uint32_t y_cu, uint8_t depth, cu_info *cur_cu);
void kvz_inter_recon_lcu(const encoder_state_t * const state, const kvz_picture * ref, int32_t xpos, int32_t ypos, int32_t width, const int16_t mv_param[2], lcu_t* lcu, hi_prec_buf_t *hi_prec_out);
void kvz_inter_recon_lcu_bipred(const encoder_state_t * const state, const kvz_picture * ref1, const kvz_picture * ref2, int32_t xpos, int32_t ypos, int32_t width, int16_t mv_param[2][2], lcu_t* lcu);

void kvz_inter_get_spatial_merge_candidates(int32_t x, int32_t y, int8_t depth, cu_info_t **b0, cu_info_t **b1,
                                        cu_info_t **b2, cu_info_t **a0, cu_info_t **a1, lcu_t *lcu);
void kvz_inter_get_mv_cand(const encoder_state_t *state, int32_t x, int32_t y, int8_t depth, int16_t mv_cand[2][2], cu_info_t* cur_cu, lcu_t *lcu, int8_t reflist);
uint8_t kvz_inter_get_merge_cand(const encoder_state_t *state, int32_t x, int32_t y, int8_t depth, inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS], lcu_t *lcu);
#endif
