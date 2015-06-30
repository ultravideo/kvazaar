#ifndef FILTER_H_
#define FILTER_H_
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
 * \brief Filtering, such as deblocking.
 */

#include "global.h"

#include "encoder.h"
#include "encoderstate.h"


//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
// Deblocking
void filter_deblock_cu(encoder_state_t *state, int32_t x_px, int32_t y_px,
                       int8_t depth, int32_t edge);
void filter_deblock_edge_luma(encoder_state_t *state,
                              int32_t x_pos, int32_t y_pos,
                              int8_t depth, int8_t dir);
void filter_deblock_edge_chroma(encoder_state_t *state,
                                int32_t xpos, int32_t ypos,
                                int8_t depth, int8_t dir);
void filter_deblock_lcu(encoder_state_t *state, int x_px, int y_px);
void filter_deblock_luma(const encoder_control_t * const encoder, kvz_pixel *src, int32_t offset, int32_t tc , int8_t sw,
                         int8_t part_p_nofilter, int8_t part_q_nofilter,
                         int32_t thr_cut,
                         int8_t filter_second_p, int8_t filter_second_q);
void filter_deblock_chroma(const encoder_control_t * const encoder, kvz_pixel *src, int32_t offset, int32_t tc,
                           int8_t part_p_nofilter, int8_t part_q_nofilter);

// SAO

//////////////////////////////////////////////////////////////////////////
// MACROS
#define EDGE_VER 0
#define EDGE_HOR 1

#endif
