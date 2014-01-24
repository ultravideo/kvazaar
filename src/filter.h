#ifndef FILTER_H_
#define FILTER_H_
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
 * \brief Filtering, such as deblocking.
 */

#include "global.h"

#include "encoder.h"


//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
// Deblocking
void filter_deblock_cu(encoder_control *encoder, int32_t x_cu, int32_t y_cu,
                       int8_t depth, int32_t edge);
void filter_deblock_edge_luma(encoder_control *encoder, 
                              int32_t x_pos, int32_t y_pos, 
                              int8_t depth, int8_t dir);
void filter_deblock_edge_chroma(encoder_control *encoder, 
                                int32_t xpos, int32_t ypos, 
                                int8_t depth, int8_t dir);
void filter_deblock(encoder_control *encoder);
void filter_deblock_luma(pixel *src, int32_t offset, int32_t tc , int8_t sw,
                         int8_t part_p_nofilter, int8_t part_q_nofilter,
                         int32_t thr_cut, 
                         int8_t filter_second_p, int8_t filter_second_q);
void filter_deblock_chroma(pixel *src, int32_t offset, int32_t tc,
                           int8_t part_p_nofilter, int8_t part_q_nofilter);

// INTERPOLATION
void filter_inter_halfpel_chroma(int16_t *src, int16_t src_stride, int width, int height,
                                 int16_t *dst, int16_t dst_stride,  int8_t hor_flag, int8_t ver_flag);

// SAO

//////////////////////////////////////////////////////////////////////////
// MACROS
#define EDGE_VER 0
#define EDGE_HOR 1

#endif
