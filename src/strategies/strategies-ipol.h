#ifndef STRATEGIES_IPOL_H_
#define STRATEGIES_IPOL_H_
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
 
#include <stdint.h>

#include "encoder.h"


typedef unsigned(ipol_func)(const encoder_control_t * encoder, pixel_t *src, int16_t src_stride, int width, int height, pixel_t *dst,
  int16_t dst_stride, int8_t hor_flag, int8_t ver_flag);

typedef unsigned(epol_func)(int xpos, int ypos, int mv_x, int mv_y, int off_x, int off_y, pixel_t *ref, int ref_width, int ref_height,
  int filterSize, int width, int height, pixel_t *dst);


// Declare function pointers.
extern ipol_func * filter_inter_quarterpel_luma;
extern ipol_func * filter_inter_halfpel_chroma;
extern ipol_func * filter_inter_octpel_chroma;
extern epol_func * extend_borders;


int strategy_register_ipol(void* opaque);


#define STRATEGIES_IPOL_EXPORTS \
  {"filter_inter_quarterpel_luma", (void**) &filter_inter_quarterpel_luma}, \
  {"filter_inter_halfpel_chroma", (void**) &filter_inter_halfpel_chroma}, \
  {"filter_inter_octpel_chroma", (void**) &filter_inter_octpel_chroma}, \
  {"extend_borders", (void**) &extend_borders}, \



#endif //STRATEGIES_IPOL_H_
