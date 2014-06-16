#ifndef STRATEGIES_PICTURE_H_
#define STRATEGIES_PICTURE_H_
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

#include "../image.h"


typedef unsigned(reg_sad_func)(const pixel *const data1, const pixel *const data2,
  const int width, const int height,
  const unsigned stride1, const unsigned stride2);
typedef unsigned (cost_pixel_nxn_func)(const pixel *block1, const pixel *block2);


// Declare function pointers.
extern reg_sad_func * reg_sad;

extern cost_pixel_nxn_func * sad_8bit_4x4;
extern cost_pixel_nxn_func * sad_8bit_8x8;
extern cost_pixel_nxn_func * sad_8bit_16x16;
extern cost_pixel_nxn_func * sad_8bit_32x32;
extern cost_pixel_nxn_func * sad_8bit_64x64;

extern cost_pixel_nxn_func * satd_8bit_4x4;
extern cost_pixel_nxn_func * satd_8bit_8x8;
extern cost_pixel_nxn_func * satd_8bit_16x16;
extern cost_pixel_nxn_func * satd_8bit_32x32;
extern cost_pixel_nxn_func * satd_8bit_64x64;


int strategy_register_picture(void* opaque);
cost_pixel_nxn_func * pixels_get_satd_func(unsigned n);
cost_pixel_nxn_func * pixels_get_sad_func(unsigned n);


#define STRATEGIES_PICTURE_EXPORTS \
  {"reg_sad", (void**) &reg_sad}, \
  {"sad_8bit_4x4", (void**) &sad_8bit_4x4}, \
  {"sad_8bit_8x8", (void**) &sad_8bit_8x8}, \
  {"sad_8bit_16x16", (void**) &sad_8bit_16x16}, \
  {"sad_8bit_32x32", (void**) &sad_8bit_32x32}, \
  {"sad_8bit_64x64", (void**) &sad_8bit_64x64}, \
  {"satd_8bit_4x4", (void**) &satd_8bit_4x4}, \
  {"satd_8bit_8x8", (void**) &satd_8bit_8x8}, \
  {"satd_8bit_16x16", (void**) &satd_8bit_16x16}, \
  {"satd_8bit_32x32", (void**) &satd_8bit_32x32}, \
  {"satd_8bit_64x64", (void**) &satd_8bit_64x64}, \



#endif //STRATEGIES_PICTURE_H_
