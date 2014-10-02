#ifndef STRATEGIES_DCT_H_
#define STRATEGIES_DCT_H_
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
#include <stdint.h>

typedef unsigned (dct_func)(int8_t bitdepth, const int16_t *input, int16_t *output);


// Declare function pointers.
extern dct_func * fast_forward_dst_4x4;

extern dct_func * dct_4x4;
extern dct_func * dct_8x8;
extern dct_func * dct_16x16;
extern dct_func * dct_32x32;

extern dct_func * fast_inverse_dst_4x4;

extern dct_func * idct_4x4;
extern dct_func * idct_8x8;
extern dct_func * idct_16x16;
extern dct_func * idct_32x32;


int strategy_register_dct(void* opaque);
dct_func * get_dct_func(int8_t width, int32_t mode);
dct_func * get_idct_func(int8_t width, int32_t mode);


#define STRATEGIES_DCT_EXPORTS \
  {"fast_forward_dst_4x4", (void**) &fast_forward_dst_4x4}, \
  \
  {"dct_4x4", (void**) &dct_4x4}, \
  {"dct_8x8", (void**) &dct_8x8}, \
  {"dct_16x16", (void**) &dct_16x16}, \
  {"dct_32x32", (void**) &dct_32x32}, \
  \
  {"fast_inverse_dst_4x4", (void**) &fast_inverse_dst_4x4}, \
  \
  {"idct_4x4", (void**)&idct_4x4}, \
  {"idct_8x8", (void**)&idct_8x8}, \
  {"idct_16x16", (void**)&idct_16x16}, \
  {"idct_32x32", (void**)&idct_32x32}, \



#endif //STRATEGIES_DCT_H_
