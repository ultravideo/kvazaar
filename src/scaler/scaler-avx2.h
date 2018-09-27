#ifndef SCALER_AVX2_H_
#define SCALER_AVX2_H_
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
#include "scaler.h"

//#define MIN(x,y) (((x) < (y)) ? (x) : (y))
//#define MAX(x,y) (((x) > (y)) ? (x) : (y))
//
//#define SHIFT(x,y) (((y) < 0) ? ((x)>>(-(y))) : ((x)<<(y)))


//void resample_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma);
//void _resample_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma);

extern resample_block_step_func *const kvz_default_block_step_resample_func_avx2;
extern resample_block_step_func *const kvz_alt1_block_step_resample_func_avx2;
extern resample_block_step_func *const kvz_alt2_block_step_resample_func_avx2;
extern resample_func *const kvz_default_resample_func_avx2;
extern resample_func *const kvz_alt_resample_func_avx2;

//int test_avx();

#endif