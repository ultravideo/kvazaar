#ifndef STRATEGIES_DCT_H_
#define STRATEGIES_DCT_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Optimization
 * \file
 * Interface for transform functions.
 */

#include "global.h" // IWYU pragma: keep
#include "cu.h"

typedef void (dct_func)(int8_t bitdepth, const int16_t *input, int16_t *output);


// Declare function pointers.
extern dct_func * kvz_fast_forward_dst_4x4;

extern dct_func * kvz_dct_4x4;
extern dct_func * kvz_dct_8x8;
extern dct_func * kvz_dct_16x16;
extern dct_func * kvz_dct_32x32;

extern dct_func * kvz_fast_inverse_dst_4x4;

extern dct_func * kvz_idct_4x4;
extern dct_func * kvz_idct_8x8;
extern dct_func * kvz_idct_16x16;
extern dct_func * kvz_idct_32x32;


int kvz_strategy_register_dct(void* opaque, uint8_t bitdepth);
dct_func * kvz_get_dct_func(int8_t width, color_t color, cu_type_t type);
dct_func * kvz_get_idct_func(int8_t width, color_t color, cu_type_t type);



#define STRATEGIES_DCT_EXPORTS \
  {"fast_forward_dst_4x4", (void**) &kvz_fast_forward_dst_4x4}, \
  \
  {"dct_4x4", (void**) &kvz_dct_4x4}, \
  {"dct_8x8", (void**) &kvz_dct_8x8}, \
  {"dct_16x16", (void**) &kvz_dct_16x16}, \
  {"dct_32x32", (void**) &kvz_dct_32x32}, \
  \
  {"fast_inverse_dst_4x4", (void**) &kvz_fast_inverse_dst_4x4}, \
  \
  {"idct_4x4", (void**)&kvz_idct_4x4}, \
  {"idct_8x8", (void**)&kvz_idct_8x8}, \
  {"idct_16x16", (void**)&kvz_idct_16x16}, \
  {"idct_32x32", (void**)&kvz_idct_32x32}, \



#endif //STRATEGIES_DCT_H_
