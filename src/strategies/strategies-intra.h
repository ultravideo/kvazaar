#ifndef STRATEGIES_INTRA_H_
#define STRATEGIES_INTRA_H_
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
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Optimization
 * \file
 * Interface for intra prediction functions.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


typedef void (angular_pred_func)(
  const int_fast8_t log2_width,
  const int_fast8_t intra_mode,
  const kvz_pixel *const in_ref_above,
  const kvz_pixel *const in_ref_left,
  kvz_pixel *const dst);

typedef void (intra_pred_planar_func)(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const dst);

typedef void (intra_pred_filtered_dc_func)(
  const int_fast8_t log2_width,
  const kvz_pixel *const ref_top,
  const kvz_pixel *const ref_left,
  kvz_pixel *const out_block);

// Declare function pointers.
extern angular_pred_func * kvz_angular_pred;
extern intra_pred_planar_func * kvz_intra_pred_planar;
extern intra_pred_filtered_dc_func * kvz_intra_pred_filtered_dc;

int kvz_strategy_register_intra(void* opaque, uint8_t bitdepth);


#define STRATEGIES_INTRA_EXPORTS \
  {"angular_pred", (void**) &kvz_angular_pred}, \
  {"intra_pred_planar", (void**) &kvz_intra_pred_planar}, \
  {"intra_pred_filtered_dc", (void**) &kvz_intra_pred_filtered_dc}, \



#endif //STRATEGIES_INTRA_H_
