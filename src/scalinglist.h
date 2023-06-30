#ifndef SCALINGLIST_H_
#define SCALINGLIST_H_
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
 * \ingroup Reconstruction
 * \file
 * Scaling list initialization.
 */

#include <stdio.h>

#include "global.h" // IWYU pragma: keep


typedef struct {
        int8_t   enable;
        int8_t   use_default_list;
        int32_t  scaling_list_dc   [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
  const coeff_t *scaling_list_coeff[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
  const coeff_t *quant_coeff[4][6][6];
  const coeff_t *de_quant_coeff  [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM];
  const double *error_scale[4][6][6];
} scaling_list_t;

extern const uint8_t kvz_g_scaling_list_num[4];
extern const uint16_t kvz_g_scaling_list_size[4];

const coeff_t*kvz_scalinglist_get_default(const uint32_t size_id, const uint32_t list_id);

void kvz_scalinglist_init(scaling_list_t * const scaling_list);
void kvz_scalinglist_destroy(scaling_list_t * const scaling_list);

int  kvz_scalinglist_parse(scaling_list_t * const scaling_list, FILE *fp);
void kvz_scalinglist_process(scaling_list_t *scaling_list, uint8_t bitdepth);






#endif
