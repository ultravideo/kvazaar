#ifndef VIDEOFRAME_H_
#define VIDEOFRAME_H_
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
 * \ingroup DataStructures
 * \file
 * \brief Container for the frame currently being encoded.
 */

#include "cu.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


/**
 * \brief Struct which contains all picture data
 */
typedef struct videoframe
{
  kvz_picture *source;         //!< \brief Source image.
  kvz_picture *rec;            //!< \brief Reconstructed image.

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.

  cu_array_t* cu_array;     //!< \brief Info for each CU at each depth.
  struct sao_info_t *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_t *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  int32_t poc;           //!< \brief Picture order count
} videoframe_t;


videoframe_t *kvz_videoframe_alloc(int32_t width, int32_t height, enum kvz_chroma_format chroma_format);
int kvz_videoframe_free(videoframe_t * const frame);

void kvz_videoframe_set_poc(videoframe_t * frame, int32_t poc);

#endif
