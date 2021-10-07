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

#include "videoframe.h"

#include <stdlib.h>

#include "image.h"
#include "sao.h"


/**
 * \brief Allocate new frame
 * \param pic picture pointer
 * \return picture pointer
 */
videoframe_t * kvz_videoframe_alloc(int32_t width,
                                    int32_t height,
                                    enum kvz_chroma_format chroma_format)
{
  videoframe_t *frame = calloc(1, sizeof(videoframe_t));
  if (!frame) return 0;

  frame->width  = width;
  frame->height = height;
  frame->width_in_lcu  = CEILDIV(frame->width,  LCU_WIDTH);
  frame->height_in_lcu = CEILDIV(frame->height, LCU_WIDTH);

  frame->sao_luma = MALLOC(sao_info_t, frame->width_in_lcu * frame->height_in_lcu);
  if (chroma_format != KVZ_CSP_400) {
    frame->sao_chroma = MALLOC(sao_info_t, frame->width_in_lcu * frame->height_in_lcu);
  }

  return frame;
}

/**
 * \brief Free memory allocated to frame
 * \param pic picture pointer
 * \return 1 on success, 0 on failure
 */
int kvz_videoframe_free(videoframe_t * const frame)
{
  kvz_image_free(frame->source);
  frame->source = NULL;
  kvz_image_free(frame->rec);
  frame->rec = NULL;

  kvz_cu_array_free(&frame->cu_array);

  FREE_POINTER(frame->sao_luma);
  FREE_POINTER(frame->sao_chroma);

  free(frame);

  return 1;
}

void kvz_videoframe_set_poc(videoframe_t * const frame, const int32_t poc) {
  frame->poc = poc;
}
