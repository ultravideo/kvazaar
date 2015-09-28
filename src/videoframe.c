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

/*
 * \file
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sao.h"
#include "threads.h"
#include "videoframe.h"

/**
 * \brief Allocate new frame
 * \param pic picture pointer
 * \return picture pointer
 */
videoframe_t *kvz_videoframe_alloc(const int32_t width, const int32_t height, const int32_t poc) {
  videoframe_t *frame = MALLOC(videoframe_t, 1);

  if (!frame) return 0;

  FILL(*frame, 0);

  frame->width  = width;
  frame->height = height;
  frame->width_in_lcu  = frame->width / LCU_WIDTH;
  if (frame->width_in_lcu * LCU_WIDTH < frame->width) frame->width_in_lcu++;
  frame->height_in_lcu = frame->height / LCU_WIDTH;
  if (frame->height_in_lcu * LCU_WIDTH < frame->height) frame->height_in_lcu++;

  {
    // Allocate height_in_scu x width_in_scu x sizeof(CU_info)
    unsigned height_in_scu = frame->height_in_lcu << MAX_DEPTH;
    unsigned width_in_scu = frame->width_in_lcu << MAX_DEPTH;
    frame->cu_array = kvz_cu_array_alloc(width_in_scu, height_in_scu);
  }

  frame->coeff_y = NULL; frame->coeff_u = NULL; frame->coeff_v = NULL;

  frame->sao_luma = MALLOC(sao_info_t, frame->width_in_lcu * frame->height_in_lcu);
  frame->sao_chroma = MALLOC(sao_info_t, frame->width_in_lcu * frame->height_in_lcu);

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

  kvz_cu_array_free(frame->cu_array);

  FREE_POINTER(frame->coeff_y);
  FREE_POINTER(frame->coeff_u);
  FREE_POINTER(frame->coeff_v);

  FREE_POINTER(frame->sao_luma);
  FREE_POINTER(frame->sao_chroma);

  free(frame);

  return 1;
}

void kvz_videoframe_set_poc(videoframe_t * const frame, const int32_t poc) {
  frame->poc = poc;
}

const cu_info_t* kvz_videoframe_get_cu_const(const videoframe_t * const frame, unsigned int x_in_scu, unsigned int y_in_scu)
{
  assert(x_in_scu < (frame->width_in_lcu << MAX_DEPTH));
  assert(y_in_scu < (frame->height_in_lcu << MAX_DEPTH));
  
  return &frame->cu_array->data[x_in_scu + y_in_scu * (frame->width_in_lcu << MAX_DEPTH)];
}

cu_info_t* kvz_videoframe_get_cu(videoframe_t * const frame, const unsigned int x_in_scu, const unsigned int y_in_scu)
{
  assert(x_in_scu < (frame->width_in_lcu << MAX_DEPTH));
  assert(y_in_scu < (frame->height_in_lcu << MAX_DEPTH));
  
  return &frame->cu_array->data[x_in_scu + y_in_scu * (frame->width_in_lcu << MAX_DEPTH)];
}

#if KVZ_BIT_DEPTH == 8
#define PSNRMAX (255.0 * 255.0)
#else
  #define PSNRMAX ((double)PIXEL_MAX * (double)PIXEL_MAX)
#endif

/**
 * \brief Calculates image PSNR value
 *
 * \param src   source picture
 * \param rec   reconstructed picture
 * \prama psnr  returns the PSNR
 */
void kvz_videoframe_compute_psnr(const kvz_picture *const src,
                                 const kvz_picture *const rec,
                                 double psnr[NUM_COLORS])
{
  assert(src->width  == rec->width);
  assert(src->height == rec->height);

  int32_t pixels = src->width * src->height;

  for (int32_t c = 0; c < NUM_COLORS; ++c) {
    int32_t num_pixels = pixels;
    if (c != COLOR_Y) {
      num_pixels >>= 2;
    }
    psnr[c] = 0;
    for (int32_t i = 0; i < num_pixels; ++i) {
      const int32_t error = src->data[c][i] - rec->data[c][i];
      psnr[c] += error * error;
    }

    // Avoid division by zero
    if (psnr[c] == 0) psnr[c] = 99.0;
    psnr[c] = 10 * log10((num_pixels * PSNRMAX) / ((double)psnr[c]));;
  }
}
