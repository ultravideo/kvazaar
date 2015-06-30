#ifndef IMAGE_H_
#define IMAGE_H_
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
 * \brief Image and pixel related functions
 */

#include "global.h"
#include "kvazaar.h"

typedef struct {
  pixel_t y[LCU_LUMA_SIZE];
  pixel_t u[LCU_CHROMA_SIZE];
  pixel_t v[LCU_CHROMA_SIZE];
} lcu_yuv_t;

typedef struct {
  int size;
  pixel_t *y;
  pixel_t *u;
  pixel_t *v;
} yuv_t;


kvz_picture *image_alloc(const int32_t width, const int32_t height);

void image_free(kvz_picture *im);

kvz_picture *image_copy_ref(kvz_picture *im);

kvz_picture *image_make_subimage(kvz_picture *const orig_image,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height);

yuv_t * yuv_t_alloc(int luma_size);
void yuv_t_free(yuv_t * yuv);

//Algorithms
unsigned image_calc_sad(const kvz_picture *pic, const kvz_picture *ref, int pic_x, int pic_y, int ref_x, int ref_y,
                        int block_width, int block_height, int max_lcu_below);


unsigned pixels_calc_ssd(const pixel_t *const ref, const pixel_t *const rec,
                  const int ref_stride, const int rec_stride,
                  const int width);


void pixels_blit(const pixel_t* orig, pixel_t *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);

#endif
