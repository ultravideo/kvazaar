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

/**
 * \brief Struct which contains all picture data
 */
typedef struct image_t
{
  pixel *fulldata;         //!< \brief Allocated buffer (only used in the base_image)

  pixel *y;                //!< \brief Pointer to luma pixel array.
  pixel *u;                //!< \brief Pointer to chroma U pixel array.
  pixel *v;                //!< \brief Pointer to chroma V pixel array.
  pixel *data[NUM_COLORS]; //!< \brief Alternate access method to same data.

  int32_t width;           //!< \brief Luma pixel array width.
  int32_t height;          //!< \brief Luma pixel array height.
  
  int32_t stride;          //!< \brief Luma pixel array width for the full picture (should be used as stride)
  
  struct image_t * base_image; //!< \brief Pointer to the image to which the pixels belong
  int32_t refcount;        //!< \brief Number of references in reflist to the picture
  
  int32_t poc;             //!< \brief Picture order count
} image_t;

typedef struct {
  pixel y[LCU_LUMA_SIZE];
  pixel u[LCU_CHROMA_SIZE];
  pixel v[LCU_CHROMA_SIZE];
} lcu_yuv_t;

typedef struct {
  int size;
  pixel *y;
  pixel *u;
  pixel *v;
} yuv_t;


image_t *image_alloc(const int32_t width, const int32_t height, const int32_t poc);
int image_free(image_t * im);
image_t *image_make_subimage(image_t * const orig_image, const unsigned int x_offset, const unsigned int y_offset, const unsigned int width, const unsigned int height);

yuv_t * yuv_t_alloc(int luma_size);
void yuv_t_free(yuv_t * yuv);

//Algorithms
unsigned image_calc_sad(const image_t *pic, const image_t *ref, int pic_x, int pic_y, int ref_x, int ref_y,
                        int block_width, int block_height, int max_lcu_below);


unsigned pixels_calc_ssd(const pixel *const ref, const pixel *const rec,
                  const int ref_stride, const int rec_stride,
                  const int width);


void pixels_blit(const pixel* orig, pixel *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);

#endif
