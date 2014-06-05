#ifndef IMAGE_H_
#define IMAGE_H_
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

/*
 * \file
 * \brief Image related functions
 */

#include "global.h"

/**
 * \brief Struct which contains all picture data
 */
typedef struct image
{
  pixel *fulldata;         //!< \brief Allocated buffer (only used in the base_image)

  pixel *y;                //!< \brief Pointer to luma pixel array.
  pixel *u;                //!< \brief Pointer to chroma U pixel array.
  pixel *v;                //!< \brief Pointer to chroma V pixel array.
  pixel *data[NUM_COLORS]; //!< \brief Alternate access method to same data.

  int32_t width;           //!< \brief Luma pixel array width.
  int32_t height;          //!< \brief Luma pixel array height.
  
  int32_t stride;          //!< \brief Luma pixel array width for the full picture (should be used as stride)
  
  struct image * base_image; //!< \brief Pointer to the image to which the pixels belong
  int32_t refcount;        //!< \brief Number of references in reflist to the picture
} image;

image *image_alloc(const int32_t width, const int32_t height);
int image_free(image * im);

#endif
