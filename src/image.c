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
 */

#include "threads.h"
#include "image.h"
#include "strategyselector.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "sao.h"

/**
 * \brief Allocate new image
 * \return image pointer
 */
image *image_alloc(const int32_t width, const int32_t height)
{
  image *im = MALLOC(image, 1);
  
  unsigned int luma_size = width * height;
  unsigned int chroma_size = luma_size / 4;
  
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  if (!im) return NULL;
  
  im->width = width;
  im->height = height;
  im->stride = width;
  
  im->base_image = im;
  
  im->refcount = 1; //We give a reference to caller
  
  //Allocate memory
  im->fulldata = MALLOC(pixel, (luma_size + 2*chroma_size));
  im->y = im->data[COLOR_Y] = &im->fulldata[0];
  im->u = im->data[COLOR_U] = &im->fulldata[luma_size];
  im->v = im->data[COLOR_V] = &im->fulldata[luma_size + chroma_size];

  return im;
}

/**
 * \brief Free memory allocated to picture (if we have no reference left)
 * \param pic picture pointer
 * \return 1 on success, 0 on failure
 */
int image_free(image * const im)
{
  //Either we are the base image, or we should have no references
  assert(im->base_image == im || im->refcount == 0);
  
  int32_t new_refcount = ATOMIC_DEC(&(im->base_image->refcount));
  if (new_refcount > 0) return 1;
  FREE_POINTER(im->fulldata);
  
  //Just to make the program crash when using those values after the free
  im->y = im->u = im->v = im->data[COLOR_Y] = im->data[COLOR_U] = im->data[COLOR_V] = NULL;
  
  free(im);

  return 1;
}


image *image_make_subimage(image * const orig_image, const unsigned int x_offset, const unsigned int y_offset, const unsigned int width, const unsigned int height)
{
  image *im = MALLOC(image, 1);
  if (!im) return NULL;
  
  im->base_image = orig_image->base_image;
  ATOMIC_INC(&(im->base_image->refcount));
  
  assert(x_offset + width <= orig_image->width);
  assert(y_offset + height <= orig_image->height);
  
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);
  
  assert((x_offset % 2) == 0);
  assert((y_offset % 2) == 0);
  
  im->stride = orig_image->stride;
  im->refcount = 0; //No references on subimages
  
  im->width = width;
  im->height = height;
  
  im->y = im->data[COLOR_Y] = &orig_image->y[x_offset + y_offset * orig_image->stride];
  im->u = im->data[COLOR_U] = &orig_image->u[x_offset/2 + y_offset/2 * orig_image->stride/2];
  im->v = im->data[COLOR_V] = &orig_image->v[x_offset/2 + y_offset/2 * orig_image->stride/2];

  return im;
}