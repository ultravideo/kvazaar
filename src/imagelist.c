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
#include "imagelist.h"
#include "strategyselector.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


/**
 * \brief Allocate memory for image_list
 * \param size  initial array size
 * \return image_list pointer, NULL on failure
 */
image_list * image_list_alloc(int size)
{
  image_list *list = (image_list *)malloc(sizeof(image_list));
  list->size = size;
  if (size > 0) {
    list->images = (image**)malloc(sizeof(image*) * size);
    list->cu_arrays = (cu_info**)malloc(sizeof(cu_info*) * size);
  }

  list->used_size = 0;

  return list;
}

/**
 * \brief Resize image_list array
 * \param list  image_list pointer
 * \param size  new array size
 * \return 1 on success, 0 on failure
 */
int image_list_resize(image_list *list, unsigned size)
{
  unsigned int i;
  image** old_images = NULL;
  cu_info** old_cu_arrays = NULL;
  
  //FIXME This could be done in a simple way using realloc...

  // No need to do anything when resizing to same size
  if (size == list->size) {
    return 1;
  }

  // Save the old list
  if (list->used_size > 0) {
    old_images = list->images;
    old_cu_arrays = list->cu_arrays;
  }

  // allocate space for the new list
  list->images = (image**)malloc(sizeof(image*)*size);
  list->cu_arrays = (cu_info**)malloc(sizeof(cu_info*)*size);

  // Copy everything from the old list to the new if needed.
  if (old_images != NULL) {
    for (i = 0; i < list->used_size; ++i) {
      list->images[i] = old_images[i];
      list->cu_arrays[i] = old_cu_arrays[i];
    }

    free(old_images);
    free(old_cu_arrays);
  }

  return 1;
}

/**
 * \brief Free memory allocated to the picture_list
 * \param list image_list pointer
 * \return 1 on success, 0 on failure
 */
int image_list_destroy(image_list *list)
{
  unsigned int i;
  if (list->used_size > 0) {
    for (i = 0; i < list->used_size; ++i) {
      image_free(list->images[i]);
      free(list->cu_arrays[i]);
      list->images[i] = NULL;
    }
  }

  if (list->size > 0) {
    free(list->images);
    free(list->cu_arrays);
  }
  free(list);
  return 1;
}

/**
 * \brief Add picture to the front of the picturelist
 * \param pic picture pointer to add
 * \param picture_list list to use
 * \return 1 on success
 */
int image_list_add(image_list *list, image* im, cu_info* cu_array)
{
  int i = 0;
  if (ATOMIC_INC(&(im->refcount)) == 1) {
    fprintf(stderr, "Tried to add an unreferenced picture. This is a bug!\n");
    assert(0); //Stop for debugging
    return 0;
  }

  if (list->size == list->used_size) {
    if (!image_list_resize(list, list->size*2)) return 0;
  }

  for (i = list->used_size; i > 0; i--) {
    list->images[i] = list->images[i - 1];
    list->cu_arrays[i] = list->cu_arrays[i - 1];
  }

  list->images[0] = im;
  //We need (only here, for malloc/memcpy) to compute the size of the image in SCU
  {
    //FIXME FIXME FIXME Do like images, use a pointer instead of copying
    unsigned int width_in_lcu, height_in_lcu, width_in_scu, height_in_scu;
    width_in_lcu  = im->width / LCU_WIDTH;
    if (width_in_lcu * LCU_WIDTH < im->width) width_in_lcu++;
    height_in_lcu = im->height / LCU_WIDTH;
    if (height_in_lcu * LCU_WIDTH < im->height) height_in_lcu++;
    height_in_scu = height_in_lcu << MAX_DEPTH;
    width_in_scu = width_in_lcu << MAX_DEPTH;
    
    list->cu_arrays[0] = (cu_info*)malloc(sizeof(cu_info) * width_in_scu * height_in_scu);
    memcpy(list->cu_arrays[0], cu_array, sizeof(cu_info) * width_in_scu * height_in_scu);
  }
  list->used_size++;
  return 1;
}

/**
 * \brief Remove picture from picturelist
 * \param list list to use
 * \param n index to remove
 * \return 1 on success
 */
int image_list_rem(image_list * const list, const unsigned n)
{
  // Must be within list boundaries
  if (n >= list->used_size)
  {
    return 0;
  }

  if (!image_free(list->images[n])) {
    fprintf(stderr, "Could not free image!\n");
    assert(0); //Stop here
    return 0;
  }
  free(list->cu_arrays[n]);

  // The last item is easy to remove
  if (n == list->used_size - 1) {
    list->images[n] = NULL;
    list->used_size--;
  } else {
    int i = n;
    // Shift all following pics one backward in the list
    for (i = n; i < list->used_size - 1; ++i) {
      list->images[i] = list->images[i + 1];
    }
    list->images[list->used_size - 1] = NULL;
    list->used_size--;
  }

  return 1;
}
