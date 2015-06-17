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

#include "threads.h"
#include "imagelist.h"
#include "strategyselector.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


/**
 * \brief Allocate memory for image_list
 * \param size  initial array size
 * \return image_list pointer, NULL on failure
 */
image_list_t * image_list_alloc(int size)
{
  image_list_t *list = (image_list_t *)malloc(sizeof(image_list_t));
  list->size = size;
  if (size > 0) {
    list->images = (image_t**)malloc(sizeof(image_t*) * size);
    list->cu_arrays = (cu_array_t**)malloc(sizeof(cu_array_t*) * size);
    list->pocs = (int32_t*)malloc(sizeof(int32_t) * size);
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
int image_list_resize(image_list_t *list, unsigned size)
{
  list->images = (image_t**)realloc(list->images, sizeof(image_t*) * size);
  list->cu_arrays = (cu_array_t**)realloc(list->cu_arrays, sizeof(cu_array_t*) * size);
  list->pocs = (int32_t*)realloc(list->pocs, sizeof(int32_t*) * size);
  list->size = size;
  return size == 0 || (list->images && list->cu_arrays && list->pocs);
}

/**
 * \brief Free memory allocated to the picture_list
 * \param list image_list pointer
 * \return 1 on success, 0 on failure
 */
int image_list_destroy(image_list_t *list)
{
  unsigned int i;
  if (list->used_size > 0) {
    for (i = 0; i < list->used_size; ++i) {
      image_free(list->images[i]);
      list->images[i] = NULL;
      cu_array_free(list->cu_arrays[i]);
      list->cu_arrays[i] = NULL;
      list->pocs[i] = 0;
    }
  }

  if (list->size > 0) {
    free(list->images);
    free(list->cu_arrays);
    free(list->pocs);
  }
  list->images = NULL;
  list->cu_arrays = NULL;
  list->pocs = NULL;
  free(list);
  return 1;
}

/**
 * \brief Add picture to the front of the picturelist
 * \param pic picture pointer to add
 * \param picture_list list to use
 * \return 1 on success
 */
int image_list_add(image_list_t *list, image_t* im, cu_array_t* cua, int32_t poc)
{
  int i = 0;
  if (ATOMIC_INC(&(im->refcount)) == 1) {
    fprintf(stderr, "Tried to add an unreferenced picture. This is a bug!\n");
    assert(0); //Stop for debugging
    return 0;
  }
  
  if (ATOMIC_INC(&(cua->refcount)) == 1) {
    fprintf(stderr, "Tried to add an unreferenced cu_array. This is a bug!\n");
    assert(0); //Stop for debugging
    return 0;
  }

  if (list->size == list->used_size) {
    if (!image_list_resize(list, list->size*2)) return 0;
  }

  for (i = list->used_size; i > 0; i--) {
    list->images[i] = list->images[i - 1];
    list->cu_arrays[i] = list->cu_arrays[i - 1];
    list->pocs[i] = list->pocs[i - 1];
  }

  list->images[0] = im;
  list->cu_arrays[0] = cua;
  list->pocs[0] = poc;
  
  list->used_size++;
  return 1;
}

/**
 * \brief Remove picture from picturelist
 * \param list list to use
 * \param n index to remove
 * \return 1 on success
 */
int image_list_rem(image_list_t * const list, const unsigned n)
{
  // Must be within list boundaries
  if (n >= list->used_size)
  {
    return 0;
  }

  image_free(list->images[n]);

  if (!cu_array_free(list->cu_arrays[n])) {
    fprintf(stderr, "Could not free cu_array!\n");
    assert(0); //Stop here
    return 0;
  }

  // The last item is easy to remove
  if (n == list->used_size - 1) {
    list->images[n] = NULL;
    list->cu_arrays[n] = NULL;
    list->pocs[n] = 0;
    list->used_size--;
  } else {
    int i = n;
    // Shift all following pics one backward in the list
    for (i = n; i < list->used_size - 1; ++i) {
      list->images[i] = list->images[i + 1];
      list->cu_arrays[i] = list->cu_arrays[i + 1];
      list->pocs[i] = list->pocs[i + 1];
    }
    list->images[list->used_size - 1] = NULL;
    list->cu_arrays[list->used_size - 1] = NULL;
    list->pocs[list->used_size - 1] = 0;
    list->used_size--;
  }

  return 1;
}

int image_list_copy_contents(image_list_t *target, image_list_t *source) {
  int i;
  while (target->used_size > 0) {
    image_list_rem(target, 0);
  }
  
  for (i = source->used_size - 1; i >= 0; --i) {
    image_list_add(target, source->images[i], source->cu_arrays[i], source->pocs[i]);
  }
  return 1;
}
