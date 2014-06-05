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
#include "picturelist.h"
#include "strategyselector.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "sao.h"


/**
 * \brief Allocate memory for picture_list
 * \param size  initial array size
 * \return picture_list pointer, NULL on failure
 */
picture_list * picture_list_init(int size)
{
  picture_list *list = (picture_list *)malloc(sizeof(picture_list));
  list->size = size;
  if (size > 0) {
    list->pics = (picture**)malloc(sizeof(picture*) * size);
  }

  list->used_size = 0;

  return list;
}

/**
 * \brief Resize picture_list array
 * \param list  picture_list pointer
 * \param size  new array size
 * \return 1 on success, 0 on failure
 */
int picture_list_resize(picture_list *list, unsigned size)
{
  unsigned int i;
  picture** old_pics = NULL;

  // No need to do anything when resizing to same size
  if (size == list->size) {
    return 1;
  }

  // Save the old list
  if (list->used_size > 0) {
    old_pics = list->pics;
  }

  // allocate space for the new list
  list->pics = (picture**)malloc(sizeof(picture*)*size);

  // Copy everything from the old list to the new if needed.
  if (old_pics != NULL) {
    for (i = 0; i < list->used_size; ++i) {
      list->pics[i] = old_pics[i];
    }

    free(old_pics);
  }

  return 1;
}

/**
 * \brief Free memory allocated to the picture_list
 * \param list picture_list pointer
 * \return 1 on success, 0 on failure
 */
int picture_list_destroy(picture_list *list)
{
  unsigned int i;
  if (list->used_size > 0) {
    for (i = 0; i < list->used_size; ++i) {
      picture_free(list->pics[i]);
      list->pics[i] = NULL;
    }
  }

  if (list->size > 0) {
    free(list->pics);
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
int picture_list_add(picture_list *list, picture* pic)
{
  int i = 0;
  if (ATOMIC_INC(&(pic->refcount)) == 1) {
    fprintf(stderr, "Tried to add an unreferenced picture. This is a bug!\n");
    assert(0); //Stop for debugging
    return 0;
  }

  if (list->size == list->used_size) {
    if (!picture_list_resize(list, list->size*2)) return 0;
  }

  for (i = list->used_size; i > 0; i--) {
    list->pics[i] = list->pics[i - 1];
  }

  list->pics[0] = pic;
  list->used_size++;
  return 1;
}

/**
 * \brief Remove picture from picturelist
 * \param list list to use
 * \param n index to remove
 * \return 1 on success
 */
int picture_list_rem(picture_list * const list, const unsigned n)
{
  // Must be within list boundaries
  if (n >= list->used_size)
  {
    return 0;
  }

  if (!picture_free(list->pics[n])) {
    fprintf(stderr, "Could not free picture!\n");
    assert(0); //Stop here
    return 0;
  }

  // The last item is easy to remove
  if (n == list->used_size - 1) {
    list->pics[n] = NULL;
    list->used_size--;
  } else {
    int i = n;
    // Shift all following pics one backward in the list
    for (i = n; i < list->used_size - 1; ++i) {
      list->pics[i] = list->pics[i + 1];
    }
    list->pics[list->used_size - 1] = NULL;
    list->used_size--;
  }

  return 1;
}
