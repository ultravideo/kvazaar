#ifndef PICTURE_LIST_H_
#define PICTURE_LIST_H_
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
 * \brief Coding Unit (CU) and picture data related functions.
 */

#include "picture.h"

/**
 * \brief Struct which contains array of picture structs
 */
typedef struct
{
  struct picture** pics;          //!< \brief Pointer to array of picture pointers.
  uint32_t size;       //!< \brief Array size.
  uint32_t used_size;
} picture_list;

picture_list * picture_list_init(int size);
int picture_list_resize(picture_list *list, unsigned size);
int picture_list_destroy(picture_list *list);
int picture_list_add(picture_list *list, picture *pic);
int picture_list_rem(picture_list *list, unsigned n);

#endif //PICTURE_LIST_H_