#ifndef PICTURE_LIST_H_
#define PICTURE_LIST_H_
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

/**
 * \ingroup DataStructures
 * \file
 *  Container for a list of reference pictures.
 */

#include "cu.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


/**
 * \brief Struct which contains array of picture structs
 */
typedef struct
{
  struct kvz_picture* *images;          //!< \brief Pointer to array of picture pointers.
  cu_array_t* *cu_arrays;
  int32_t *pocs;
  uint8_t (*ref_LXs)[2][16]; //!< L0 and L1 reference index list for each image
  uint32_t size;       //!< \brief Array size.
  uint32_t used_size;


} image_list_t;

image_list_t * kvz_image_list_alloc(int size);
int kvz_image_list_resize(image_list_t *list, unsigned size);
int kvz_image_list_destroy(image_list_t *list);
int kvz_image_list_add(image_list_t *list, kvz_picture *im, cu_array_t* cua, int32_t poc, uint8_t ref_LX[2][16]);
int kvz_image_list_rem(image_list_t *list, unsigned n);

int kvz_image_list_copy_contents(image_list_t *target, image_list_t *source);

enum { REF_PIC_LIST_0 = 0, REF_PIC_LIST_1 = 1, REF_PIC_LIST_X = 100 };

#endif //PICTURE_LIST_H_
