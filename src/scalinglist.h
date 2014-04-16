#ifndef SCALINGLIST_H_
#define SCALINGLIST_H_
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

#include "global.h"

typedef struct {
        int8_t   enable;
        int32_t  scaling_list_dc   [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
  const int32_t *scaling_list_coeff[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
  const int32_t *quant_coeff[4][6][6];
  const int32_t *de_quant_coeff  [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM][SCALING_LIST_REM_NUM];
  const double *error_scale[4][6][6];
} scaling_list;

extern const uint8_t g_scaling_list_num[4];
extern const uint16_t g_scaling_list_size[4];

const int32_t *scalinglist_get_default(const uint32_t size_id, const uint32_t list_id);

void scalinglist_init(scaling_list * const scaling_list);
void scalinglist_destroy(scaling_list * const scaling_list);

int  scalinglist_parse(scaling_list * const scaling_list, FILE *fp);
void scalinglist_process(scaling_list * const scaling_list);

//void scalinglist_set(scaling_list * const scaling_list, const int32_t * const coeff, uint32_t listId, uint32_t sizeId, uint32_t qp);
//void scalinglist_set_err_scale(scaling_list * const scaling_list, uint32_t list, uint32_t size, uint32_t qp);






#endif