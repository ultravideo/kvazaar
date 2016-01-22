#ifndef VIDEOFRAME_H_
#define VIDEOFRAME_H_
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

/**
 * \ingroup DataStructures
 * \file
 * \brief Container for the frame currently being encoded.
 */

#include "global.h"

#include "cu.h"
#include "image.h"

struct sao_info_t;

/**
 * \brief Struct which contains all picture data
 */
typedef struct videoframe
{
  kvz_picture *source;         //!< \brief Source image.
  kvz_picture *rec;            //!< \brief Reconstructed image.

  coeff_t* coeff_y;   //!< \brief coefficient pointer Y
  coeff_t* coeff_u;   //!< \brief coefficient pointer U
  coeff_t* coeff_v;   //!< \brief coefficient pointer V

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.

  cu_array_t* cu_array;     //!< \brief Info for each CU at each depth.
  struct sao_info_t *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_t *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  int32_t poc;           //!< \brief Picture order count
} videoframe_t;


videoframe_t *kvz_videoframe_alloc(int32_t width, int32_t height, int32_t poc);
int kvz_videoframe_free(videoframe_t * const frame);

void kvz_videoframe_set_poc(videoframe_t * frame, int32_t poc);

const cu_info_t* kvz_videoframe_get_cu_const(const videoframe_t * const frame, unsigned int x_in_scu, unsigned int y_in_scu);
cu_info_t* kvz_videoframe_get_cu(videoframe_t * const frame, const unsigned int x_in_scu, const unsigned int y_in_scu);

#endif
