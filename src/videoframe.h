#ifndef VIDEOFRAME_H_
#define VIDEOFRAME_H_
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
 * \brief Video frame stuff
 */

#include "global.h"
#include "cu.h"
#include "image.h"

struct sao_info_struct;

/**
 * \brief Struct which contains all picture data
 */
typedef struct videoframe
{
  image* source;         //!< \brief Source image.
  image* rec;            //!< \brief Reconstructed image.
  
  coefficient* coeff_y;   //!< \brief coefficient pointer Y
  coefficient* coeff_u;   //!< \brief coefficient pointer U
  coefficient* coeff_v;   //!< \brief coefficient pointer V

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.

  cu_info* cu_array;     //!< \brief Info for each CU at each depth.
  struct sao_info_struct *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_struct *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  int32_t poc;           //!< \brief Picture order count
} videoframe;


videoframe *videoframe_alloc(int32_t width, int32_t height, int32_t poc);
int videoframe_free(videoframe * const frame);

void videoframe_set_poc(videoframe * frame, int32_t poc);

const cu_info* videoframe_get_cu_const(const videoframe * const frame, unsigned int x_in_scu, unsigned int y_in_scu);
cu_info* videoframe_get_cu(videoframe * const frame, const unsigned int x_in_scu, const unsigned int y_in_scu);
void videoframe_compute_psnr(const videoframe * const frame, double psnr[NUM_COLORS]);


#endif
