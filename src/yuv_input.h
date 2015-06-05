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

#ifndef YUV_INPUT_H_
#define YUV_INPUT_H_

/*
 * \file
 * \brief Functions related to reading YUV input.
 */

#include "global.h"

int read_yuv_frame(FILE* file,
                   unsigned input_width, unsigned input_height,
                   unsigned array_width, unsigned array_height,
                   image_t *img_out);

#endif // YUV_INPUT_H_
