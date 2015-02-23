#ifndef STRATEGIES_NAL_H_
#define STRATEGIES_NAL_H_
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

#include "../nal.h"

//Function pointer to array_checksum
/**
 * \brief Calculate checksum for one color of the picture.
 * \param data Beginning of the pixel data for the picture.
 * \param height Height of the picture.
 * \param width Width of the picture.
 * \param stride Width of one row in the pixel array.
 */
typedef void (*array_checksum_func)(const pixel* data,
                                    const int height, const int width,
                                    const int stride,
                                    unsigned char checksum_out[SEI_HASH_MAX_LENGTH]);
extern array_checksum_func array_checksum;


int strategy_register_nal(void* opaque);


#define STRATEGIES_NAL_EXPORTS \
  {"array_checksum", (void**) &array_checksum},

#endif //STRATEGIES_NAL_H_
