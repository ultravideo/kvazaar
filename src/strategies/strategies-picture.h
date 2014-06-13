#ifndef STRATEGIES_PICTURE_H_
#define STRATEGIES_PICTURE_H_
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

#include "../image.h"

//Function pointer to reg_sad
typedef unsigned(*reg_sad_func)(const pixel *const data1, const pixel *const data2,
                                const int width, const int height,
                                const unsigned stride1, const unsigned stride2);
extern reg_sad_func reg_sad;


int strategy_register_picture(void* opaque);


#define STRATEGIES_PICTURE_EXPORTS {"reg_sad", (void**) &reg_sad}

#endif //STRATEGIES_PICTURE_H_
