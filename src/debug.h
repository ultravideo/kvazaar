#ifndef DEBUG_H_
#define DEBUG_H_
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
 * \brief Tools for visualizing and debugging the encoder.
 */

#include "global.h"

#include <stdio.h>

#include "encoder.h"


FILE * open_cu_file(char *filename);
void close_cu_file(FILE *fp);
unsigned render_cu_file(encoder_control *encoder, picture *pic, unsigned depth, uint16_t x_cu, uint16_t y_cu, FILE *fp);

#endif
