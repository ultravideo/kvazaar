#ifndef SEARCH_H_
#define SEARCH_H_
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
 * \brief Searching of parameters for intra and inter frames.
 */

#include "global.h"

#include "encoder.h"


void search_slice_data(encoder_control *encoder);
void search_tree(encoder_control *encoder, int x, int y, uint8_t depth);
uint32_t search_best_mode(encoder_control* encoder,uint16_t x_cu,uint16_t y_cu, uint8_t depth);

#endif
