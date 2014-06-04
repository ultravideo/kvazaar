#ifndef ENCODER_STATE_CTORS_DTORS_H_
#define ENCODER_STATE_CTORS_DTORS_H_
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
 */

#include "global.h"


// Forward declare because including the header would lead  to a cyclic
// dependency.
struct encoder_state;


int encoder_state_init(struct encoder_state * child_state, struct encoder_state * parent_state);
void encoder_state_finalize(struct encoder_state *encoder_state);
void encoder_state_init_lambda(struct encoder_state *encoder_state);


#endif // ENCODER_STATE_CTORS_DTORS_H_
