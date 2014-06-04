#ifndef ENCODER_STATE_BITSTREAM_H_
#define ENCODER_STATE_BITSTREAM_H_
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


void encoder_state_write_bitstream_slice_header(struct encoder_state * const encoder_state);
void encoder_state_write_bitstream(struct encoder_state * const main_state);
void encoder_state_write_bitstream_leaf(struct encoder_state * const encoder_state);
void encoder_state_worker_write_bitstream_leaf(void * opaque);


#endif // ENCODER_STATE_BITSTREAM_H_
