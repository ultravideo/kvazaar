#ifndef KVAZAAR_INTERNAL_H_
#define KVAZAAR_INTERNAL_H_
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

#include "kvazaar.h"

// Forward declarations.
struct encoder_state_t;
struct encoder_control_t;

struct kvz_encoder {
  struct encoder_control_t* control;
  struct encoder_state_t* states;
  unsigned num_encoder_states;
  unsigned cur_state_num;
  unsigned frames_started;
  unsigned frames_done;

  size_t bitstream_length;
};

#endif // KVAZAAR_INTERNAL_H_
