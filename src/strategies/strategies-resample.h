#ifndef STRATEGIES_RESAMPLE_H_
#define STRATEGIES_RESAMPLE_H_
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

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "scaler/scaler.h"

extern resample_block_step_func * kvz_resample_block_step;
extern resample_func * kvz_resample;

int kvz_strategy_register_resample(void *opaque);

#define STRATEGIES_RESAMPLE_EXPORTS \
  {"resample_block_step", (void**) &kvz_resample_block_step}, \
  {"resample", (void**) &kvz_resample}, \

#endif