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

/*
* \file
*/

#include "strategies/avx2/resample-avx2.h"

#include "strategyselector.h"

#include "scaler/scaler-avx2.h"


int kvz_strategy_register_resample_avx2(void * opaque)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "resample_block_step", "avx2_2", 41, kvz_alt1_block_step_resample_func_avx2);
  success &= kvz_strategyselector_register(opaque, "resample_block_step", "avx2", 40, kvz_default_block_step_resample_func_avx2);
  
  success &= kvz_strategyselector_register(opaque, "resample", "avx2", 0, kvz_default_resample_func_avx2);

  return success;
}
