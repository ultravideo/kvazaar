#ifndef NAL_H_
#define NAL_H_
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
 * \brief Network Abstraction Layer (NAL) messages.
 */

#include "global.h"

#include <stdio.h>

#include "image.h"
#include "bitstream.h"


//////////////////////////////////////////////////////////////////////////
// TYPES

/**
 * \brief NAL unit type codes.
 *
 * These are the nal_unit_type codes from Table 7-1 ITU-T H.265 v1.0.
 * The type codes have been prefixed with "NAL_".
 */
enum {
  NAL_TRAIL_N = 0,
  NAL_TRAIL_R = 1,

  NAL_TSA_N = 2,
  NAL_TSA_R = 3,

  NAL_STSA_N = 4,
  NAL_STSA_R = 5,

  NAL_RADL_N = 6,
  NAL_RADL_R = 7,

  NAL_RASL_N = 8,
  NAL_RASL_R = 9,

  // Reserved RSV_VCL_ N/R 10-15

  NAL_BLA_W_LP = 16,
  NAL_BLA_W_RADL = 17,
  NAL_BLA_N_LP = 18,

  NAL_IDR_W_RADL = 19,
  NAL_IDR_N_LP = 20,

  NAL_CRA_NUT = 21,

  // Reserved RSV_IRAP_VCL 22-23
  NAL_RSV_IRAP_VCL23 = 23,

  // Reserved RSV_VCL 24-31

  NAL_VPS_NUT = 32,
  NAL_SPS_NUT = 33,
  NAL_PPS_NUT = 34,

  AUD_NUT = 35,
  EOS_NUT = 36,
  EOB_NUT = 37,
  FD_NUT = 38,

  PREFIX_SEI_NUT = 39,
  NAL_SUFFIT_SEI_NUT = 40,

  // Reserved RSV_NVCL 41-47
  // Unspecified UNSPEC 48-63
};

#define SEI_HASH_MAX_LENGTH 4

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
void kvz_nal_write(bitstream_t * const bitstream, const uint8_t nal_type,
               const uint8_t temporal_id, const int long_start_code);
void kvz_image_checksum(const kvz_picture *im,
                      unsigned char checksum_out[][SEI_HASH_MAX_LENGTH], const uint8_t bitdepth);



#endif
