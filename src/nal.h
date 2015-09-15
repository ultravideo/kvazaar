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
 */
enum kvz_nal_unit_type {

  // Trailing pictures

  KVZ_NAL_TRAIL_N = 0,
  KVZ_NAL_TRAIL_R = 1,

  KVZ_NAL_TSA_N = 2,
  KVZ_NAL_TSA_R = 3,

  KVZ_NAL_STSA_N = 4,
  KVZ_NAL_STSA_R = 5,

  // Leading pictures

  KVZ_NAL_RADL_N = 6,
  KVZ_NAL_RADL_R = 7,

  KVZ_NAL_RASL_N = 8,
  KVZ_NAL_RASL_R = 9,

  // Reserved non-IRAP RSV_VCL_N/R 10-15

  // Intra random access point pictures

  KVZ_NAL_BLA_W_LP   = 16,
  KVZ_NAL_BLA_W_RADL = 17,
  KVZ_NAL_BLA_N_LP   = 18,

  KVZ_NAL_IDR_W_RADL = 19,
  KVZ_NAL_IDR_N_LP   = 20,

  KVZ_NAL_CRA_NUT    = 21,

  // Reserved IRAP

  KVZ_NAL_RSV_IRAP_VCL22 = 22,
  KVZ_NAL_RSV_IRAP_VCL23 = 23,

  // Reserved non-IRAP RSV_VCL 24-32

  // non-VCL

  KVZ_NAL_VPS_NUT = 32,
  KVZ_NAL_SPS_NUT = 33,
  KVZ_NAL_PPS_NUT = 34,

  KVZ_NAL_AUD_NUT = 35,
  KVZ_NAL_EOS_NUT = 36,
  KVZ_NAL_EOB_NUT = 37,
  KVZ_NAL_FD_NUT  = 38,

  KVZ_NAL_PREFIX_SEI_NUT = 39,
  KVZ_NAL_SUFFIX_SEI_NUT = 40,

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
