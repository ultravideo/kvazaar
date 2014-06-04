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

#include "nal.h"
#include "strategyselector.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "bitstream.h"
#include "cabac.h"
#include "encoder.h"

/**
 * \brief Write a Network Abstraction Layer (NAL) packet to the output.
 */
void nal_write(bitstream * const bitstream, const uint8_t nal_type,
               const uint8_t temporal_id, const int long_start_code)
{
  uint8_t byte;

  // Some useful constants
  const uint8_t start_code_prefix_one_3bytes = 0x01;
  const uint8_t zero = 0x00;

  // zero_byte (0x00) shall be present in the byte stream NALU of VPS, SPS
  // and PPS, or the first NALU of an access unit
  if(long_start_code)
    bitstream_writebyte(bitstream, zero);

  // start_code_prefix_one_3bytes
  bitstream_writebyte(bitstream, zero);
  bitstream_writebyte(bitstream, zero);
  bitstream_writebyte(bitstream, start_code_prefix_one_3bytes);

  // Handle header bits with full bytes instead of using bitstream
  // forbidden_zero_flag(1) + nal_unit_type(6) + 1bit of nuh_layer_id
  byte = nal_type << 1;
  bitstream_writebyte(bitstream, byte);

  // 5bits of nuh_layer_id + nuh_temporal_id_plus1(3)
  byte = (temporal_id + 1) & 7;
  bitstream_writebyte(bitstream, byte);
}

/*!
 \brief Calculate checksums for all colors of the picture.
 \param pic The picture that checksum is calculated for.
 \param checksum_out Result of the calculation.
 \returns Void
*/
void picture_checksum(const picture* pic, unsigned char checksum_out[][SEI_HASH_MAX_LENGTH])
{
  array_checksum(pic->y_recdata, pic->height, pic->width, pic->width, checksum_out[0]);

  /* The number of chroma pixels is half that of luma. */
  array_checksum(pic->u_recdata, pic->height >> 1, pic->width >> 1, pic->width >> 1, checksum_out[1]);
  array_checksum(pic->v_recdata, pic->height >> 1, pic->width >> 1, pic->width >> 1, checksum_out[2]);
}
