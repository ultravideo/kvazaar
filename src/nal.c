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
void nal_write(FILE *output, uint8_t nal_ref, uint8_t nal_type,
               uint8_t temporal_id,int long_start_code)
{
  uint8_t byte;

  // Some useful constants
  const uint8_t start_code_prefix_one_3bytes = 0x01;
  const uint8_t zero = 0x00;

  // zero_byte (0x00) shall be present in the byte stream NALU of VPS, SPS
  // and PPS, or the first NALU of an access unit
  if(long_start_code)
    fwrite(&zero, 1, 1, output);

  // start_code_prefix_one_3bytes
  fwrite(&zero, 1, 1, output);
  fwrite(&zero, 1, 1, output);
  fwrite(&start_code_prefix_one_3bytes, 1, 1, output);

  // Handle header bits with full bytes instead of using bitstream
  // forbidden_zero_flag(1) + nal_unit_type(6) + 1bit of nuh_layer_id
  byte = nal_type << 1;
  fwrite(&byte, 1, 1, output);

  // 5bits of nuh_layer_id + nuh_temporal_id_plus1(3)
  byte = (temporal_id + 1) & 7;
  fwrite(&byte, 1, 1, output);
}


/**
 * \brief Calculate checksum for one color of the picture.
 * \param data Beginning of the pixel data for the picture.
 * \param height Height of the picture.
 * \param width Width of the picture.
 * \param stride Width of one row in the pixel array.
 */
static void array_checksum(const pixel* data,
                           const int height, const int width,
                           const int stride,
                           unsigned char checksum_out[SEI_HASH_MAX_LENGTH])
{
  uint8_t mask;
  uint32_t checksum = 0;
  int y, x;

  assert(SEI_HASH_MAX_LENGTH >= 4);

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      mask = (uint8_t)((x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8));
      checksum += (data[(y * stride) + x] & 0xff) ^ mask;
    }
  }

  // Unpack uint into byte-array.
  checksum_out[0] = (checksum >> 24) & 0xff;
  checksum_out[1] = (checksum >> 16) & 0xff;
  checksum_out[2] = (checksum >> 8) & 0xff;
  checksum_out[3] = (checksum) & 0xff;
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
