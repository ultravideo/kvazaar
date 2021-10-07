/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "nal.h"

#include "bitstream.h"
#include "strategies/strategies-nal.h"


/**
 * \brief Write a Network Abstraction Layer (NAL) packet to the output.
 */
void kvz_nal_write(bitstream_t * const bitstream, const uint8_t nal_type,
               const uint8_t temporal_id, const int long_start_code)
{
  uint8_t byte;

  // Some useful constants
  const uint8_t start_code_prefix_one_3bytes = 0x01;
  const uint8_t zero = 0x00;

  // zero_byte (0x00) shall be present in the byte stream NALU of VPS, SPS
  // and PPS, or the first NALU of an access unit
  if(long_start_code)
    kvz_bitstream_writebyte(bitstream, zero);

  // start_code_prefix_one_3bytes
  kvz_bitstream_writebyte(bitstream, zero);
  kvz_bitstream_writebyte(bitstream, zero);
  kvz_bitstream_writebyte(bitstream, start_code_prefix_one_3bytes);

  // Handle header bits with full bytes instead of using bitstream
  // forbidden_zero_flag(1) + nal_unit_type(6) + 1bit of nuh_layer_id
  byte = nal_type << 1;
  kvz_bitstream_writebyte(bitstream, byte);

  // 5bits of nuh_layer_id + nuh_temporal_id_plus1(3)
  byte = (temporal_id + 1) & 7;
  kvz_bitstream_writebyte(bitstream, byte);
}

/*!
 \brief Calculate checksums for all colors of the picture.
 \param im The image that checksum is calculated for.
 \param checksum_out Result of the calculation.
 \returns Void
*/
void kvz_image_checksum(const kvz_picture *im, unsigned char checksum_out[][SEI_HASH_MAX_LENGTH], const uint8_t bitdepth)
{
  kvz_array_checksum(im->y, im->height, im->width, im->width, checksum_out[0], bitdepth);

  /* The number of chroma pixels is half that of luma. */
  if (im->chroma_format != KVZ_CSP_400) {
    kvz_array_checksum(im->u, im->height >> 1, im->width >> 1, im->width >> 1, checksum_out[1], bitdepth);
    kvz_array_checksum(im->v, im->height >> 1, im->width >> 1, im->width >> 1, checksum_out[2], bitdepth);
  }
}

/*!
\brief Calculate md5 for all colors of the picture.
\param im The image that md5 is calculated for.
\param checksum_out Result of the calculation.
\returns Void
*/
void kvz_image_md5(const kvz_picture *im, unsigned char checksum_out[][SEI_HASH_MAX_LENGTH], const uint8_t bitdepth)
{
  kvz_array_md5(im->y, im->height, im->width, im->width, checksum_out[0], bitdepth);

  /* The number of chroma pixels is half that of luma. */
  if (im->chroma_format != KVZ_CSP_400) {
    kvz_array_md5(im->u, im->height >> 1, im->width >> 1, im->width >> 1, checksum_out[1], bitdepth);
    kvz_array_md5(im->v, im->height >> 1, im->width >> 1, im->width >> 1, checksum_out[2], bitdepth);
  }
}
