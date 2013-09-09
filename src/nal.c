/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file nal.c
    \brief NAL
    \author Marko Viitanen
    \date 2013-06
    
    NAL functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "cabac.h"
#include "encoder.h"


#include "nal.h"

void nal_write(FILE* output, uint8_t* buffer, uint32_t buffer_len, uint8_t nal_ref, uint8_t nal_type, uint8_t temporal_id)
{
  uint8_t byte;
  uint32_t i;
  uint8_t zerocount=0;

  /* Some useful constants */
  const uint8_t emulation_prevention_three_byte = 0x03;
  const uint8_t start_code_prefix_one_3bytes = 0x01;
  const uint8_t zero = 0x00;

  /*start_code_prefix_one_3bytes */  
  
  /*
  if(temporal_id == 0)
  {
    fwrite(&zero, 1, 1, output);
  }
  */  
  fwrite(&zero, 1, 1, output);
  fwrite(&zero, 1, 1, output);
  fwrite(&start_code_prefix_one_3bytes, 1, 1, output);

  /* Handle header bits with full bytes instead of using bitstream */
  /* forbidden_zero_flag(1) + nal_unit_type(6) + 1bit of nuh_layer_id*/
  byte = nal_type<<1;
  fwrite(&byte, 1, 1, output);

  /* 5bits of nuh_layer_id + nuh_temporal_id_plus1(3) */
  byte = (temporal_id+1)&7;
  fwrite(&byte, 1, 1, output);

  /* Write out bytes and add emulation_prevention_three_byte when needed */
  for(i = 0; i < buffer_len; i++)
  {
    if(zerocount == 2 && buffer[i] < 4) /* Prevent 0x0000 + 00/01/02/03 */
    {
      /* Inserting 0x03 */
      fwrite(&emulation_prevention_three_byte, 1, 1, output);
      zerocount = 0;
    }
    if(buffer[i] == 0)
    {
      zerocount++;
    }
    else
    {
      zerocount = 0;
    }

    /* Write the actual data */
    fwrite(&buffer[i], 1, 1, output);
  }

  /* If last byte was 0, add emulation_prevention_three_byte */
  if(buffer[buffer_len-1] == 0)
  {
    fwrite(&emulation_prevention_three_byte, 1, 1, output);
  }

}


/*!
 \brief Calculate checksum for one color of the picture.
 \param data Beginning of the pixel data for the picture.
 \param height Height of the picture.
 \param width Width of the picture.
 \param stride Width of one row in the pixel array.
 \returns Void
*/
static void array_checksum(const uint8_t* data, const int height, const int width, const int stride, unsigned char checksum_out[])
{
	unsigned char mask;
	unsigned int checksum = 0;
	int y, x;
	for (y = 0; y < height; ++y) {
		for (x = 0; x < width; ++x) {
			mask = (x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8);
			checksum += (data[(y * stride) + x] & 0xff) ^ mask;
			checksum &= 0xffffffff;
		}
	}

  /* Unpack uint into byte-array.*/
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
void picture_checksum(const picture* pic, unsigned char checksum_out[][16])
{
	int stride = pic->width; /* ToDo: != width, if there is a luma margin. */
	array_checksum(pic->yRecData, pic->height, pic->width, pic->width, checksum_out[0]);

  /* The number of chroma pixels is half that of luma. */
	array_checksum(pic->uRecData, pic->height >> 1, pic->width >> 1, pic->width >> 1, checksum_out[1]);
	array_checksum(pic->vRecData, pic->height >> 1, pic->width >> 1, pic->width >> 1, checksum_out[2]);
}