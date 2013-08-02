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