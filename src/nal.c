/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file nal.c
    \brief NAL
    \author Marko Viitanen
    \date 2012-06
    
    NAL functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "cabac.h"
#include "picture.h"
#include "nal.h"

void nal_write(FILE* output, uint8_t* buffer, uint32_t buffer_len, uint8_t nal_ref, uint8_t nal_type, uint8_t temporal_id)
{
  uint8_t byte;
  uint32_t i;
  uint8_t zerocount=0;
  uint8_t emulation_prevention_three_byte = 0x03;
  uint8_t start_code_prefix_one_3bytes = 0x01;
  uint8_t zero = 0x00;
  
  //start_code_prefix_one_3bytes
  fwrite(&zero, 1, 1, output);
  fwrite(&zero, 1, 1, output);
  fwrite(&zero, 1, 1, output);
  fwrite(&start_code_prefix_one_3bytes, 1, 1, output);

  //forbidden_zero_flag(1) + nal_ref_flag(1) + nal_unit_type(6)
  byte = nal_ref<<6 | nal_type;
  fwrite(&byte, 1, 1, output);

  //Temporal_id(3) + reserved_one_5bits(5)
  byte = temporal_id << 5 | 1;
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
      zerocount++;
    else
      zerocount = 0;

    fwrite(&buffer[i], 1, 1, output);
  }
  //If last byte was 0, add emulation_prevention_three_byte
  if(buffer[buffer_len-1] == 0)
    fwrite(&emulation_prevention_three_byte, 1, 1, output);

}