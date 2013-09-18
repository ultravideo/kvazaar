/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing 2013.
 */

/*! \file nal.h
    \brief NAL
    \author Marko Viitanen
    \date 2013-02
    
    NAL function headers
*/
#ifndef __NAL_H
#define __NAL_H

#include "global.h"

#include <stdio.h>

#include "picture.h"


/*!
 * \brief NAL unit type codes
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

void nal_write(FILE* output, uint8_t* buffer, uint32_t buffer_len, uint8_t nal_ref, uint8_t nal_type, uint8_t temporal_id);
void picture_checksum(const picture* pic, unsigned char checksum_out[][16]);

#endif
