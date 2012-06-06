/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems 2012.
 */

/*! \file nal.h
    \brief NAL
    \author Marko Viitanen
    \date 2012-06
    
    NAL function headers
*/
#ifndef __NAL_H
#define __NAL_H

enum { NAL_IDR_SLICE = 0x45, NAL_SEQ_PARAMETER_SET = 7, NAL_PIC_PARAMETER_SET = 8 };

void nal_write(FILE* output, uint8_t* buffer, uint32_t buffer_len, uint8_t nal_ref, uint8_t nal_type, uint8_t temporal_id);

#endif
