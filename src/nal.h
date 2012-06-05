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

#define NAL_SEQ_PARAMETER_SET 7
#define NAL_PIC_PARAMETER_SET 8

void nal_write(FILE* output, uint8_t* buffer, uint32_t buffer_len, uint8_t nal_ref, uint8_t nal_type, uint8_t temporal_id);