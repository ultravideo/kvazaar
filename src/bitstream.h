/*! \file bitstream.h
    \brief Bitstream related functions
    \author Marko Viitanen
    \date 2012-05
    
    This file has all bitstream headers
*/
#ifndef _BITSTREAM_H
#define _BITSTREAM_H
 
 
typedef struct
{
    uint32_t data[32];
    uint8_t  cur_byte;
    uint8_t  cur_bit; 
 
} bitstream;
 
void bitstream_init(bitstream* stream);
 
void bitstream_put(bitstream* stream, uint32_t* data, uint8_t bits);
 
void bitstream_align(bitstream* stream);
 
void bitstream_flush(bitstream* stream);
 
#endif