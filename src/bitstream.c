/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file bitstream.c
    \brief Bitstream related functions
    \author Marko Viitanen
    \date 2013-02
    
    This file has all bitstream functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
/* for hton */
#ifdef WIN32
#include <Winsock2.h>
#else
#include <net/hton.h>
#endif
 
#include "global.h"
#include "bitstream.h"

//#define VERBOSE

#ifdef VERBOSE
void printf_bitstream(char *msg, ...)
{ 
  va_list fmtargs;
  char buffer[1024]; 
  va_start(fmtargs,msg);
  vsnprintf(buffer,sizeof(buffer)-1,msg,fmtargs);
  va_end(fmtargs);
  printf("%s",buffer); 
} 
#endif

bitTable *g_exp_table;

//From wikipedia
//http://en.wikipedia.org/wiki/Binary_logarithm#Algorithm
int floorLog2(unsigned int n) {
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return ((n == 0) ? (-1) : pos);  
}

//Initialize the Exp Golomb code table with desired number of values
void init_exp_golomb(uint32_t len)
{
    uint32_t code_num;
    uint32_t M;
    uint32_t info;
    g_exp_table=(bitTable*)malloc(len*sizeof(bitTable));    
    
    for(code_num=0;code_num<len;code_num++)
    {
        M=(uint32_t)floorLog2(code_num+1);
        info=code_num+1-(uint32_t)pow(2,M);        
        g_exp_table[code_num].len=M*2+1;
        g_exp_table[code_num].value=(1<<M)|info;
        //printf_cavlc("Len: %i %x\n", M*2+1, (1<<M)|info);
    }
}

/*
 * Clear bitstream
 */
void bitstream_init(bitstream* stream)
{
    stream->cur_byte=0;
    stream->cur_bit=0;
    memset(stream->data, 0, sizeof(uint32_t)*32); 
}

/*
 *  Allocate buffer
 */
void bitstream_alloc(bitstream* stream, uint32_t alloc)
{
  stream->buffer = (uint8_t*)malloc(alloc);
  stream->bufferlen = alloc;
  //Clear just to be sure
  bitstream_clear_buffer(stream);
}

void bitstream_clear_buffer(bitstream* stream)
{
  memset(stream->buffer,0,stream->bufferlen);
  stream->buffer_pos = 0;
}
 
 
/*
 * Put bits to bitstream
 * Input:
 *          stream = pointer bitstream to put the data
 *          data   = pointer to actual data
 *          bits   = number of bits to write      
 */
 
void bitstream_put(bitstream* stream, uint32_t data, uint8_t bits)
{
  uint32_t bitsleft=32-stream->cur_bit;
  #ifdef VERBOSE
  uint8_t i=0;
  printf_bitstream("put: ");
  for(i=0;i<bits;i++)
  {
      printf("%i",(data&(1<<(bits-i-1)))?1:0);
  }
  printf_bitstream("\n");
  //printf_bitstream(" count: %i\n",bits);
  #endif
 
  //Theres space for all the bits
  if(bits<=bitsleft)
  {
    stream->data[stream->cur_byte] |= (data<<((bitsleft-bits)));
    stream->cur_bit+=bits;
    bits=0;
  }
  //No space for everything, store the bits we can and continue later
  else
  {
    stream->data[stream->cur_byte] |= (data>>(bits-bitsleft));
    stream->cur_bit=32;
    bits-=bitsleft;
  }
 
  //Check if the buffer is full
  if(stream->cur_bit==32)
  {
    bitsleft=32;
    stream->cur_byte++;
    stream->cur_bit = 0;
    if(stream->cur_byte==32)
    {
        //Flush data out
        bitstream_flush(stream);
    }
  }
 
  //..still some writing to do
  if(bits!=0)
  {
    stream->data[stream->cur_byte] |= (data<<(bitsleft-bits));
    stream->cur_bit+=bits;
  } 
}
 
/*
 *  \brief Align the bitstream
 */
void bitstream_align(bitstream* stream)
{  
  if((stream->cur_bit&7) != 0)
  {
    bitstream_put(stream,0, 8-(stream->cur_bit&7));
  }
}
 
void bitstream_flush(bitstream* stream)
{
   /*
    *  SAVE DATA TO OUTPUT
    */
  int i;
  uint32_t correct_endian;
  if(stream->output)
  {
    if(stream->cur_byte)
    {
      fwrite(&stream->data[0], stream->cur_byte*4, 1, stream->output);
    }
   
    if(stream->cur_bit>>3)
    {
      fwrite(&stream->data[stream->cur_byte], stream->cur_bit>>3, 1, stream->output);
    }
  }
  /* No file open, write to buffer */   
  else
  {
    if(stream->cur_byte)
    {
      /* Handle endianness issue */
      for(i = 0; i < stream->cur_byte; i++)
      {
        /* "network" is big-endian */
        correct_endian = htonl(stream->data[i]);
        memcpy((uint8_t*)&stream->buffer[stream->buffer_pos],&correct_endian,4);
        stream->buffer_pos += 4;
      }
    }
   
    if(stream->cur_bit>>3)
    {
      correct_endian = htonl(stream->data[stream->cur_byte]);
      memcpy((uint8_t*)&stream->buffer[stream->buffer_pos],&correct_endian,stream->cur_bit>>3);
      stream->buffer_pos += stream->cur_bit>>3;
    }
  }
  //Stream flushed, zero out the values
  bitstream_init(stream);
}

