#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
 
#include "global.h"
#include "bitstream.h"
 
 
void printf_bitstream(char *msg, ...)
{
 
va_list fmtargs;
char buffer[1024];
 
  va_start(fmtargs,msg);
  vsnprintf(buffer,sizeof(buffer)-1,msg,fmtargs);
  va_end(fmtargs);
  printf("%s",buffer);
 
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
 * Put bits to bitstream
 * Input:
 *          stream = pointer bitstream to put the data
 *          data   = pointer to actual data
 *          bits   = number of bits to write      
 */
 
void bitstream_put(bitstream* stream, uint32_t data, uint8_t bits)
{
    uint8_t i=0;
    uint32_t bitsleft=32-stream->cur_bit;
    printf_bitstream("put: ");
    for(i=0;i<bits;i++)
    {
        printf("%i",(data&(1<<(bits-i-1)))?1:0);
    }
    printf_bitstream("\n");
    //printf_bitstream(" count: %i\n",bits);
 
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
        stream->cur_byte++;
        bitsleft=32;
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
 *  Align the bitstream
 */
void bitstream_align(bitstream* stream)
{  
    if(stream->cur_byte==32)
    {
        //Stream flushed, zero out the values
        bitstream_init(stream);
    }
    else
    {
        stream->cur_byte++;
    }
}
 
void bitstream_flush(bitstream* stream)
{
   
    /*
     *  SAVE DATA TO OUTPUT
     */
 
    //Stream flushed, zero out the values
    bitstream_init(stream);
}
 