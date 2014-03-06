/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 */

#include "bitstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
//for hton
#ifdef _WIN32
#include <Winsock2.h>
#else
#include <arpa/inet.h>
#endif


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

bit_table *g_exp_table;

//From wikipedia
//http://en.wikipedia.org/wiki/Binary_logarithm#Algorithm
int floor_log2(unsigned int n) {
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return ((n == 0) ? (-1) : pos);
}

/**
 * \brief Initialize the Exp Golomb code table with desired number of values
 * \param len table length to init
 * \return 1 on success, 0 on failure
 *
 * Allocates g_exp_table with len*sizeof(bit_table) and fills it with exponential golomb codes
 */
int init_exp_golomb(uint32_t len)
{
  uint32_t code_num;
  uint8_t M;
  uint32_t info;
  g_exp_table = (bit_table*)malloc(len*sizeof(bit_table));
  if(!g_exp_table)
    return 0;

  for (code_num = 0; code_num < len; code_num++) {
    M = (uint8_t)floor_log2(code_num + 1);
    info = code_num + 1 - (uint32_t)pow(2, M);
    g_exp_table[code_num].len = M * 2 + 1;
    g_exp_table[code_num].value = (1<<M) | info;
  }

  return 1;
}

/**
 * \brief Create and initialize a new bitstream
 */
bitstream *create_bitstream(int32_t width)
{
  bitstream *stream = malloc(sizeof(bitstream));
  if (!stream) {
    fprintf(stderr, "Failed to allocate the bitstream object!\n");
    return stream;
  }

  // Initialize the bitstream
  bitstream_reinit(stream);

  // Initialize buffer-related values
  stream->output     = NULL;
  stream->buffer     = NULL;
  stream->buffer_pos = 0;
  stream->bufferlen  = 0;

  // Alloc 2kB*width for bitstream buffer (for one coded frame)
  bitstream_alloc(stream, 1024*2*width);
  if (!stream->buffer) {
    fprintf(stderr, "Failed to allocate the bitstream buffer!\n");
    goto creation_failure;
  }

  //Clear buffer just to be sure
  bitstream_clear_buffer(stream);

  // Return the created bitstream
  return stream;

creation_failure:
  // In case of failure, free whatever was allocated
  free(stream->buffer);
  free(stream);

  return NULL;
}

/**
 * \brief Reinitialize bitstream
 */
void bitstream_reinit(bitstream *stream)
{
  stream->cur_byte = 0;
  stream->cur_bit = 0;
  memset(stream->data, 0, sizeof(uint32_t)*32);
}

/**
 *  \brief Allocate buffer
 *  \param stream pointer bitstream to put the data
 *  \param alloc size to allocate
 */
void bitstream_alloc(bitstream *stream, uint32_t alloc)
{
  stream->buffer = (uint8_t*)malloc(alloc);
  if (stream->buffer)
    stream->bufferlen = alloc;
}

/**
 *  \brief Free bitstream buffer
 */
void bitstream_free(bitstream *stream)
{
  free(stream->buffer);
  stream->bufferlen = 0;
}

/**
 *  \brief clear output buffer
 */
void bitstream_clear_buffer(bitstream *stream)
{
  memset(stream->buffer, 0, stream->bufferlen);
  stream->buffer_pos = 0;
}

/**
 * \brief Put bits to bitstream
 * \param stream pointer bitstream to put the data
 * \param data input data
 * \param bits number of bits to write from data to stream
 */
void bitstream_put(bitstream *stream, uint32_t data, uint8_t bits)
{
  uint8_t bitsleft = 32 - stream->cur_bit;
  #ifdef VERBOSE
  uint8_t i=0;
  printf_bitstream("put: ");
  for (i = 0; i < bits; i++) {
      printf("%i",(data&(1<<(bits-i-1)))?1:0);
  }
  printf_bitstream("\n");
  //printf_bitstream(" count: %i\n",bits);
  #endif

  //There's space for all the bits
  if (bits <= bitsleft) {
    stream->data[stream->cur_byte] |= (data<<((bitsleft-bits)));
    stream->cur_bit += bits;
    bits = 0;
  } else { //No space for everything, store the bits we can and continue later
    stream->data[stream->cur_byte] |= (data>>(bits-bitsleft));
    stream->cur_bit = 32;
    bits -= bitsleft;
  }

  //Check if the buffer is full, and flush to output if it is
  if (stream->cur_bit == 32) {
    bitsleft = 32;
    stream->cur_byte++;
    stream->cur_bit = 0;
    if (stream->cur_byte == 32) {
      //Flush data out
      bitstream_flush(stream);
    }
  }

  //Write the last of the bits (if buffer was full and flushed before)
  if (bits != 0) {
    stream->data[stream->cur_byte] |= (data<<(bitsleft-bits));
    stream->cur_bit += bits;
  }
}

/**
 * \brief Align the bitstream with one-bit padding
 */
void bitstream_align(bitstream *stream)
{
  bitstream_put(stream, 1, 1);
  if ((stream->cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}

/**
 * \brief Align the bitstream with zero
 */
void bitstream_align_zero(bitstream *stream)
{
  if ((stream->cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}

/**
 * \brief Flush bitstream to output
 */
void bitstream_flush(bitstream *stream)
{
  int i;
  uint32_t correct_endian;
  //If output open, write to output
  if (stream->output) {
    if (stream->cur_byte) fwrite(&stream->data[0], stream->cur_byte * 4, 1, stream->output);
    if (stream->cur_bit>>3) fwrite(&stream->data[stream->cur_byte], stream->cur_bit>>3, 1, stream->output);
  } else { //No file open, write to buffer
    if (stream->cur_byte) {
      //Handle endianness issue
      for (i = 0; i < stream->cur_byte; i++) {
        //"network" is big-endian
        correct_endian = htonl(stream->data[i]);
        memcpy((uint8_t*)&stream->buffer[stream->buffer_pos], &correct_endian, 4);
        stream->buffer_pos += 4;
      }
    }

    if (stream->cur_bit>>3) {
      correct_endian = htonl(stream->data[stream->cur_byte]);
      memcpy((uint8_t*)&stream->buffer[stream->buffer_pos], &correct_endian, stream->cur_bit>>3);
      stream->buffer_pos += stream->cur_bit>>3;
    }
  }
  //Stream flushed, zero out the values
  bitstream_reinit(stream);
}
