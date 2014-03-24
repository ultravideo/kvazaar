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

const uint32_t bit_set_mask[] =
{
0x00000001,0x00000002,0x00000004,0x00000008,
0x00000010,0x00000020,0x00000040,0x00000080,
0x00000100,0x00000200,0x00000400,0x00000800,
0x00001000,0x00002000,0x00004000,0x00008000,
0x00010000,0x00020000,0x00040000,0x00080000,
0x00100000,0x00200000,0x00400000,0x00800000,
0x01000000,0x02000000,0x04000000,0x08000000,
0x10000000,0x20000000,0x40000000,0x80000000
};


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
bitstream *create_bitstream()
{
  bitstream *stream = malloc(sizeof(bitstream));
  if (!stream) {
    fprintf(stderr, "Failed to allocate the bitstream object!\n");
    return stream;
  }

  // Initialize buffer-related values
  stream->data = 0;
  stream->cur_bit = 0;
  stream->zerocount = 0;

  // Return the created bitstream
  return stream;
}

/**
 * \brief Put bits to bitstream
 * \param stream pointer bitstream to put the data
 * \param data input data
 * \param bits number of bits to write from data to stream
 */
void bitstream_put(bitstream *stream, uint32_t data, uint8_t bits)
{
  const uint8_t emulation_prevention_three_byte = 0x03;
  while(bits--) {
    stream->data <<= 1;

    if (data & bit_set_mask[bits]) {
      stream->data |= 1;
    }
    stream->cur_bit++;

  // write byte to output
    if (stream->cur_bit==8) {
      if((stream->zerocount == 2) && (stream->data < 4)) {
        fwrite(&emulation_prevention_three_byte, 1, 1, stream->output);
        stream->zerocount = 0;
      }
      if(stream->data == 0) {
        stream->zerocount++;
      } else {
        stream->zerocount = 0;
      }
      fwrite(&stream->data, 1, 1, stream->output);
      stream->cur_bit = 0;
    }
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