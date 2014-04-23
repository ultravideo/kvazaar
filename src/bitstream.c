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
#include <stdlib.h>
#include <assert.h>
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

const bit_table *g_exp_table;

//From wikipedia
//http://en.wikipedia.org/wiki/Binary_logarithm#Algorithm
static int floor_log2(unsigned int n) {
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
int init_exp_golomb(const uint32_t len)
{
  uint32_t code_num;
  uint8_t M;
  uint32_t info;
  bit_table *exp_table;
  exp_table = (bit_table*)malloc(len*sizeof(bit_table));
  if(!exp_table)
    return 0;

  for (code_num = 0; code_num < len; code_num++) {
    M = (uint8_t)floor_log2(code_num + 1);
    info = code_num + 1 - (uint32_t)pow(2, M);
    exp_table[code_num].len = M * 2 + 1;
    exp_table[code_num].value = (1<<M) | info;
  }
  
  g_exp_table = exp_table;

  return 1;
}

/**
 * \brief Free Exp Golomb tables
 */
void free_exp_golomb()
{
  FREE_POINTER(g_exp_table);
}

/**
 * \brief Initialize a new bitstream
 */
int bitstream_init(bitstream * const stream, const bitstream_type type) {
  switch (type) {
    case BITSTREAM_TYPE_MEMORY:
      stream->mem.allocated_length = 0;
      stream->mem.output_data = NULL;
      stream->mem.output_length = 0;
      break;
      
    case BITSTREAM_TYPE_FILE:
      stream->file.output = NULL;
      break;
      
    default:
      fprintf(stderr, "Unknown type for bitstream!\n");
      return 0;
  }
  
  stream->base.cur_bit = 0;
  stream->base.data = 0;
  stream->base.zerocount = 0;
  stream->base.type = type;
  
  return 1;
}

/**
 * \brief Finalize bitstream internal structures
 */

int bitstream_finalize(bitstream * const stream) {
  switch (stream->base.type) {
    case BITSTREAM_TYPE_MEMORY:
      FREE_POINTER(stream->mem.output_data);
      stream->mem.allocated_length = 0;
      stream->mem.output_length = 0;
      break;
      
    case BITSTREAM_TYPE_FILE:
      //FIXME: if we fix create_bitstream, we would maybe have to do something here
      stream->file.output = NULL;
      break;
      
    default:
      fprintf(stderr, "Unknown type for bitstream!\n");
      return 0;
  }
  
  return 1;
}


/**
 * \brief Write a byte to bitstream
 * \param stream pointer bitstream to put the data
 * \param byte byte to write
 * \return 1 on success, 0 on failure
 */
int bitstream_writebyte(bitstream * const stream, const uint8_t byte) {
  switch (stream->base.type) {
    case BITSTREAM_TYPE_FILE:
      if (fwrite(&byte, 1, 1, stream->file.output) != 1) {
        fprintf(stderr, "Could not write byte to bitstream_file object.");
        return 0;
      }
      break;
      
    case BITSTREAM_TYPE_MEMORY:
      if (stream->mem.allocated_length==stream->mem.output_length) {
        //Need to reallocate
        uint32_t new_size = stream->mem.allocated_length + BITSTREAM_MEMORY_CHUNK_SIZE;
        uint8_t* new_data = realloc(stream->mem.output_data, new_size);
        if (!new_data) {
          fprintf(stderr, "Failed to allocate memory for bitstream_mem object");
          return 0;
        }
        stream->mem.output_data = new_data;
        stream->mem.allocated_length = new_size;
      }
      //Write byte
      stream->mem.output_data[stream->mem.output_length++] = byte;
      break;
      
    default:
      fprintf(stderr, "Unknown stream type!\n");
      assert(0);
      return 0;
  }
  return 1;
}

/**
 * \brief Put bits to bitstream
 * \param stream pointer bitstream to put the data
 * \param data input data
 * \param bits number of bits to write from data to stream
 */
void bitstream_put(bitstream * const stream, const uint32_t data, uint8_t bits)
{
  const uint8_t emulation_prevention_three_byte = 0x03;
  while(bits--) {
    stream->base.data <<= 1;

    if (data & bit_set_mask[bits]) {
      stream->base.data |= 1;
    }
    stream->base.cur_bit++;

  // write byte to output
    if (stream->base.cur_bit==8) {
      if((stream->base.zerocount == 2) && (stream->base.data < 4)) {
        bitstream_writebyte(stream, emulation_prevention_three_byte);
        stream->base.zerocount = 0;
      }
      if(stream->base.data == 0) {
        stream->base.zerocount++;
      } else {
        stream->base.zerocount = 0;
      }
      bitstream_writebyte(stream, stream->base.data);
      stream->base.cur_bit = 0;
    }
  }
}

/**
 * \brief Align the bitstream with one-bit padding
 */
void bitstream_align(bitstream * const stream)
{
  bitstream_put(stream, 1, 1);
  if ((stream->base.cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->base.cur_bit & 7));
  }
}

/**
 * \brief Align the bitstream with zero
 */
void bitstream_align_zero(bitstream * const stream)
{
  if ((stream->base.cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->base.cur_bit & 7));
  }
}
