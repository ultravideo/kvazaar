/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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

bit_table_t g_exp_table[EXP_GOLOMB_TABLE_SIZE];


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

static int floor_log2(unsigned int n)
{
  assert(n != 0);

  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return pos;
}

/**
 * \brief Initialize the Exp Golomb code table.
 *
 * Fills g_exp_table with exponential golomb codes.
 */
void init_exp_golomb()
{
  static int exp_table_initialized = 0;
  if (exp_table_initialized) return;

  uint32_t code_num;
  uint8_t M;
  uint32_t info;
  for (code_num = 0; code_num < EXP_GOLOMB_TABLE_SIZE; code_num++) {
    M = (uint8_t)floor_log2(code_num + 1);
    info = code_num + 1 - (uint32_t)pow(2, M);
    g_exp_table[code_num].len = M * 2 + 1;
    g_exp_table[code_num].value = (1<<M) | info;
  }

  exp_table_initialized = 1;
}

/**
 * \brief Initialize a new bitstream.
 */
void bitstream_init(bitstream_t *const stream)
{
  memset(stream, 0, sizeof(bitstream_t));
}

/**
 * \brief Take chunks from a bitstream.
 *
 * Move ownership of the chunks to the caller and clear the bitstream.
 *
 * The bitstream must be byte-aligned.
 */
kvz_data_chunk * bitstream_take_chunks(bitstream_t *const stream)
{
  assert(stream->cur_bit == 0);
  kvz_data_chunk *chunks = stream->first;
  stream->first = stream->last = NULL;
  stream->len = 0;
  return chunks;
}

/**
 * \brief Allocates a new bitstream chunk.
 *
 * \return Pointer to the new chunk, or NULL.
 */
kvz_data_chunk * bitstream_alloc_chunk()
{
    kvz_data_chunk *chunk = malloc(sizeof(kvz_data_chunk));
    if (chunk) {
      chunk->len = 0;
      chunk->next = NULL;
    }
    return chunk;
}

/**
 * \brief Free a list of chunks.
 */
void bitstream_free_chunks(kvz_data_chunk *chunk)
{
  while (chunk != NULL) {
    kvz_data_chunk *next = chunk->next;
    free(chunk);
    chunk = next;
  }
}

/**
 * \brief Free resources used by a bitstream.
 */
void bitstream_finalize(bitstream_t *const stream)
{
  bitstream_clear(stream);
}

/**
 * \brief Get the number of bits written.
 * \param stream  bitstream
 * \return        position
 */
uint64_t bitstream_tell(const bitstream_t *const stream)
{
  uint64_t position = stream->len;
  return position * 8 + stream->cur_bit;
}

/**
 * \brief Write a byte to bitstream
 *
 * The stream must be byte-aligned.
 *
 * \param stream  pointer bitstream to put the data
 * \param byte    byte to write
 */
void bitstream_writebyte(bitstream_t *const stream, const uint8_t byte)
{
  assert(stream->cur_bit == 0);

  if (stream->last == NULL || stream->last->len == KVZ_DATA_CHUNK_SIZE) {
    // Need to allocate a new chunk.
    kvz_data_chunk *new_chunk = bitstream_alloc_chunk();
    assert(new_chunk);

    if (!stream->first) stream->first = new_chunk;
    if (stream->last)   stream->last->next = new_chunk;
    stream->last = new_chunk;
  }
  assert(stream->last->len < KVZ_DATA_CHUNK_SIZE);

  stream->last->data[stream->last->len] = byte;
  stream->last->len += 1;
  stream->len += 1;
}

/**
 * \brief Move data from one stream to another.
 *
 * Destination stream must be byte-aligned.
 *
 * Equivalent to bitstream_append(dst, src) followed by
 * bitstream_clear(src).
 */
void bitstream_move(bitstream_t *const dst, bitstream_t *const src)
{
  assert(dst->cur_bit == 0);

  if (src->len > 0) {
    if (dst->first == NULL) {
      dst->first = src->first;
      dst->last = src->last;
      dst->len = src->len;
    } else {
      dst->last->next = src->first;
      dst->last = src->last;
      dst->len += src->len;
    }
  }

  // Move the leftover bits.
  dst->data = src->data;
  dst->cur_bit = src->cur_bit;
  dst->zerocount = src->zerocount;

  src->first = src->last = NULL;
  bitstream_clear(src);
}

/**
 * \brief Copy data from one stream to another.
 *
 * Destination stream must be byte-aligned.
 */
void bitstream_append(bitstream_t *const dst, const bitstream_t *const src)
{
  assert(dst->cur_bit == 0);

  for (const kvz_data_chunk *chunk = src->first; chunk != NULL; chunk = chunk->next) {
    for (uint32_t i = 0; i < chunk->len; ++i) {
      bitstream_writebyte(dst, chunk->data[i]);
    }
  }

  // Copy the leftover bits.
  dst->data = src->data;
  dst->cur_bit = src->cur_bit;
  dst->zerocount = src->zerocount;
}

/**
 * Reset stream.
 */
void bitstream_clear(bitstream_t *const stream)
{
  bitstream_free_chunks(stream->first);
  bitstream_init(stream);
}

/**
 * \brief Write bits to bitstream
 * \param stream pointer bitstream to put the data
 * \param data input data
 * \param bits number of bits to write from data to stream
 */
void bitstream_put(bitstream_t *const stream, const uint32_t data, uint8_t bits)
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
      stream->cur_bit = 0;
      if((stream->zerocount == 2) && (stream->data < 4)) {
        bitstream_writebyte(stream, emulation_prevention_three_byte);
        stream->zerocount = 0;
      }
      if(stream->data == 0) {
        stream->zerocount++;
      } else {
        stream->zerocount = 0;
      }
      bitstream_writebyte(stream, stream->data);
    }
  }
}

/**
 * \brief Align the bitstream with one-bit padding
 */
void bitstream_align(bitstream_t * const stream)
{
  bitstream_put(stream, 1, 1);
  if ((stream->cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}

/**
 * \brief Align the bitstream with zero
 */
void bitstream_align_zero(bitstream_t * const stream)
{
  if ((stream->cur_bit & 7) != 0) {
    bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}
