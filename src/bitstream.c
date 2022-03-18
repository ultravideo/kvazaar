/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "bitstream.h"

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "kvz_math.h"


const uint32_t kvz_bit_set_mask[] =
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

/**
 * \brief Initialize a new bitstream.
 */
void kvz_bitstream_init(bitstream_t *const stream)
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
kvz_data_chunk * kvz_bitstream_take_chunks(bitstream_t *const stream)
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
kvz_data_chunk * kvz_bitstream_alloc_chunk()
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
void kvz_bitstream_free_chunks(kvz_data_chunk *chunk)
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
void kvz_bitstream_finalize(bitstream_t *const stream)
{
  kvz_bitstream_clear(stream);
}

/**
 * \brief Get the number of bits written.
 * \param stream  bitstream
 * \return        position
 */
uint64_t kvz_bitstream_tell(const bitstream_t *const stream)
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
void kvz_bitstream_writebyte(bitstream_t *const stream, const uint8_t byte)
{
  assert(stream->cur_bit == 0);

  if (stream->last == NULL || stream->last->len == KVZ_DATA_CHUNK_SIZE) {
    // Need to allocate a new chunk.
    kvz_data_chunk *new_chunk = kvz_bitstream_alloc_chunk();
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
 * Destination stream must be byte-aligned. Source stream will be cleared.
 */
void kvz_bitstream_move(bitstream_t *const dst, bitstream_t *const src)
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
  kvz_bitstream_clear(src);
}

/**
 * Reset stream.
 */
void kvz_bitstream_clear(bitstream_t *const stream)
{
  kvz_bitstream_free_chunks(stream->first);
  kvz_bitstream_init(stream);
}

/**
 * \brief Write a byte to a byte aligned bitstream
 * \param stream  stream the data is to be appended to
 * \param data  input data
 */
void kvz_bitstream_put_byte(bitstream_t *const stream, uint32_t data)
{
  assert(stream->cur_bit == 0);
  const uint8_t emulation_prevention_three_byte = 0x03;

  if ((stream->zerocount == 2) && (data < 4)) {
    kvz_bitstream_writebyte(stream, emulation_prevention_three_byte);
    stream->zerocount = 0;
  }
  stream->zerocount = data == 0 ? stream->zerocount + 1 : 0;
  kvz_bitstream_writebyte(stream, data);
}

/**
 * \brief Write bits to bitstream
 *        Buffers individual bits untill they make a full byte.
 * \param stream  stream the data is to be appended to
 * \param data  input data
 * \param bits  number of bits to write from data to stream
 */
void kvz_bitstream_put(bitstream_t *const stream, const uint32_t data, uint8_t bits)
{
  while (bits--) {
    stream->data <<= 1;

    if (data & kvz_bit_set_mask[bits]) {
      stream->data |= 1;
    }
    stream->cur_bit++;

    // write byte to output
    if (stream->cur_bit == 8) {
      stream->cur_bit = 0;
      kvz_bitstream_put_byte(stream, stream->data);
    }
  }
}

/**
 * \brief Write unsigned Exp-Golomb bit string
 */
void kvz_bitstream_put_ue(bitstream_t *stream, uint32_t code_num)
{
  unsigned code_num_log2 = kvz_math_floor_log2(code_num + 1);
  unsigned prefix = 1 << code_num_log2;
  unsigned suffix = code_num + 1 - prefix;
  unsigned num_bits = code_num_log2 * 2 + 1;
  unsigned value = prefix | suffix;

  kvz_bitstream_put(stream, value, num_bits);
}

/**
 * \brief Write signed Exp-Golomb bit string
 */
void kvz_bitstream_put_se(bitstream_t *stream, int32_t data)
{
  // Map positive values to even and negative to odd values.
  uint32_t code_num = data <= 0 ? (-data) << 1 : (data << 1) - 1;
  kvz_bitstream_put_ue(stream, code_num);
}

/**
 * \brief Add rbsp_trailing_bits syntax element, which aligns the bitstream.
 */
void kvz_bitstream_add_rbsp_trailing_bits(bitstream_t * const stream)
{
  kvz_bitstream_put(stream, 1, 1);
  if ((stream->cur_bit & 7) != 0) {
    kvz_bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}

/**
* \brief Align the bitstream, unless it's already aligned.
*/
void kvz_bitstream_align(bitstream_t * const stream)
{
  if ((stream->cur_bit & 7) != 0) {
    kvz_bitstream_add_rbsp_trailing_bits(stream);
  }
}

/**
 * \brief Align the bitstream with zero
 */
void kvz_bitstream_align_zero(bitstream_t * const stream)
{
  if ((stream->cur_bit & 7) != 0) {
    kvz_bitstream_put(stream, 0, 8 - (stream->cur_bit & 7));
  }
}
