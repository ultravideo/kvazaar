#ifndef BITSTREAM_H_
#define BITSTREAM_H_
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
 * \brief Bitstream can be written to one or several bits at a time.
 */

#include "global.h"

//Size of the allocation for a memory bitstream in bytes
#define BITSTREAM_MEMORY_CHUNK_SIZE 4096

/**
 * \brief A list of chunks of data in a bitstream.
 */
typedef struct bitstream_chunk_t {
  /// \brief Buffer for the data.
  uint8_t data[BITSTREAM_MEMORY_CHUNK_SIZE];

  /// \brief Number of bytes filled in this chunk.
  uint32_t len;

  /// \brief Next chunk in the list.
  struct bitstream_chunk_t *next;
} bitstream_chunk_t;

/**
 * A stream of bits.
 */
typedef struct bitstream_t
{
  /// \brief Total number of complete bytes.
  uint32_t len;

  /// \brief Pointer to the first chunk, or NULL.
  bitstream_chunk_t *first;

  /// \brief Pointer to the last chunk, or NULL.
  bitstream_chunk_t *last;

  /// \brief The incomplete byte.
  uint8_t data;

  /// \brief Number of bits in the incomplete byte.
  uint8_t cur_bit;

  uint8_t zerocount;
} bitstream_t;

typedef struct
{
  uint8_t len;
  uint32_t value;
} bit_table_t;

extern bit_table_t g_exp_table[EXP_GOLOMB_TABLE_SIZE];

void init_exp_golomb();

void bitstream_init(bitstream_t * stream);
bitstream_chunk_t * bitstream_alloc_chunk();
bitstream_chunk_t * bitstream_take_chunks(bitstream_t *stream);
void bitstream_free_chunks(bitstream_chunk_t *chunk);
void bitstream_finalize(bitstream_t * stream);

uint64_t bitstream_tell(const bitstream_t * stream);

void bitstream_writebyte(bitstream_t *stream, uint8_t byte);
void bitstream_move(bitstream_t *dst, bitstream_t *src);
void bitstream_append(bitstream_t *dst, const bitstream_t *src);
void bitstream_clear(bitstream_t *stream);

void bitstream_put(bitstream_t *stream, uint32_t data, uint8_t bits);
/* Use macros to force inlining */
#define bitstream_put_ue(stream, data) { bitstream_put(stream,g_exp_table[data].value,g_exp_table[data].len); }
#define bitstream_put_se(stream, data) { uint32_t index=(uint32_t)(((data)<=0)?(-(data))<<1:((data)<<1)-1);    \
                                         bitstream_put(stream,g_exp_table[index].value,g_exp_table[index].len); }

void bitstream_align(bitstream_t *stream);
void bitstream_align_zero(bitstream_t *stream);

/* In debug mode print out some extra info */
#ifdef NOTDEFINED//_DEBUG
/* Counter to keep up with bits written */
#define WRITE_U(stream, data, bits, name) { printf("%-40s u(%d) : %d\n", name,bits,data); bitstream_put(stream,data,bits);}
#define WRITE_UE(stream, data, name) { printf("%-40s ue(v): %d\n", name,data); bitstream_put_ue(stream,data);}
#define WRITE_SE(stream, data, name) { printf("%-40s se(v): %d\n", name,data); bitstream_put_se(stream,(data));}
#else
#define WRITE_U(stream, data, bits, name) { bitstream_put(stream,data,bits); }
#define WRITE_UE(stream, data, name) { bitstream_put_ue(stream,data); }
#define WRITE_SE(stream, data, name) { bitstream_put_se(stream,data); }
#endif


#endif
