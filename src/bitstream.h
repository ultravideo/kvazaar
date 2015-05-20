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

typedef enum {BITSTREAM_TYPE_FILE, BITSTREAM_TYPE_MEMORY} bitstream_type;

#define BASE_BITSTREAM uint8_t data; uint8_t cur_bit;  uint8_t zerocount; bitstream_type type;

//Size of the allocation for a memory bitstream in bytes
#define BITSTREAM_MEMORY_CHUNK_SIZE 4096

typedef struct
{
  BASE_BITSTREAM
} bitstream_base_t;

typedef struct
{
  BASE_BITSTREAM
  FILE*    output;
} bitstream_file_t;

typedef struct
{
  BASE_BITSTREAM
  uint8_t* output_data;
  uint32_t output_length;
  uint32_t allocated_length;
} bitstream_mem_t;

typedef union bitstream_t
{
  bitstream_base_t base;
  bitstream_file_t file;
  bitstream_mem_t mem;
} bitstream_t;

typedef struct
{
  uint8_t len;
  uint32_t value;
} bit_table_t;

extern const bit_table_t *g_exp_table;

int bitstream_init(bitstream_t * stream, bitstream_type type);
int bitstream_finalize(bitstream_t * stream);
void bitstream_put(bitstream_t *stream, uint32_t data, uint8_t bits);
int bitstream_writebyte(bitstream_t *stream_abstract, uint8_t byte);
long long unsigned int bitstream_tell(const bitstream_t * stream);

int bitstream_append(bitstream_t *dst, const bitstream_t *src);
int bitstream_clear(bitstream_t *stream);

/* Use macros to force inlining */
#define bitstream_put_ue(stream, data) { bitstream_put(stream,g_exp_table[data].value,g_exp_table[data].len); }
#define bitstream_put_se(stream, data) { uint32_t index=(uint32_t)(((data)<=0)?(-(data))<<1:((data)<<1)-1);    \
                                         bitstream_put(stream,g_exp_table[index].value,g_exp_table[index].len); }

void bitstream_align(bitstream_t *stream);
void bitstream_align_zero(bitstream_t *stream);
int init_exp_golomb(uint32_t len);
void free_exp_golomb();


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
