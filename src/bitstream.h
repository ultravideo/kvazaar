#ifndef BITSTREAM_H_
#define BITSTREAM_H_
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

/**
 * \ingroup CABAC
 * \file
 * Appending bits into an Annex-B coded bitstream.
 */

#include "global.h" // IWYU pragma: keep

#include "kvazaar.h"


/**
 * A stream of bits.
 */
typedef struct bitstream_t
{
  /// \brief Total number of complete bytes.
  uint32_t len;

  /// \brief Pointer to the first chunk, or NULL.
  kvz_data_chunk *first;

  /// \brief Pointer to the last chunk, or NULL.
  kvz_data_chunk *last;

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

void kvz_bitstream_init(bitstream_t * stream);
kvz_data_chunk * kvz_bitstream_alloc_chunk();
kvz_data_chunk * kvz_bitstream_take_chunks(bitstream_t *stream);
void kvz_bitstream_free_chunks(kvz_data_chunk *chunk);
void kvz_bitstream_finalize(bitstream_t * stream);

uint64_t kvz_bitstream_tell(const bitstream_t * stream);

void kvz_bitstream_writebyte(bitstream_t *stream, uint8_t byte);
void kvz_bitstream_move(bitstream_t *dst, bitstream_t *src);
void kvz_bitstream_clear(bitstream_t *stream);

void kvz_bitstream_put(bitstream_t *stream, uint32_t data, uint8_t bits);
void kvz_bitstream_put_byte(bitstream_t *const stream, const uint32_t data);

void kvz_bitstream_put_ue(bitstream_t *stream, uint32_t data);
void kvz_bitstream_put_se(bitstream_t *stream, int32_t data);

void kvz_bitstream_add_rbsp_trailing_bits(bitstream_t *stream);
void kvz_bitstream_align(bitstream_t *stream);
void kvz_bitstream_align_zero(bitstream_t *stream);

/* In debug mode print out some extra info */
#ifdef KVZ_DEBUG_PRINT_CABAC
/* Counter to keep up with bits written */
#define WRITE_U(stream, data, bits, name) { printf("%-40s u(%d) : %d\n", name,bits,data); kvz_bitstream_put(stream,data,bits);}
#define WRITE_UE(stream, data, name) { printf("%-40s ue(v): %d\n", name,data); kvz_bitstream_put_ue(stream,data);}
#define WRITE_SE(stream, data, name) { printf("%-40s se(v): %d\n", name,data); kvz_bitstream_put_se(stream,(data));}
#else
#define WRITE_U(stream, data, bits, name) { kvz_bitstream_put(stream,data,bits); }
#define WRITE_UE(stream, data, name) { kvz_bitstream_put_ue(stream,data); }
#define WRITE_SE(stream, data, name) { kvz_bitstream_put_se(stream,data); }
#endif


#endif
