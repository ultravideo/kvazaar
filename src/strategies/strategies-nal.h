#ifndef STRATEGIES_NAL_H_
#define STRATEGIES_NAL_H_
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
 * \ingroup Optimization
 * \file
 * Interface for hash functions.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "nal.h"


//Function pointer to kvz_array_checksum
/**
 * \brief Calculate checksum for one color of the picture.
 * \param data Beginning of the pixel data for the picture.
 * \param height Height of the picture.
 * \param width Width of the picture.
 * \param stride Width of one row in the pixel array.
 */
typedef void (*array_checksum_func)(const kvz_pixel* data,
                                    const int height, const int width,
                                    const int stride,
                                    unsigned char checksum_out[SEI_HASH_MAX_LENGTH], const uint8_t bitdepth);
extern array_checksum_func kvz_array_checksum;
extern array_checksum_func kvz_array_md5;


int kvz_strategy_register_nal(void* opaque, uint8_t bitdepth);


#define STRATEGIES_NAL_EXPORTS \
  {"array_checksum", (void**) &kvz_array_checksum},\
  {"array_md5", (void**) &kvz_array_md5},

#endif //STRATEGIES_NAL_H_
