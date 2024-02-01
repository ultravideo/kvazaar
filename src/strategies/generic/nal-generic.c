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

#include "strategies/generic/nal-generic.h"

#include "extras/libmd5.h"
#include "kvazaar.h"
#include "nal.h"
#include "strategyselector.h"


static void array_md5_generic(const kvz_pixel* data,
                              const int height, const int width,
                              const int stride,
                              unsigned char checksum_out[SEI_HASH_MAX_LENGTH], const uint8_t bitdepth)
{
  assert(SEI_HASH_MAX_LENGTH >= 16);

  context_md5_t md5_ctx;
  kvz_md5_init(&md5_ctx);
  
  unsigned bytes = width * height * sizeof(kvz_pixel);
  kvz_md5_update(&md5_ctx, (const unsigned char *)data, bytes);

  kvz_md5_final(checksum_out, &md5_ctx);
}

static void array_checksum_generic(const kvz_pixel* data,
                                   const int height, const int width,
                                   const int stride,
                                   unsigned char checksum_out[SEI_HASH_MAX_LENGTH], const uint8_t bitdepth) {
  int x, y;
  int checksum = 0;
  
  assert(SEI_HASH_MAX_LENGTH >= 4);
  
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      const uint8_t mask = (uint8_t)((x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8));
      checksum += (data[(y * stride) + x] & 0xff) ^ mask;
#if KVZ_BIT_DEPTH > 8
      checksum += ((data[(y * stride) + x] >> 8) & 0xff) ^ mask;
#endif
    }
  }

  // Unpack uint into byte-array.
  checksum_out[0] = (checksum >> 24) & 0xff;
  checksum_out[1] = (checksum >> 16) & 0xff;
  checksum_out[2] = (checksum >> 8) & 0xff;
  checksum_out[3] = (checksum) & 0xff;
}

static void array_checksum_generic4(const kvz_pixel* data,
                                   const int height, const int width,
                                   const int stride,
                                   unsigned char checksum_out[SEI_HASH_MAX_LENGTH], const uint8_t bitdepth) {
  uint32_t checksum = 0;
  int y, x, xp;
  
  static uint8_t ckmap_initialized = 0;
  static uint32_t ckmap[64*256];

  //TODO: add 10-bit support
  if(bitdepth != 8) {
    array_checksum_generic(data, height, width, stride,checksum_out, bitdepth);
    return;
  }
  if (!ckmap_initialized) {
    uint8_t * const ckmap_uint8 = (uint8_t*)&ckmap;
    int x, y;
    for (y = 0; y < 256; ++y) {
      for (x = 0; x < 256; ++x) {
        ckmap_uint8[y*256+x] = x^y;
      }
    }
    ckmap_initialized = 1;
  }

  assert(SEI_HASH_MAX_LENGTH >= 4);

  for (y = 0; y < height; ++y) {
    for (xp = 0; xp < width/4; ++xp) {
      const int x = xp * 4;
      const uint32_t mask = ckmap[(xp&63)+64*(y&255)] ^ (((x >> 8) ^ (y >> 8)) * 0x1010101);
      const uint32_t cksumbytes = (*((uint32_t*)(&data[(y * stride) + x]))) ^ mask;
      checksum += ((cksumbytes >> 24) & 0xff) + ((cksumbytes >> 16) & 0xff) + ((cksumbytes >> 8) & 0xff) + (cksumbytes & 0xff);
    }
    for (x = xp*4; x < width; ++x) {
      uint8_t mask = (uint8_t)((x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8));
      checksum += (data[(y * stride) + x] & 0xff) ^ mask;
    }
  }

  // Unpack uint into byte-array.
  checksum_out[0] = (checksum >> 24) & 0xff;
  checksum_out[1] = (checksum >> 16) & 0xff;
  checksum_out[2] = (checksum >> 8) & 0xff;
  checksum_out[3] = (checksum) & 0xff;
}

static void array_checksum_generic8(const kvz_pixel* data,
                                   const int height, const int width,
                                   const int stride,
                                   unsigned char checksum_out[SEI_HASH_MAX_LENGTH], const uint8_t bitdepth) {
  uint32_t checksum = 0;
  int y, x, xp;
  
  static uint8_t ckmap_initialized = 0;
  static uint64_t ckmap[32*256];
  
  //TODO: add 10-bit support
  if(bitdepth != 8) {
    array_checksum_generic(data, height, width, stride,checksum_out, bitdepth);
    return;
  }
  if (!ckmap_initialized) {
    uint8_t * const ckmap_uint8 = (uint8_t*)&ckmap;
    int x, y;
    for (y = 0; y < 256; ++y) {
      for (x = 0; x < 256; ++x) {
        ckmap_uint8[y*256+x] = x^y;
      }
    }
    ckmap_initialized = 1;
  }

  assert(SEI_HASH_MAX_LENGTH >= 4);

  for (y = 0; y < height; ++y) {
    if (y*stride % 8 != 0) {
      for (x = 0; x < width; ++x) {
        uint8_t mask = (uint8_t)((x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8));
        checksum += (data[(y * stride) + x] & 0xff) ^ mask;
      }
      continue;
    }
    for (xp = 0; xp < width/8; ++xp) {
      const int x = xp * 8;
      const uint64_t mask = ckmap[(xp&31)+32*(y&255)] ^ ((uint64_t)((x >> 8) ^ (y >> 8)) * 0x101010101010101);
      const uint64_t cksumbytes = (*((uint64_t*)(&data[(y * stride) + x]))) ^ mask;
      checksum += ((cksumbytes >> 56) & 0xff) + ((cksumbytes >> 48) & 0xff) + ((cksumbytes >> 40) & 0xff) + ((cksumbytes >> 32) & 0xff) + ((cksumbytes >> 24) & 0xff) + ((cksumbytes >> 16) & 0xff) + ((cksumbytes >> 8) & 0xff) + (cksumbytes & 0xff);
    }
    for (x = xp*8; x < width; ++x) {
      uint8_t mask = (uint8_t)((x & 0xff) ^ (y & 0xff) ^ (x >> 8) ^ (y >> 8));
      checksum += (data[(y * stride) + x] & 0xff) ^ mask;
    }
  }

  // Unpack uint into byte-array.
  checksum_out[0] = (checksum >> 24) & 0xff;
  checksum_out[1] = (checksum >> 16) & 0xff;
  checksum_out[2] = (checksum >> 8) & 0xff;
  checksum_out[3] = (checksum) & 0xff;
}

int kvz_strategy_register_nal_generic(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "array_md5", "generic", 0, &array_md5_generic);
  success &= kvz_strategyselector_register(opaque, "array_checksum", "generic", 0, &array_checksum_generic);
  success &= kvz_strategyselector_register(opaque, "array_checksum", "generic4", 1, &array_checksum_generic4);
  success &= kvz_strategyselector_register(opaque, "array_checksum", "generic8", 2, &array_checksum_generic8);
  
  return success;
}
