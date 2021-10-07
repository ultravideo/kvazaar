#ifndef TABLES_H_
#define TABLES_H_
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
 * \ingroup Reconstruction
 * \file
 * Various tables.
 */

#include "global.h" // IWYU pragma: keep


/***
  * List of indices for 4x4 coefficient groups within 8x8 transform block.
  * First index: 0 = diagonal, 1 = vertical, 2 horizontal scan pattern.
  * Second index: (log2 - 2) size of transform block. 4x4 .. 32x32
  */
static const uint32_t g_sig_last_scan_8x8[3][4] =
{ {0, 2, 1, 3},
  {0, 1, 2, 3},
  {0, 2, 1, 3}
};

/***
  * List of indices for 4x4 coefficient groups within 16x16 transform block.
  */
static const uint32_t g_sig_last_scan_16x16[16] = {
  0,  4,  1,  8,
  5,  2, 12,  9,
  6,  3, 13, 10,
  7, 14, 11, 15
};

/***
  * List of indices for 4x4 coefficient groups within 32x32 transform block.
  */
static const uint32_t g_sig_last_scan_32x32[64] = {
  0,   8,  1, 16,  9,  2, 24, 17,
  10,  3, 32, 25, 18, 11,  4, 40,
  33, 26, 19, 12,  5, 48, 41, 34,
  27, 20, 13,  6, 56, 49, 42, 35,
  28, 21, 14,  7, 57, 50, 43, 36,
  29, 22, 15, 58, 51, 44, 37, 30,
  23, 59, 52, 45, 38, 31, 60, 53,
  46, 39, 61, 54, 47, 62, 55, 63
};

/**
 * List of pointers to coefficient group mappings.
 * First index: (log2 - 2) of transform block size
 * Second index: scan pattern 0 = diagonal, 1 = horizontal, 2 = vertical
 */
static const uint32_t * const g_sig_last_scan_cg[4][3] = {
  { g_sig_last_scan_8x8[0], g_sig_last_scan_8x8[1], g_sig_last_scan_8x8[2] },  // 4x4, only first element is used
  { g_sig_last_scan_8x8[0], g_sig_last_scan_8x8[1], g_sig_last_scan_8x8[2] },
  { g_sig_last_scan_16x16, 0, 0 },
  { g_sig_last_scan_32x32, 0, 0 }
};


typedef enum
{
  SCAN_DIAG = 0, // up-right diagonal scan
  SCAN_HOR,      // horizontal first scan
  SCAN_VER       // vertical first scan
} coeff_scan_order_t;


/**
 * List of mappings for coefficients within a transform block.
 * First index: scan pattern 0 = diagonal, 1 = horizontal, 2 = vertical
 * Second index: (log2 - 1) size of transform block. 2x2 .. 32x32
 */
extern const uint32_t* const kvz_g_sig_last_scan[3][5];
extern const int8_t kvz_g_convert_to_bit[LCU_WIDTH + 1];

#endif //TABLES_H_
