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

#include "transform.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "nal.h"

//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
//

const int16_t g_t4[4][4] =
{
  { 64, 64, 64, 64},
  { 83, 36,-36,-83},
  { 64,-64,-64, 64},
  { 36,-83, 83,-36}
};

const int16_t g_t8[8][8] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64},
  { 89, 75, 50, 18,-18,-50,-75,-89},
  { 83, 36,-36,-83,-83,-36, 36, 83},
  { 75,-18,-89,-50, 50, 89, 18,-75},
  { 64,-64,-64, 64, 64,-64,-64, 64},
  { 50,-89, 18, 75,-75,-18, 89,-50},
  { 36,-83, 83,-36,-36, 83,-83, 36},
  { 18,-50, 75,-89, 89,-75, 50,-18}
};

const int16_t g_t16[16][16] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

const int16_t g_t32[32][32] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4, -4,-13,-22,-31,-38,-46,-54,-61,-67,-73,-78,-82,-85,-88,-90,-90},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90,-90,-87,-80,-70,-57,-43,-25, -9,  9, 25, 43, 57, 70, 80, 87, 90},
  { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31,  4,-22,-46,-67,-82,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89, 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22,-22,-61,-85,-90,-73,-38,  4, 46, 78, 90, 82, 54, 13,-31,-67,-88},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87,-87,-57, -9, 43, 80, 90, 70, 25,-25,-70,-90,-80,-43,  9, 57, 87},
  { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31, 31, 78, 90, 61,  4,-54,-88,-82,-38, 22, 73, 90, 67, 13,-46,-85},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38,-38,-88,-73, -4, 67, 90, 46,-31,-85,-78,-13, 61, 90, 54,-22,-82},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80,-80, -9, 70, 87, 25,-57,-90,-43, 43, 90, 57,-25,-87,-70,  9, 80},
  { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46, 46, 90, 38,-54,-90,-31, 61, 88, 22,-67,-85,-13, 73, 82,  4,-78},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75, 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54,-54,-85,  4, 88, 46,-61,-82, 13, 90, 38,-67,-78, 22, 90, 31,-73},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70,-70, 43, 87, -9,-90,-25, 80, 57,-57,-80, 25, 90,  9,-87,-43, 70},
  { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61, 61, 73,-46,-82, 31, 88,-13,-90, -4, 90, 22,-85,-38, 78, 54,-67},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67,-67,-54, 78, 38,-85,-22, 90,  4,-90, 13, 88,-31,-82, 46, 73,-61},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57,-57, 80, 25,-90,  9, 87,-43,-70, 70, 43,-87, -9, 90,-25,-80, 57},
  { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73, 73, 31,-90, 22, 78,-67,-38, 90,-13,-82, 61, 46,-88,  4, 85,-54},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50, 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78,-78, -4, 82,-73,-13, 85,-67,-22, 88,-61,-31, 90,-54,-38, 90,-46},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43,-43, 90,-57,-25, 87,-70, -9, 80,-80,  9, 70,-87, 25, 57,-90, 43},
  { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82, 82,-22,-54, 90,-61,-13, 78,-85, 31, 46,-90, 67,  4,-73, 88,-38},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85,-85, 46, 13,-67, 90,-73, 22, 38,-82, 88,-54, -4, 61,-90, 78,-31},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25,-25, 70,-90, 80,-43, -9, 57,-87, 87,-57,  9, 43,-80, 90,-70, 25},
  { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88, 88,-67, 31, 13,-54, 82,-90, 78,-46,  4, 38,-73, 90,-85, 61,-22},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18, 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90,-90, 82,-67, 46,-22, -4, 31,-54, 73,-85, 90,-88, 78,-61, 38,-13},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9, -9, 25,-43, 57,-70, 80,-87, 90,-90, 87,-80, 70,-57, 43,-25,  9},
  {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90, 90,-90, 88,-85, 82,-78, 73,-67, 61,-54, 46,-38, 31,-22, 13, -4}
};


const uint8_t g_chroma_scale[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//

/**
 * \brief Get scaled QP used in quantization
 *
 */
int32_t get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset)
{
  int32_t qp_scaled = 0;
  if(type == 0) {
    qp_scaled = qp + qp_offset;
  } else {
    qp_scaled = CLIP(-qp_offset, 57, qp);
    if(qp_scaled < 0) {
      qp_scaled = qp_scaled + qp_offset;
    } else {
      qp_scaled = g_chroma_scale[qp_scaled] + qp_offset;
    }
  }
  return qp_scaled;
}




static void partial_butterfly_4(short *src, short *dst,
                                int32_t shift, int32_t line)
{
  int32_t j;
  int32_t e[2],o[2];
  int32_t add = 1<<(shift - 1);

  for (j = 0; j < line; j++) {
    // E and O
    e[0] = src[0] + src[3];
    o[0] = src[0] - src[3];
    e[1] = src[1] + src[2];
    o[1] = src[1] - src[2];

    dst[0]      = (short)((g_t4[0][0]*e[0] + g_t4[0][1]*e[1] + add) >> shift);
    dst[2*line] = (short)((g_t4[2][0]*e[0] + g_t4[2][1]*e[1] + add) >> shift);
    dst[line]   = (short)((g_t4[1][0]*o[0] + g_t4[1][1]*o[1] + add) >> shift);
    dst[3*line] = (short)((g_t4[3][0]*o[0] + g_t4[3][1]*o[1] + add) >> shift);

    src += 4;
    dst ++;
  }
}

static void partial_butterfly_inverse_4(short *src,short *dst,
                                        int shift, int line)
{
  int j;
  int e[2],o[2];
  int add = 1<<(shift - 1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    o[0] = g_t4[1][0]*src[line] + g_t4[3][0]*src[3*line];
    o[1] = g_t4[1][1]*src[line] + g_t4[3][1]*src[3*line];
    e[0] = g_t4[0][0]*src[0]    + g_t4[2][0]*src[2*line];
    e[1] = g_t4[0][1]*src[0]    + g_t4[2][1]*src[2*line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    dst[0] = (short)CLIP(-32768, 32767, (e[0] + o[0] + add) >> shift);
    dst[1] = (short)CLIP(-32768, 32767, (e[1] + o[1] + add) >> shift);
    dst[2] = (short)CLIP(-32768, 32767, (e[1] - o[1] + add) >> shift);
    dst[3] = (short)CLIP(-32768, 32767, (e[0] - o[0] + add) >> shift);

    src++;
    dst += 4;
  }
}

// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm
// gives identical results
static void fast_forward_dst(short *block, short *coeff, int32_t shift)  // input block, output coeff
{
  int32_t i, c[4];
  int32_t rnd_factor = 1<<(shift - 1);
  for (i = 0; i < 4; i++) {
    // int32_termediate Variables
    c[0] = block[4*i + 0] + block[4*i + 3];
    c[1] = block[4*i + 1] + block[4*i + 3];
    c[2] = block[4*i + 0] - block[4*i + 1];
    c[3] = 74* block[4*i + 2];

    coeff[   i] =  (short)(( 29*c[0] + 55*c[1]         + c[3]                     + rnd_factor ) >> shift);
    coeff[ 4+i] =  (short)(( 74*(block[4*i + 0]+ block[4*i + 1] - block[4*i + 3]) + rnd_factor ) >> shift);
    coeff[ 8+i] =  (short)(( 29*c[2] + 55*c[0]         - c[3]                     + rnd_factor ) >> shift);
    coeff[12+i] =  (short)(( 55*c[2] - 29*c[1]         + c[3]                     + rnd_factor ) >> shift);
  }
}

static void fast_inverse_dst(short *tmp,short *block,int shift)  // input tmp, output block
{
  int i, c[4];
  int rnd_factor = 1<<(shift-1);
  for (i = 0; i < 4; i++) {
    // Intermediate Variables
    c[0] = tmp[    i] + tmp[ 8 + i];
    c[1] = tmp[8 + i] + tmp[12 + i];
    c[2] = tmp[    i] - tmp[12 + i];
    c[3] = 74 * tmp[4 + i];

    block[4*i + 0] = (short)CLIP(-32768, 32767, ( 29*c[0] + 55*c[1]     + c[3]            + rnd_factor ) >> shift);
    block[4*i + 1] = (short)CLIP(-32768, 32767, ( 55*c[2] - 29*c[1]     + c[3]            + rnd_factor ) >> shift);
    block[4*i + 2] = (short)CLIP(-32768, 32767, ( 74*(tmp[i] - tmp[8 + i]  + tmp[12 + i]) + rnd_factor ) >> shift);
    block[4*i + 3] = (short)CLIP(-32768, 32767, ( 55*c[0] + 29*c[2]     - c[3]            + rnd_factor ) >> shift);
  }
}


static void partial_butterfly_8(short *src, short *dst,
                                int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[4],o[4];
  int32_t ee[2],eo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {
    // E and O
    for (k = 0; k < 4; k++) {
      e[k] = src[k] + src[7 - k];
      o[k] = src[k] - src[7 - k];
    }
    // EE and EO
    ee[0] = e[0] + e[3];
    eo[0] = e[0] - e[3];
    ee[1] = e[1] + e[2];
    eo[1] = e[1] - e[2];

    dst[0]      = (short)((g_t8[0][0]*ee[0] + g_t8[0][1]*ee[1] + add) >> shift);
    dst[4*line] = (short)((g_t8[4][0]*ee[0] + g_t8[4][1]*ee[1] + add) >> shift);
    dst[2*line] = (short)((g_t8[2][0]*eo[0] + g_t8[2][1]*eo[1] + add) >> shift);
    dst[6*line] = (short)((g_t8[6][0]*eo[0] + g_t8[6][1]*eo[1] + add) >> shift);

    dst[line]   = (short)((g_t8[1][0]*o[0] + g_t8[1][1]*o[1] + g_t8[1][2]*o[2] + g_t8[1][3]*o[3] + add) >> shift);
    dst[3*line] = (short)((g_t8[3][0]*o[0] + g_t8[3][1]*o[1] + g_t8[3][2]*o[2] + g_t8[3][3]*o[3] + add) >> shift);
    dst[5*line] = (short)((g_t8[5][0]*o[0] + g_t8[5][1]*o[1] + g_t8[5][2]*o[2] + g_t8[5][3]*o[3] + add) >> shift);
    dst[7*line] = (short)((g_t8[7][0]*o[0] + g_t8[7][1]*o[1] + g_t8[7][2]*o[2] + g_t8[7][3]*o[3] + add) >> shift);

    src += 8;
    dst++;
  }
}

static void partial_butterfly_inverse_8(int16_t *src,int16_t *dst,
                                        int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[4],o[4];
  int32_t ee[2],eo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 4; k++) {
      o[k] = g_t8[ 1][k]*src[line] + g_t8[ 3][k]*src[3*line] + g_t8[ 5][k]*src[5*line] + g_t8[ 7][k]*src[7*line];
    }

    eo[0] = g_t8[2][0]*src[ 2*line ] + g_t8[6][0]*src[ 6*line ];
    eo[1] = g_t8[2][1]*src[ 2*line ] + g_t8[6][1]*src[ 6*line ];
    ee[0] = g_t8[0][0]*src[ 0      ] + g_t8[4][0]*src[ 4*line ];
    ee[1] = g_t8[0][1]*src[ 0      ] + g_t8[4][1]*src[ 4*line ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    e[0] = ee[0] + eo[0];
    e[3] = ee[0] - eo[0];
    e[1] = ee[1] + eo[1];
    e[2] = ee[1] - eo[1];
    for (k = 0; k < 4; k++) {
      dst[ k   ] = (int16_t)MAX(-32768, MIN(32767, (e[k] + o[k] + add)>>shift));
      dst[ k+4 ] = (int16_t)MAX(-32768, MIN(32767, (e[3-k] - o[3-k] + add)>>shift));
    }
    src++;
    dst += 8;
  }
}


static void partial_butterfly_16(short *src,short *dst,
                                 int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[8],o[8];
  int32_t ee[4],eo[4];
  int32_t eee[2],eeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {
    // E and O
    for (k = 0; k < 8; k++) {
      e[k] = src[k] + src[15 - k];
      o[k] = src[k] - src[15 - k];
    }
    // EE and EO
    for (k = 0; k < 4; k++) {
      ee[k] = e[k] + e[7 - k];
      eo[k] = e[k] - e[7 - k];
    }
    // EEE and EEO
    eee[0] = ee[0] + ee[3];
    eeo[0] = ee[0] - ee[3];
    eee[1] = ee[1] + ee[2];
    eeo[1] = ee[1] - ee[2];

    dst[0      ] = (short)((g_t16[ 0][0]*eee[0] + g_t16[ 0][1]*eee[1] + add) >> shift);
    dst[8*line ] = (short)((g_t16[ 8][0]*eee[0] + g_t16[ 8][1]*eee[1] + add) >> shift);
    dst[4*line ] = (short)((g_t16[ 4][0]*eeo[0] + g_t16[ 4][1]*eeo[1] + add) >> shift);
    dst[12*line] = (short)((g_t16[12][0]*eeo[0] + g_t16[12][1]*eeo[1] + add) >> shift);

    for (k = 2; k < 16; k += 4) {
      dst[k*line] = (short)((g_t16[k][0]*eo[0] + g_t16[k][1]*eo[1] + g_t16[k][2]*eo[2] + g_t16[k][3]*eo[3] + add) >> shift);
    }

    for (k = 1; k < 16; k += 2) {
      dst[k*line] = (short)((g_t16[k][0]*o[0] + g_t16[k][1]*o[1] + g_t16[k][2]*o[2] + g_t16[k][3]*o[3] +
                     g_t16[k][4]*o[4] + g_t16[k][5]*o[5] + g_t16[k][6]*o[6] + g_t16[k][7]*o[7] + add) >> shift);
    }

    src += 16;
    dst++;
  }
}


static void partial_butterfly_inverse_16(int16_t *src, int16_t *dst,
                                         int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[8],o[8];
  int32_t ee[4],eo[4];
  int32_t eee[2],eeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 8; k++)  {
      o[k] = g_t16[ 1][k]*src[  line] + g_t16[ 3][k]*src[ 3*line] + g_t16[ 5][k]*src[ 5*line] + g_t16[ 7][k]*src[ 7*line] +
             g_t16[ 9][k]*src[9*line] + g_t16[11][k]*src[11*line] + g_t16[13][k]*src[13*line] + g_t16[15][k]*src[15*line];
    }
    for (k = 0; k < 4; k++) {
      eo[k] = g_t16[ 2][k]*src[ 2*line] + g_t16[ 6][k]*src[ 6*line] + g_t16[10][k]*src[10*line] + g_t16[14][k]*src[14*line];
    }
    eeo[0] = g_t16[4][0]*src[ 4*line ] + g_t16[12][0]*src[ 12*line ];
    eee[0] = g_t16[0][0]*src[ 0      ] + g_t16[ 8][0]*src[ 8*line  ];
    eeo[1] = g_t16[4][1]*src[ 4*line ] + g_t16[12][1]*src[ 12*line ];
    eee[1] = g_t16[0][1]*src[ 0      ] + g_t16[ 8][1]*src[ 8*line  ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    for (k = 0; k < 2; k++) {
      ee[k]   = eee[k]   + eeo[k];
      ee[k+2] = eee[1-k] - eeo[1-k];
    }
    for (k = 0; k < 4; k++) {
      e[k]   = ee[k] + eo[k];
      e[k+4] = ee[3-k] - eo[3-k];
    }
    for (k = 0; k < 8; k++) {
      dst[k]   = (short)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k+8] = (short)MAX(-32768, MIN(32767, (e[7-k] - o[7-k] + add) >> shift));
    }
    src++;
    dst += 16;
  }
}



static void partial_butterfly_32(short *src, short *dst,
                                 int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[16],o[16];
  int32_t ee[8],eo[8];
  int32_t eee[4],eeo[4];
  int32_t eeee[2],eeeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {
    // E and O
    for (k = 0; k < 16; k++) {
      e[k] = src[k] + src[31-k];
      o[k] = src[k] - src[31-k];
    }
    // EE and EO
    for (k = 0; k < 8; k++) {
      ee[k] = e[k] + e[15-k];
      eo[k] = e[k] - e[15-k];
    }
    // EEE and EEO
    for (k = 0; k < 4; k++) {
      eee[k] = ee[k] + ee[7-k];
      eeo[k] = ee[k] - ee[7-k];
    }
    // EEEE and EEEO
    eeee[0] = eee[0] + eee[3];
    eeeo[0] = eee[0] - eee[3];
    eeee[1] = eee[1] + eee[2];
    eeeo[1] = eee[1] - eee[2];

    dst[0      ] = (short)((g_t32[ 0][0]*eeee[0] + g_t32[ 0][1]*eeee[1] + add) >> shift);
    dst[16*line] = (short)((g_t32[16][0]*eeee[0] + g_t32[16][1]*eeee[1] + add) >> shift);
    dst[ 8*line] = (short)((g_t32[ 8][0]*eeeo[0] + g_t32[ 8][1]*eeeo[1] + add) >> shift);
    dst[24*line] = (short)((g_t32[24][0]*eeeo[0] + g_t32[24][1]*eeeo[1] + add) >> shift);
    for (k = 4; k < 32; k += 8) {
      dst[k*line] = (short)((g_t32[k][0]*eeo[0] + g_t32[k][1]*eeo[1] + g_t32[k][2]*eeo[2] + g_t32[k][3]*eeo[3] + add) >> shift);
    }
    for (k = 2; k < 32; k += 4) {
      dst[k*line] = (short)((g_t32[k][0]*eo[0] + g_t32[k][1]*eo[1] + g_t32[k][2]*eo[2] + g_t32[k][3]*eo[3] +
                             g_t32[k][4]*eo[4] + g_t32[k][5]*eo[5] + g_t32[k][6]*eo[6] + g_t32[k][7]*eo[7] + add) >> shift);
    }
    for (k = 1; k < 32; k += 2) {
      dst[k*line] = (short)((g_t32[k][ 0]*o[ 0] + g_t32[k][ 1]*o[ 1] + g_t32[k][ 2]*o[ 2] + g_t32[k][ 3]*o[ 3] +
                             g_t32[k][ 4]*o[ 4] + g_t32[k][ 5]*o[ 5] + g_t32[k][ 6]*o[ 6] + g_t32[k][ 7]*o[ 7] +
                             g_t32[k][ 8]*o[ 8] + g_t32[k][ 9]*o[ 9] + g_t32[k][10]*o[10] + g_t32[k][11]*o[11] +
                             g_t32[k][12]*o[12] + g_t32[k][13]*o[13] + g_t32[k][14]*o[14] + g_t32[k][15]*o[15] + add) >> shift);
    }
    src += 32;
    dst++;
  }
}


static void partial_butterfly_inverse_32(int16_t *src, int16_t *dst,
                                         int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[16],o[16];
  int32_t ee[8],eo[8];
  int32_t eee[4],eeo[4];
  int32_t eeee[2],eeeo[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 16; k++) {
      o[k] = g_t32[ 1][k]*src[ line  ]   + g_t32[ 3][k]*src[ 3*line  ] + g_t32[ 5][k]*src[ 5*line  ] + g_t32[ 7][k]*src[ 7*line  ] +
             g_t32[ 9][k]*src[ 9*line  ] + g_t32[11][k]*src[ 11*line ] + g_t32[13][k]*src[ 13*line ] + g_t32[15][k]*src[ 15*line ] +
             g_t32[17][k]*src[ 17*line ] + g_t32[19][k]*src[ 19*line ] + g_t32[21][k]*src[ 21*line ] + g_t32[23][k]*src[ 23*line ] +
             g_t32[25][k]*src[ 25*line ] + g_t32[27][k]*src[ 27*line ] + g_t32[29][k]*src[ 29*line ] + g_t32[31][k]*src[ 31*line ];
    }
    for (k = 0; k < 8; k++) {
      eo[k] = g_t32[ 2][k]*src[ 2*line  ] + g_t32[ 6][k]*src[ 6*line  ] + g_t32[10][k]*src[ 10*line ] + g_t32[14][k]*src[ 14*line ] +
              g_t32[18][k]*src[ 18*line ] + g_t32[22][k]*src[ 22*line ] + g_t32[26][k]*src[ 26*line ] + g_t32[30][k]*src[ 30*line ];
    }
    for (k = 0; k < 4; k++) {
      eeo[k] = g_t32[4][k]*src[ 4*line ] + g_t32[12][k]*src[ 12*line ] + g_t32[20][k]*src[ 20*line ] + g_t32[28][k]*src[ 28*line ];
    }
    eeeo[0] = g_t32[8][0]*src[ 8*line ] + g_t32[24][0]*src[ 24*line ];
    eeeo[1] = g_t32[8][1]*src[ 8*line ] + g_t32[24][1]*src[ 24*line ];
    eeee[0] = g_t32[0][0]*src[ 0      ] + g_t32[16][0]*src[ 16*line ];
    eeee[1] = g_t32[0][1]*src[ 0      ] + g_t32[16][1]*src[ 16*line ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    eee[0] = eeee[0] + eeeo[0];
    eee[3] = eeee[0] - eeeo[0];
    eee[1] = eeee[1] + eeeo[1];
    eee[2] = eeee[1] - eeeo[1];
    for (k = 0; k < 4; k++) {
      ee[k]   = eee[k]   + eeo[k];
      ee[k+4] = eee[3-k] - eeo[3-k];
    }
    for (k = 0; k < 8; k++) {
      e[k]   = ee[k]   + eo[k];
      e[k+8] = ee[7-k] - eo[7-k];
    }
    for (k=0;k<16;k++) {
      dst[k]    = (short)MAX( -32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k+16] = (short)MAX( -32768, MIN(32767, (e[15-k] - o[15-k] + add) >> shift));
    }
    src++;
    dst += 32;
  }
}

/**
 * \brief NxN inverse transform (2D)
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size input data (width of transform)
 */
void transformskip(const encoder_control * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  g_convert_to_bit[block_size] + 2;
  int32_t  shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size;
  int32_t  j,k;
  for (j = 0; j < block_size; j++) {
    for(k = 0; k < block_size; k ++) {
      coeff[j * block_size + k] = block[j * block_size + k] << shift;
    }
  }
}

/**
 * \brief inverse transform skip
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size width of transform
 */
void itransformskip(const encoder_control * const encoder, int16_t *block,int16_t *coeff, int8_t block_size)
{
  uint32_t log2_tr_size =  g_convert_to_bit[block_size] + 2;
  int32_t  shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size;
  int32_t  j,k;
  int32_t offset;
  offset = (1 << (shift -1)); // For rounding
  for ( j = 0; j < block_size; j++ ) {
    for(k = 0; k < block_size; k ++) {
      block[j * block_size + k] =  (coeff[j * block_size + k] + offset) >> shift;
    }
  }
}

/**
 * \brief forward transform (2D)
 * \param block input residual
 * \param coeff transform coefficients
 * \param block_size width of transform
 */
void transform2d(const encoder_control * const encoder, int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode)
{
  int32_t shift_1st = g_convert_to_bit[block_size]  + 1 + (encoder->bitdepth - 8);
  int32_t shift_2nd = g_convert_to_bit[block_size]  + 8;

  int16_t tmp[LCU_WIDTH * LCU_WIDTH];

  if(block_size== 4) {
    if (mode != 65535) {
      // Forward DST BY FAST ALGORITHM
      fast_forward_dst(block,tmp,shift_1st);
      fast_forward_dst(tmp,coeff,shift_2nd);
    } else {
      partial_butterfly_4(block, tmp, shift_1st, block_size);
      partial_butterfly_4(tmp, coeff, shift_2nd, block_size);
    }
  } else {
    switch(block_size) {
      case 8:
        partial_butterfly_8( block, tmp, shift_1st, block_size );
        partial_butterfly_8( tmp, coeff, shift_2nd, block_size );
        break;
      case 16:
        partial_butterfly_16( block, tmp, shift_1st, block_size );
        partial_butterfly_16( tmp, coeff, shift_2nd, block_size );
        break;
      case 32:
        partial_butterfly_32( block, tmp, shift_1st, block_size );
        partial_butterfly_32( tmp, coeff, shift_2nd, block_size );
        break;
    }
  }
}

/**
 * \brief NxN inverse transform (2D)
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size input data (width of transform)
 * \param mode
 */
void itransform2d(const encoder_control * const encoder,int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (encoder->bitdepth - 8);
  int16_t tmp[LCU_WIDTH*LCU_WIDTH];

  if( block_size == 4) {
    if (mode != 65535) {
      // Inverse DST by FAST Algorithm
      fast_inverse_dst(coeff, tmp, shift_1st);
      fast_inverse_dst(tmp, block, shift_2nd);
    } else {
      partial_butterfly_inverse_4(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_4(tmp, block, shift_2nd, block_size);
    }
  } else {
    switch(block_size) {
    case 8:
      partial_butterfly_inverse_8(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_8(tmp, block, shift_2nd, block_size);
      break;
    case 16:
      partial_butterfly_inverse_16(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_16(tmp, block, shift_2nd, block_size);
      break;
    case 32:
      partial_butterfly_inverse_32(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_32(tmp, block, shift_2nd, block_size);
      break;
    }
  }
}


#define QUANT_SHIFT 14
/**
 * \brief quantize transformed coefficents
 *
 */
void quant(const encoder_state * const encoder_state, int16_t *coef, int16_t *q_coef, int32_t width,
           int32_t height, uint32_t *ac_sum, int8_t type, int8_t scan_idx, int8_t block_type )
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  const uint32_t log2_block_size = g_convert_to_bit[ width ] + 2;
  const uint32_t * const scan = g_sig_last_scan[ scan_idx ][ log2_block_size - 1 ];

  #if ENABLE_SIGN_HIDING == 1
  int32_t delta_u[LCU_WIDTH*LCU_WIDTH>>2];
  #endif

  int32_t qp_scaled = get_scaled_qp(type, encoder_state->global->QP, 0);

  //New block for variable definitions
  {
  int32_t n;
  uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
  int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

  const int32_t *quant_coeff = encoder->scaling_list.quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];

  int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  int32_t q_bits = QUANT_SHIFT + qp_scaled/6 + transform_shift;
  int32_t add = ((encoder_state->tile->cur_pic->slicetype == SLICE_I) ? 171 : 85) << (q_bits - 9);

  int32_t q_bits8 = q_bits - 8;
  for (n = 0; n < width * height; n++) {
    int32_t level;
    int32_t  sign;

    level = coef[n];
    sign  = (level < 0 ? -1: 1);

    level = ((int64_t)abs(level) * quant_coeff[n] + add ) >> q_bits;

    #if ENABLE_SIGN_HIDING == 1
    delta_u[n] = (int32_t)( ((int64_t)abs(coef[n]) * quant_coeff[n] - (level<<q_bits) )>> q_bits8 );
    *ac_sum += level;
    #endif

    level *= sign;
    q_coef[n] = (int16_t)(CLIP( -32768, 32767, level));
  }

  #if ENABLE_SIGN_HIDING == 1
  if(*ac_sum >= 2) {
    #define SCAN_SET_SIZE 16
    #define LOG2_SCAN_SET_SIZE 4
    int32_t n,last_cg = -1, abssum = 0, subset, subpos;
    for(subset = (width*height - 1)>>LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg=-1;
      subpos = subset<<LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n = SCAN_SET_SIZE - 1; n >= 0; n--)  {
        if (q_coef[scan[n + subpos]])  {
          last_nz_pos_in_cg = n;
          break;
        }
      }

      // First coeff pos
      for (n = 0; n <SCAN_SET_SIZE; n++) {
        if (q_coef[scan[n + subpos]]) {
          first_nz_pos_in_cg = n;
          break;
        }
      }

      // Sum all quant coeffs between first and last
      for(n = first_nz_pos_in_cg; n <= last_nz_pos_in_cg; n++) {
        abssum += q_coef[scan[n + subpos]];
      }

      if(last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if(last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1) ;
        if(signbit != (abssum&0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff,  min_pos =-1, cur_cost=0x7fffffff;
          int16_t final_change = 0, cur_change=0;
          for(n = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n >= 0; n--) {
            uint32_t blkPos  = scan[n + subpos];
            if(q_coef[blkPos] != 0) {
              if(delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change=1;
              } else if(n == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost=0x7fffffff;
              } else {
                cur_cost = delta_u[blkPos];
                cur_change =-1;
              }
            } else if(n < first_nz_pos_in_cg && ((coef[blkPos] >= 0)?0:1) != signbit) {
              cur_cost = 0x7fffffff;
            } else {
              cur_cost   = -delta_u[blkPos];
              cur_change = 1;
            }

            if(cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos      = blkPos;
            }
          } // CG loop

          if(q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if(coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;
        } // Hide
      }
      if (last_cg == 1) last_cg=0;
    }

    #undef SCAN_SET_SIZE
    #undef LOG2_SCAN_SET_SIZE
  }
  #endif
  }
}

/**
 * \brief inverse quantize transformed and quantized coefficents
 *
 */
void dequant(const encoder_state * const encoder_state, int16_t *q_coef, int16_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  int32_t shift,add,coeff_q;
  int32_t n;
  int32_t transform_shift = 15 - encoder->bitdepth - (g_convert_to_bit[ width ] + 2);

  int32_t qp_scaled = get_scaled_qp(type, encoder_state->global->QP, 0);

  shift = 20 - QUANT_SHIFT - transform_shift;

  if (encoder->scaling_list.enable)
  {
    uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    const int32_t *dequant_coef = encoder->scaling_list.de_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
    shift += 4;

    if (shift >qp_scaled / 6) {
      add = 1 << (shift - qp_scaled/6 - 1);

      for (n = 0; n < width * height; n++) {
        coeff_q = ((q_coef[n] * dequant_coef[n]) + add ) >> (shift -  qp_scaled/6);
        coef[n] = (int16_t)CLIP(-32768,32767,coeff_q);
      }
    } else {
      for (n = 0; n < width * height; n++) {
        // Clip to avoid possible overflow in following shift left operation
        coeff_q   = CLIP(-32768, 32767, q_coef[n] * dequant_coef[n]);
        coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q << (qp_scaled/6 - shift));
      }
    }
  } else {
    int32_t scale = g_inv_quant_scales[qp_scaled%6] << (qp_scaled/6);
    add = 1 << (shift-1);

    for (n = 0; n < width*height; n++) {
      coeff_q   = (q_coef[n] * scale + add) >> shift;
      coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q);
    }
  }
}

