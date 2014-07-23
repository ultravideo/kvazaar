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

#include <stdlib.h>

#include "strategyselector.h"

extern const int16_t g_t4[4][4];
extern const int16_t g_t8[8][8];
extern const int16_t g_t16[16][16];
extern const int16_t g_t32[32][32];

/**
 * \brief Generic partial butterfly functions
 *
 * TODO: description
 *
 * \param TODO   
 *
 * \returns TODO
 */


static void partial_butterfly_4_generic(short *src, short *dst,
  int32_t shift, int32_t line)
{
  int32_t j;
  int32_t e[2], o[2];
  int32_t add = 1 << (shift - 1);

  for (j = 0; j < line; j++) {
    // E and O
    e[0] = src[0] + src[3];
    o[0] = src[0] - src[3];
    e[1] = src[1] + src[2];
    o[1] = src[1] - src[2];

    dst[0] = (short)((g_t4[0][0] * e[0] + g_t4[0][1] * e[1] + add) >> shift);
    dst[2 * line] = (short)((g_t4[2][0] * e[0] + g_t4[2][1] * e[1] + add) >> shift);
    dst[line] = (short)((g_t4[1][0] * o[0] + g_t4[1][1] * o[1] + add) >> shift);
    dst[3 * line] = (short)((g_t4[3][0] * o[0] + g_t4[3][1] * o[1] + add) >> shift);

    src += 4;
    dst++;
  }
}


static void partial_butterfly_inverse_4_generic(short *src, short *dst,
  int shift, int line)
{
  int j;
  int e[2], o[2];
  int add = 1 << (shift - 1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    o[0] = g_t4[1][0] * src[line] + g_t4[3][0] * src[3 * line];
    o[1] = g_t4[1][1] * src[line] + g_t4[3][1] * src[3 * line];
    e[0] = g_t4[0][0] * src[0] + g_t4[2][0] * src[2 * line];
    e[1] = g_t4[0][1] * src[0] + g_t4[2][1] * src[2 * line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    dst[0] = (short)CLIP(-32768, 32767, (e[0] + o[0] + add) >> shift);
    dst[1] = (short)CLIP(-32768, 32767, (e[1] + o[1] + add) >> shift);
    dst[2] = (short)CLIP(-32768, 32767, (e[1] - o[1] + add) >> shift);
    dst[3] = (short)CLIP(-32768, 32767, (e[0] - o[0] + add) >> shift);

    src++;
    dst += 4;
  }
}


static void partial_butterfly_8_generic(short *src, short *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[4], o[4];
  int32_t ee[2], eo[2];
  int32_t add = 1 << (shift - 1);

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

    dst[0] = (short)((g_t8[0][0] * ee[0] + g_t8[0][1] * ee[1] + add) >> shift);
    dst[4 * line] = (short)((g_t8[4][0] * ee[0] + g_t8[4][1] * ee[1] + add) >> shift);
    dst[2 * line] = (short)((g_t8[2][0] * eo[0] + g_t8[2][1] * eo[1] + add) >> shift);
    dst[6 * line] = (short)((g_t8[6][0] * eo[0] + g_t8[6][1] * eo[1] + add) >> shift);

    dst[line] = (short)((g_t8[1][0] * o[0] + g_t8[1][1] * o[1] + g_t8[1][2] * o[2] + g_t8[1][3] * o[3] + add) >> shift);
    dst[3 * line] = (short)((g_t8[3][0] * o[0] + g_t8[3][1] * o[1] + g_t8[3][2] * o[2] + g_t8[3][3] * o[3] + add) >> shift);
    dst[5 * line] = (short)((g_t8[5][0] * o[0] + g_t8[5][1] * o[1] + g_t8[5][2] * o[2] + g_t8[5][3] * o[3] + add) >> shift);
    dst[7 * line] = (short)((g_t8[7][0] * o[0] + g_t8[7][1] * o[1] + g_t8[7][2] * o[2] + g_t8[7][3] * o[3] + add) >> shift);

    src += 8;
    dst++;
  }
}


static void partial_butterfly_inverse_8_generic(int16_t *src, int16_t *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[4], o[4];
  int32_t ee[2], eo[2];
  int32_t add = 1 << (shift - 1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 4; k++) {
      o[k] = g_t8[1][k] * src[line] + g_t8[3][k] * src[3 * line] + g_t8[5][k] * src[5 * line] + g_t8[7][k] * src[7 * line];
    }

    eo[0] = g_t8[2][0] * src[2 * line] + g_t8[6][0] * src[6 * line];
    eo[1] = g_t8[2][1] * src[2 * line] + g_t8[6][1] * src[6 * line];
    ee[0] = g_t8[0][0] * src[0] + g_t8[4][0] * src[4 * line];
    ee[1] = g_t8[0][1] * src[0] + g_t8[4][1] * src[4 * line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    e[0] = ee[0] + eo[0];
    e[3] = ee[0] - eo[0];
    e[1] = ee[1] + eo[1];
    e[2] = ee[1] - eo[1];
    for (k = 0; k < 4; k++) {
      dst[k] = (int16_t)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k + 4] = (int16_t)MAX(-32768, MIN(32767, (e[3 - k] - o[3 - k] + add) >> shift));
    }
    src++;
    dst += 8;
  }
}


static void partial_butterfly_16_generic(short *src, short *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[8], o[8];
  int32_t ee[4], eo[4];
  int32_t eee[2], eeo[2];
  int32_t add = 1 << (shift - 1);

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

    dst[0] = (short)((g_t16[0][0] * eee[0] + g_t16[0][1] * eee[1] + add) >> shift);
    dst[8 * line] = (short)((g_t16[8][0] * eee[0] + g_t16[8][1] * eee[1] + add) >> shift);
    dst[4 * line] = (short)((g_t16[4][0] * eeo[0] + g_t16[4][1] * eeo[1] + add) >> shift);
    dst[12 * line] = (short)((g_t16[12][0] * eeo[0] + g_t16[12][1] * eeo[1] + add) >> shift);

    for (k = 2; k < 16; k += 4) {
      dst[k*line] = (short)((g_t16[k][0] * eo[0] + g_t16[k][1] * eo[1] + g_t16[k][2] * eo[2] + g_t16[k][3] * eo[3] + add) >> shift);
    }

    for (k = 1; k < 16; k += 2) {
      dst[k*line] = (short)((g_t16[k][0] * o[0] + g_t16[k][1] * o[1] + g_t16[k][2] * o[2] + g_t16[k][3] * o[3] +
        g_t16[k][4] * o[4] + g_t16[k][5] * o[5] + g_t16[k][6] * o[6] + g_t16[k][7] * o[7] + add) >> shift);
    }

    src += 16;
    dst++;
  }
}


static void partial_butterfly_inverse_16_generic(int16_t *src, int16_t *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[8], o[8];
  int32_t ee[4], eo[4];
  int32_t eee[2], eeo[2];
  int32_t add = 1 << (shift - 1);

  for (j = 0; j < line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 8; k++)  {
      o[k] = g_t16[1][k] * src[line] + g_t16[3][k] * src[3 * line] + g_t16[5][k] * src[5 * line] + g_t16[7][k] * src[7 * line] +
        g_t16[9][k] * src[9 * line] + g_t16[11][k] * src[11 * line] + g_t16[13][k] * src[13 * line] + g_t16[15][k] * src[15 * line];
    }
    for (k = 0; k < 4; k++) {
      eo[k] = g_t16[2][k] * src[2 * line] + g_t16[6][k] * src[6 * line] + g_t16[10][k] * src[10 * line] + g_t16[14][k] * src[14 * line];
    }
    eeo[0] = g_t16[4][0] * src[4 * line] + g_t16[12][0] * src[12 * line];
    eee[0] = g_t16[0][0] * src[0] + g_t16[8][0] * src[8 * line];
    eeo[1] = g_t16[4][1] * src[4 * line] + g_t16[12][1] * src[12 * line];
    eee[1] = g_t16[0][1] * src[0] + g_t16[8][1] * src[8 * line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    for (k = 0; k < 2; k++) {
      ee[k] = eee[k] + eeo[k];
      ee[k + 2] = eee[1 - k] - eeo[1 - k];
    }
    for (k = 0; k < 4; k++) {
      e[k] = ee[k] + eo[k];
      e[k + 4] = ee[3 - k] - eo[3 - k];
    }
    for (k = 0; k < 8; k++) {
      dst[k] = (short)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k + 8] = (short)MAX(-32768, MIN(32767, (e[7 - k] - o[7 - k] + add) >> shift));
    }
    src++;
    dst += 16;
  }
}


static void partial_butterfly_32_generic(short *src, short *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[16], o[16];
  int32_t ee[8], eo[8];
  int32_t eee[4], eeo[4];
  int32_t eeee[2], eeeo[2];
  int32_t add = 1 << (shift - 1);

  for (j = 0; j < line; j++) {
    // E and O
    for (k = 0; k < 16; k++) {
      e[k] = src[k] + src[31 - k];
      o[k] = src[k] - src[31 - k];
    }
    // EE and EO
    for (k = 0; k < 8; k++) {
      ee[k] = e[k] + e[15 - k];
      eo[k] = e[k] - e[15 - k];
    }
    // EEE and EEO
    for (k = 0; k < 4; k++) {
      eee[k] = ee[k] + ee[7 - k];
      eeo[k] = ee[k] - ee[7 - k];
    }
    // EEEE and EEEO
    eeee[0] = eee[0] + eee[3];
    eeeo[0] = eee[0] - eee[3];
    eeee[1] = eee[1] + eee[2];
    eeeo[1] = eee[1] - eee[2];

    dst[0] = (short)((g_t32[0][0] * eeee[0] + g_t32[0][1] * eeee[1] + add) >> shift);
    dst[16 * line] = (short)((g_t32[16][0] * eeee[0] + g_t32[16][1] * eeee[1] + add) >> shift);
    dst[8 * line] = (short)((g_t32[8][0] * eeeo[0] + g_t32[8][1] * eeeo[1] + add) >> shift);
    dst[24 * line] = (short)((g_t32[24][0] * eeeo[0] + g_t32[24][1] * eeeo[1] + add) >> shift);
    for (k = 4; k < 32; k += 8) {
      dst[k*line] = (short)((g_t32[k][0] * eeo[0] + g_t32[k][1] * eeo[1] + g_t32[k][2] * eeo[2] + g_t32[k][3] * eeo[3] + add) >> shift);
    }
    for (k = 2; k < 32; k += 4) {
      dst[k*line] = (short)((g_t32[k][0] * eo[0] + g_t32[k][1] * eo[1] + g_t32[k][2] * eo[2] + g_t32[k][3] * eo[3] +
        g_t32[k][4] * eo[4] + g_t32[k][5] * eo[5] + g_t32[k][6] * eo[6] + g_t32[k][7] * eo[7] + add) >> shift);
    }
    for (k = 1; k < 32; k += 2) {
      dst[k*line] = (short)((g_t32[k][0] * o[0] + g_t32[k][1] * o[1] + g_t32[k][2] * o[2] + g_t32[k][3] * o[3] +
        g_t32[k][4] * o[4] + g_t32[k][5] * o[5] + g_t32[k][6] * o[6] + g_t32[k][7] * o[7] +
        g_t32[k][8] * o[8] + g_t32[k][9] * o[9] + g_t32[k][10] * o[10] + g_t32[k][11] * o[11] +
        g_t32[k][12] * o[12] + g_t32[k][13] * o[13] + g_t32[k][14] * o[14] + g_t32[k][15] * o[15] + add) >> shift);
    }
    src += 32;
    dst++;
  }
}


static void partial_butterfly_inverse_32_generic(int16_t *src, int16_t *dst,
  int32_t shift, int32_t line)
{
  int32_t j, k;
  int32_t e[16], o[16];
  int32_t ee[8], eo[8];
  int32_t eee[4], eeo[4];
  int32_t eeee[2], eeeo[2];
  int32_t add = 1 << (shift - 1);

  for (j = 0; j<line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 16; k++) {
      o[k] = g_t32[1][k] * src[line] + g_t32[3][k] * src[3 * line] + g_t32[5][k] * src[5 * line] + g_t32[7][k] * src[7 * line] +
        g_t32[9][k] * src[9 * line] + g_t32[11][k] * src[11 * line] + g_t32[13][k] * src[13 * line] + g_t32[15][k] * src[15 * line] +
        g_t32[17][k] * src[17 * line] + g_t32[19][k] * src[19 * line] + g_t32[21][k] * src[21 * line] + g_t32[23][k] * src[23 * line] +
        g_t32[25][k] * src[25 * line] + g_t32[27][k] * src[27 * line] + g_t32[29][k] * src[29 * line] + g_t32[31][k] * src[31 * line];
    }
    for (k = 0; k < 8; k++) {
      eo[k] = g_t32[2][k] * src[2 * line] + g_t32[6][k] * src[6 * line] + g_t32[10][k] * src[10 * line] + g_t32[14][k] * src[14 * line] +
        g_t32[18][k] * src[18 * line] + g_t32[22][k] * src[22 * line] + g_t32[26][k] * src[26 * line] + g_t32[30][k] * src[30 * line];
    }
    for (k = 0; k < 4; k++) {
      eeo[k] = g_t32[4][k] * src[4 * line] + g_t32[12][k] * src[12 * line] + g_t32[20][k] * src[20 * line] + g_t32[28][k] * src[28 * line];
    }
    eeeo[0] = g_t32[8][0] * src[8 * line] + g_t32[24][0] * src[24 * line];
    eeeo[1] = g_t32[8][1] * src[8 * line] + g_t32[24][1] * src[24 * line];
    eeee[0] = g_t32[0][0] * src[0] + g_t32[16][0] * src[16 * line];
    eeee[1] = g_t32[0][1] * src[0] + g_t32[16][1] * src[16 * line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    eee[0] = eeee[0] + eeeo[0];
    eee[3] = eeee[0] - eeeo[0];
    eee[1] = eeee[1] + eeeo[1];
    eee[2] = eeee[1] - eeeo[1];
    for (k = 0; k < 4; k++) {
      ee[k] = eee[k] + eeo[k];
      ee[k + 4] = eee[3 - k] - eeo[3 - k];
    }
    for (k = 0; k < 8; k++) {
      e[k] = ee[k] + eo[k];
      e[k + 8] = ee[7 - k] - eo[7 - k];
    }
    for (k = 0; k<16; k++) {
      dst[k] = (short)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k + 16] = (short)MAX(-32768, MIN(32767, (e[15 - k] - o[15 - k] + add) >> shift));
    }
    src++;
    dst += 32;
  }
}


int strategy_register_partial_butterfly_generic(void* opaque)
{
  bool success = true;

  success &= strategyselector_register(opaque, "partial_butterfly_4", "generic", 0, &partial_butterfly_4_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_8", "generic", 0, &partial_butterfly_8_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_16", "generic", 0, &partial_butterfly_16_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_32", "generic", 0, &partial_butterfly_32_generic);

  success &= strategyselector_register(opaque, "partial_butterfly_inverse_4", "generic", 0, &partial_butterfly_inverse_4_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_inverse_8", "generic", 0, &partial_butterfly_inverse_8_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_inverse_16", "generic", 0, &partial_butterfly_inverse_16_generic);
  success &= strategyselector_register(opaque, "partial_butterfly_inverse_32", "generic", 0, &partial_butterfly_inverse_32_generic);
  return success;
}
