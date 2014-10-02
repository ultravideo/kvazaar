#ifndef STRATEGIES_DCT_GENERIC_H_
#define STRATEGIES_DCT_GENERIC_H_
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

extern const int16_t g_dst[4][4];
extern const int16_t g_t4[4][4];
extern const int16_t g_t8[8][8];
extern const int16_t g_t16[16][16];
extern const int16_t g_t32[32][32];

extern const int16_t g_dst_t[4][4];
extern const int16_t g_t4_t[4][4];
extern const int16_t g_t8_t[8][8];
extern const int16_t g_t16_t[16][16];
extern const int16_t g_t32_t[32][32];

int strategy_register_dct_generic(void* opaque);

#endif //STRATEGIES_DCT_GENERIC_H_
