#ifndef STRATEGIES_PARTIAL_BUTTERFLY_H_
#define STRATEGIES_PARTIAL_BUTTERFLY_H_
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

typedef unsigned (partial_butterfly_func)(short *src, short *dst, int shift, int line);


// Declare function pointers.
extern partial_butterfly_func * partial_butterfly_4;
extern partial_butterfly_func * partial_butterfly_8;
extern partial_butterfly_func * partial_butterfly_16;
extern partial_butterfly_func * partial_butterfly_32;

extern partial_butterfly_func * partial_butterfly_inverse_4;
extern partial_butterfly_func * partial_butterfly_inverse_8;
extern partial_butterfly_func * partial_butterfly_inverse_16;
extern partial_butterfly_func * partial_butterfly_inverse_32;


int strategy_register_partial_butterfly(void* opaque);
partial_butterfly_func * get_partial_butterfly_func(unsigned n);
partial_butterfly_func * get_partial_butterfly_inverse_func(unsigned n);


#define STRATEGIES_PARTIAL_BUTTERFLY_EXPORTS \
  {"partial_butterfly_4", (void**) &partial_butterfly_4}, \
  {"partial_butterfly_8", (void**) &partial_butterfly_8}, \
  {"partial_butterfly_16", (void**) &partial_butterfly_16}, \
  {"partial_butterfly_32", (void**) &partial_butterfly_32}, \
  \
  {"partial_butterfly_inverse_4", (void**)&partial_butterfly_inverse_4}, \
  {"partial_butterfly_inverse_8", (void**)&partial_butterfly_inverse_8}, \
  {"partial_butterfly_inverse_16", (void**)&partial_butterfly_inverse_16}, \
  {"partial_butterfly_inverse_32", (void**)&partial_butterfly_inverse_32}, \



#endif //STRATEGIES_PARTIAL_BUTTERFLY_H_
