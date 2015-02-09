#ifndef GLOBAL_H_
#define GLOBAL_H_
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
 * \brief Header that is included in every other header.
 *
 * This file contains global constants that can be referred to from any header
 * or source file. It also contains some helper macros and includes stdint.h
 * so that any file can refer to integer types with exact widths.
 */

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>

#if defined(_MSC_VER) && defined(_M_AMD64)
  #define X86_64
#endif

#if defined(__GNUC__) && defined(__x86_64__)
  #define X86_64
#endif

#define BIT_DEPTH 8
#define PIXEL_MIN 0
#define PIXEL_MAX ((1 << BIT_DEPTH) - 1)

#if BIT_DEPTH == 8
typedef uint8_t pixel;
#else
typedef uint16_t pixel;
#endif
typedef int16_t coefficient;

//#define VERBOSE 1

/* CONFIG VARIABLES */

//spec: references to variables defined in Rec. ITU-T H.265 (04/2013)

// Limits for prediction block sizes. 0 = 64x64, 4 = 4x4.
#define PU_DEPTH_INTER_MIN 0
#define PU_DEPTH_INTER_MAX 3
#define PU_DEPTH_INTRA_MIN 0
#define PU_DEPTH_INTRA_MAX 4

// Maximum CU depth when descending form LCU level.
#define MAX_DEPTH 3  /*!< spec: log2_diff_max_min_luma_coding_block_size */
// Minimum log2 size of CUs.
#define MIN_SIZE 3   /*!< spec: MinCbLog2SizeY */
// Minimum log2 size of PUs.
#define MAX_PU_DEPTH 4 /*!< Search is started at depth 0 and goes in Z-order to MAX_PU_DEPTH, see search_cu() */

// Minimum log2 transform sizes.
#define TR_DEPTH_INTER 2 /*!< spec: max_transform_hierarchy_depth_inter */

#define ENABLE_PCM 0 /*!< spec: pcm_enabled_flag, Setting to 1 will enable using PCM blocks (current intra-search does not consider PCM) */

#define ENABLE_TEMPORAL_MVP 0 /*!< Enable usage of temporal Motion Vector Prediction */

#define OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD 0 /*!< skip residual coding when it's under _some_ threshold */

/* END OF CONFIG VARIABLES */

#define CU_MIN_SIZE_PIXELS (1 << MIN_SIZE) /*!< pow(2, MIN_SIZE) */
#define LCU_WIDTH (1 << (MIN_SIZE + MAX_DEPTH)) /*!< spec: CtbSizeY */
#define LCU_WIDTH_C (LCU_WIDTH / 2) /*!< spec: CtbWidthC and CtbHeightC */

#define TR_MAX_LOG2_SIZE 5 /*!< spec: Log2MaxTrafoSize <= Min(CtbLog2SizeY, 5) */
#define TR_MAX_WIDTH (1 << TR_MAX_LOG2_SIZE)
#define TR_MIN_LOG2_SIZE 2 /*!< spec: Log2MinTrafoSize */
#define TR_MIN_WIDTH (1 << TR_MIN_LOG2_SIZE)

#if LCU_WIDTH != 64
  #error "Kvazaar only support LCU_WIDTH == 64"
#endif

#define LCU_LUMA_SIZE (LCU_WIDTH * LCU_WIDTH)
#define LCU_CHROMA_SIZE (LCU_WIDTH * LCU_WIDTH >> 2)

#define MAX_REF_PIC_COUNT 16
#define DEFAULT_REF_PIC_COUNT 3

#define AMVP_MAX_NUM_CANDS 2
#define AMVP_MAX_NUM_CANDS_MEM 3
#define MRG_MAX_NUM_CANDS 5

/* Some tools */
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define CLIP(low,high,value) MAX((low),MIN((high),(value)))
#define SWAP(a,b,swaptype) { swaptype tempval; tempval = a; a = b; b = tempval; }
#define CU_WIDTH_FROM_DEPTH(depth) (LCU_WIDTH >> depth)
#define WITHIN(val, min_val, max_val) ((min_val) <= (val) && (val) <= (max_val))
#define PU_INDEX(x_pu, y_pu) (((x_pu) % 2)  + 2 * ((y_pu) % 2))
#define CEILDIV(x,y) (((x) + (y) - 1) / (y))

#define LOG2_LCU_WIDTH 6
// CU_TO_PIXEL = y * lcu_width * pic_width + x * lcu_width
#define CU_TO_PIXEL(x, y, depth, stride) (((y) << (LOG2_LCU_WIDTH - (depth))) * (stride) \
                                         + ((x) << (LOG2_LCU_WIDTH - (depth))))
//#define SIGN3(x) ((x) > 0) ? +1 : ((x) == 0 ? 0 : -1)
#define SIGN3(x) (((x) > 0) - ((x) < 0))

#define VERSION_STRING "0.4.0"

//#define VERBOSE 1

#define SAO_ABS_OFFSET_MAX ((1 << (MIN(BIT_DEPTH, 10) - 5)) - 1)


#define SIZE_2Nx2N 0
#define SIZE_2NxN  1
#define SIZE_Nx2N  2
#define SIZE_NxN   3
#define SIZE_NONE  15

#define MAX_TILES_PER_DIM 48
#define MAX_SLICES 16

/* Inlining functions */
#ifdef _MSC_VER /* Visual studio */
  #define INLINE __forceinline
  #pragma inline_recursion(on)
#else /* others */
  #define INLINE inline
#endif

// Return the next aligned address for *p. Result is at most alignment larger than p.
#define ALIGNED_POINTER(p, alignment) (void*)((intptr_t)(p) + (alignment) - ((intptr_t)(p) % (alignment)))
// 32 bytes is enough for AVX2
#define SIMD_ALIGNMENT 32

#ifdef _MSC_VER
// Buggy VS2010 throws intellisense warnings if void* is not casted.
  #define MALLOC(type, num) (type *)malloc(sizeof(type) * (num))
#else
  #define MALLOC(type, num) malloc(sizeof(type) * (num))
#endif

#define FREE_POINTER(pointer) { free((void*)pointer); pointer = NULL; }
#define MOVE_POINTER(dst_pointer,src_pointer) { dst_pointer = src_pointer; src_pointer = NULL; }

#ifndef MAX_INT
#define MAX_INT 0x7FFFFFFF
#endif
#ifndef MAX_INT64
#define MAX_INT64 0x7FFFFFFFFFFFFFFFLL
#endif
#ifndef MAX_DOUBLE
#define MAX_DOUBLE 1.7e+308
#endif

//For transform.h and encoder.h
#define SCALING_LIST_4x4      0
#define SCALING_LIST_8x8      1
#define SCALING_LIST_16x16    2
#define SCALING_LIST_32x32    3
#define SCALING_LIST_SIZE_NUM 4
#define SCALING_LIST_NUM      6
#define MAX_MATRIX_COEF_NUM   64
#define SCALING_LIST_REM_NUM  6

#define MAX_TR_DYNAMIC_RANGE 15


//DEBUG BITMASK
#define _DEBUG_PERF_FRAME_LEVEL 0x0001
#define _DEBUG_PERF_JOB 0x0002
#define _DEBUG_PERF_ENCODE_LCU 0x0004
#define _DEBUG_PERF_SAO_RECONSTRUCT_FRAME 0x0008
#define _DEBUG_PERF_WRITE_BITSTREAM_LEAF 0x0010
#define _DEBUG_PERF_SEARCH_CU 0x0020
#define _DEBUG_PERF_SEARCH_PIXELS 0x0040

//Constants
typedef enum { COLOR_Y = 0, COLOR_U, COLOR_V, NUM_COLORS } color_index;
enum { SLICE_B = 0, SLICE_P = 1, SLICE_I = 2 };

#endif
