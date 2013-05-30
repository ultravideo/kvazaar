/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file global.h
    \brief Contains global includes
    \author Marko Viitanen
    \date 2013-05
  
    This file should be included in every C-file.
*/
#ifndef __GLOBAL_H
#define __GLOBAL_H

/* CONFIG VARIABLES */
#define LCU_WIDTH 64 /*!< Largest Coding Unit */

#define MAX_SEARCH_DEPTH 3
#define MIN_SEARCH_DEPTH 2

#define MAX_DEPTH 3  /*!< smallest CU is LCU_WIDTH>>MAX_DEPTH */
#define MIN_SIZE 3   /*!< log2_min_coding_block_size */

#define ENABLE_PCM 0 /*!< Setting to 1 will enable using PCM blocks (current intra-search does not consider PCM) */
#define ENABLE_SIGN_HIDING 0 /*!< NEED QUANT CHANGES! */

/* END OF CONFIG VARIABLES */


//Including stdint.h, 
#ifdef _MSC_VER
  #include "../include/stdint.h"
#else
  #include <stdint.h>
#endif

/* Some tools */
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define CLIP(low,high,value) MAX((low),MIN((high),(value)))
#define SWAP(a,b,swaptype) { swaptype tempval; tempval = a; a = b; b = tempval; }

#define VERSION_STRING "0.2               "
#define VERSION 0.2

//#define VERBOSE 1


#define SIZE_2Nx2N 0
#define SIZE_2NxN  1
#define SIZE_Nx2N  2
#define SIZE_NxN   3
#define SIZE_NONE  15

/*
#define MODE_SKIP  0
#define MODE_INTER 1
#define MODE_INTRA 2
#define MODE_NONE  15
*/


/* Inlining functions */
#ifdef _MSC_VER /* Visual studio */
  #define INLINE __forceinline
  #pragma inline_recursion(on)
#else /* others */
  #define INLINE inline
#endif

#define free_pointer(pointer) { free(pointer); pointer = NULL; }

#endif