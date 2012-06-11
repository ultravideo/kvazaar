/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file global.h
    \brief Contains global includes
    \author Marko Viitanen
    \date 2012-05
  
    This file should be included to every C-file.
*/
#ifndef __GLOBAL_H
#define __GLOBAL_H

/* CONFIG VARIABLES */
#define LCU_WIDTH 64 /*!< Largest Coding Unit */
#define MAX_DEPTH 2

#define ENABLE_PCM 1

/* END OF CONFIG VARIABLES */


//Including stdint.h, 
#ifdef _MSC_VER
  #include "../include/stdint.h"
#else
  #include <stdint.h>
#endif

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define VERSION_STRING "0.1               "
#define VERSION 0.1

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

#endif