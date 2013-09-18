/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing 2013.
 */

/*! \file filter.h
    \brief filter
    \author Marko Viitanen
    \date 2013-05
    
    Filtering function headers
*/
#ifndef __FILTER_H
#define __FILTER_H

#include "global.h"

#include "encoder.h"


#define EDGE_VER 0
#define EDGE_HOR 1

/* DEBLOCKING */
void filter_deblock_CU(encoder_control* encoder, int32_t xCtb, int32_t yCtb, int8_t depth, int32_t edge);
void filter_deblock_edge_luma(encoder_control* encoder, int32_t xpos, int32_t ypos, int8_t depth, int8_t dir);
void filter_deblock_edge_chroma(encoder_control* encoder,int32_t xpos, int32_t ypos, int8_t depth, int8_t dir);
void filter_deblock(encoder_control* encoder);
void filter_deblock_luma( uint8_t* src, int32_t offset, int32_t tc , int8_t sw, int8_t part_P_nofilter, int8_t part_Q_nofilter, int32_t thr_cut, int8_t filter_second_P, int8_t filter_second_Q);
void filter_deblock_chroma( uint8_t* src, int32_t offset, int32_t tc ,int8_t part_P_nofilter, int8_t part_Q_nofilter);

/* SAO */

#endif
