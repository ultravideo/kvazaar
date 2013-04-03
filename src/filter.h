/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems 2012.
 */

/*! \file filter.h
    \brief filter
    \author Marko Viitanen
    \date 2013-04
    
    Filtering function headers
*/
#ifndef __FILTER_H
#define __FILTER_H

#define EDGE_VER 0
#define EDGE_HOR 1


void filter_deblock_edge_luma(encoder_control* encoder, int32_t xpos, int32_t ypos, int8_t depth, int32_t edge, int8_t chroma, int8_t dir);

#endif
