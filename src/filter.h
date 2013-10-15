#ifndef FILTER_H_
#define FILTER_H_
/**
 * \file
 * \brief Filtering, such as deblocking.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "encoder.h"


//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
// Deblocking
void filter_deblock_cu(encoder_control *encoder, int32_t x_cu, int32_t y_cu,
                       int8_t depth, int32_t edge);
void filter_deblock_edge_luma(encoder_control *encoder, 
                              int32_t x_pos, int32_t y_pos, 
                              int8_t depth, int8_t dir);
void filter_deblock_edge_chroma(encoder_control *encoder, 
                                int32_t xpos, int32_t ypos, 
                                int8_t depth, int8_t dir);
void filter_deblock(encoder_control *encoder);
void filter_deblock_luma(pixel *src, int32_t offset, int32_t tc , int8_t sw,
                         int8_t part_p_nofilter, int8_t part_q_nofilter,
                         int32_t thr_cut, 
                         int8_t filter_second_p, int8_t filter_second_q);
void filter_deblock_chroma(pixel *src, int32_t offset, int32_t tc,
                           int8_t part_p_nofilter, int8_t part_q_nofilter);

// INTERPOLATION
void filter_inter_halfpel_chroma(int16_t *src, int16_t src_stride, int width, int height,
                                 int16_t *dst, int16_t dst_stride,  int8_t hor_flag, int8_t ver_flag);

// SAO

//////////////////////////////////////////////////////////////////////////
// MACROS
#define EDGE_VER 0
#define EDGE_HOR 1

#endif
