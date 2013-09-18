/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file inter.h
    \brief Inter function headers
    \author Marko Viitanen
    \date 2013-04
    
    Inter functions
*/
#ifndef __INTER_H
#define __INTER_H

#include "global.h"

#include "picture.h"
#include "encoder.h"


void inter_setBlockMode(picture* pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, CU_info* cur_cu);
void inter_recon(picture *ref,int32_t xpos, int32_t ypos,int32_t width, int16_t mv[2], picture* dst);

void inter_get_mv_cand(encoder_control *encoder, int32_t x_cu, int32_t y_cu, int8_t depth, int16_t mv_cand[2][2]);

#endif
