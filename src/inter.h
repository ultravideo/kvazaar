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

void inter_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, CU_info* cur_cu);
void inter_recon(picture* ref,int32_t xpos, int32_t ypos,int32_t width, int16_t mv[2], picture* dst);
#endif
