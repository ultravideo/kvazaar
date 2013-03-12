/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file intra.h
    \brief Intra function headers
    \author Marko Viitanen
    \date 2013-03
    
    Intra functions
*/
#ifndef __INTRA_H
#define __INTRA_H

void intra_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, uint8_t mode);
int8_t intra_getBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth);
int16_t intra_getDCPred(uint8_t* pic, uint16_t picwidth,uint32_t xpos, uint32_t ypos, uint8_t width);
int8_t intra_getDirLumaPredictor(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t* preds);
void intra_DCPredFiltering(uint8_t* pSrc, int32_t iSrcStride, uint8_t* rpDst, int32_t iDstStride, int32_t iWidth, int32_t iHeight );
void intra_getPlanarPred(uint8_t* src,uint32_t srcstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride);

#endif
