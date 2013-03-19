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

int8_t intra_getDirLumaPredictor(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t* preds);
void intra_DCPredFiltering(uint8_t* pSrc, int32_t iSrcStride, uint8_t* rpDst, int32_t iDstStride, int32_t iWidth, int32_t iHeight );

void intra_buildReferenceBorder(picture* pic, int32_t xCtb, int32_t yCtb,int8_t outwidth, int16_t* dst, int32_t dststride, int8_t chroma);

/* Predictions */
int16_t intra_prediction(uint8_t* orig,uint32_t origstride,int16_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, uint32_t *sad);

int16_t intra_getDCPred(int16_t* pic, uint16_t picwidth,uint32_t xpos, uint32_t ypos, uint8_t width);
void intra_getPlanarPred(int16_t* src,int32_t srcstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride);
void intra_getAngularPred(int16_t* pSrc, int32_t srcStride, int16_t* rpDst, int32_t dstStride, int32_t width, int32_t height, int32_t dirMode, int8_t leftAvail,int8_t topAvail, int8_t filter);

void intra_recon(int16_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, int8_t mode, int8_t chroma);


#endif
