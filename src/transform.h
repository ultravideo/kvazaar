/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file transform.h
    \brief Transform functions
    \author Marko Viitanen
    \date 2012-09
    
    Transform functions
*/
#ifndef __TRANSFORM_H
#define __TRANSFORM_H

extern int32_t* g_quant_coeff[4][6][6];
extern const int32_t g_quantIntraDefault8x8[64];

void quant(encoder_control* encoder, int16_t* pSrc, int16_t* pDes,int32_t iWidth,int32_t iHeight, int8_t eTType );
void dequant(encoder_control* encoder, int16_t* piQCoef, int16_t* piCoef, int32_t iWidth, int32_t iHeight );

void transform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int8_t uiMode);
void itransform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int8_t uiMode);

void scalinglist_init();
void scalinglist_processEnc( int32_t *coeff, int32_t *quantcoeff, int32_t quantScales, uint32_t height,uint32_t width, uint32_t ratio, int32_t sizuNum, uint32_t dc, uint8_t flat);
void scalinglist_process();
void scalinglist_set(int32_t *coeff, uint32_t listId, uint32_t sizeId, uint32_t qp);
void scalinglist_destroy();

#endif
