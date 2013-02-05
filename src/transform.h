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

extern int32_t* g_quant_coeff[4][6][6][3];

void quant(encoder_control* encoder, int16_t* pSrc, int16_t* pDes, int32_t iWidth,
           int32_t iHeight, uint32_t *uiAcSum, int8_t eTType);

void transform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int8_t uiMode);
void scalinglist_init();
void scalinglist_destroy();

#endif
