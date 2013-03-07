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
int16_t intra_getDCPred(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth);
int8_t intra_getDirLumaPredictor(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t* preds);

#endif
