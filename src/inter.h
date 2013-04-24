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

void inter_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, uint8_t mode);

#endif
