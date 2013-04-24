/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file inter.c
    \brief Inter functions
    \author Marko Viitanen
    \date 2013-04
    
    Inter functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "picture.h"
#include "inter.h"


/*!
 \brief Set block mode (and init typedata)
 \param pic picture to use
 \param xCtb x CU position (smallest CU)
 \param yCtb y CU position (smallest CU)
 \param depth current CU depth
 \param mode mode to set
 \returns Void
*/
void inter_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, uint8_t mode)
{
  uint32_t x,y,d;
  /* Width in smallest CU */
  int width_in_SCU = pic->width/(LCU_WIDTH>>MAX_DEPTH);
  int block_SCU_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  for(y = yCtb; y < yCtb+block_SCU_width; y++)
  {
    int CUpos = y*width_in_SCU;
    for(x = xCtb; x < xCtb+block_SCU_width; x++)
    {      
      for(d = 0; d < MAX_DEPTH+1; d++)
      {
        pic->CU[d][CUpos+x].depth = depth;
        pic->CU[d][CUpos+x].type = CU_INTER;
        pic->CU[d][CUpos+x].inter.mode = mode;        
      }
    }
  }
}

