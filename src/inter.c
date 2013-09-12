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
void inter_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, CU_info* cur_cu)
{
  uint32_t x,y,d;
  /* Width in smallest CU */
  int width_in_SCU = pic->width_in_LCU<<MAX_DEPTH;
  int block_SCU_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  for(y = yCtb; y < yCtb+block_SCU_width; y++)
  {
    int CUpos = y*width_in_SCU;
    for(x = xCtb; x < xCtb+block_SCU_width; x++)
    {
      for(d = 0; d < MAX_DEPTH+1; d++)
      {
        pic->CU[d][CUpos+x].depth = depth;
        pic->CU[d][CUpos+x].type  = CU_INTER;
        pic->CU[d][CUpos+x].inter.mode = cur_cu->inter.mode;
        pic->CU[d][CUpos+x].inter.mv[0] = cur_cu->inter.mv[0];
        pic->CU[d][CUpos+x].inter.mv[1] = cur_cu->inter.mv[1];
        pic->CU[d][CUpos+x].inter.mv_dir = cur_cu->inter.mv_dir;
      }
    }
  }
}

/*!
 \brief Reconstruct inter block
 \param ref picture to copy the data from
 \param xpos block x position
 \param ypos block y position
 \param width block width
 \param mv[2] motion vector
 \param dst destination picture
 \returns Void
*/
void inter_recon(picture* ref,int32_t xpos, int32_t ypos,int32_t width, int16_t mv[2], picture* dst)
{
  int x,y,coord_x,coord_y;

  /* negative overflow present */
  int8_t overflow_neg_x = (xpos+mv[0] < 0)?1:0;
  int8_t overflow_neg_y = (ypos+mv[1] < 0)?1:0;

  /* positive overflow present */
  int8_t overflow_pos_x = (xpos+mv[0]+width > ref->width )?1:0;
  int8_t overflow_pos_y = (ypos+mv[1]+width > ref->height)?1:0;

  /* TODO: Fractional pixel support */
  mv[0] = mv[0]>>2;
  mv[1] = mv[1]>>2;

  /* With overflow present, more checking */
  if(overflow_neg_x || overflow_neg_y || overflow_pos_x || overflow_pos_y)
  {
    /* Copy Luma with boundary checking */
    for(y = ypos; y < ypos+width; y++)
    {
      for(x = xpos; x < xpos+width; x++)
      {
        coord_x = x;
        coord_y = y;
        overflow_neg_x = (x < 0)?1:0;
        overflow_neg_y = (y < 0)?1:0;

        overflow_pos_x = (x >= ref->width )?1:0;
        overflow_pos_y = (y >= ref->height)?1:0;

        if(overflow_neg_x) {
          coord_x = 0;
        } else if(overflow_pos_x) {
          coord_x = ref->width-1;
        }

        if(overflow_neg_y) {
          coord_y = 0;
        } else if(overflow_pos_y) {
          coord_y = ref->height-1;
        }

        dst->yRecData[y*dst->width+x] = ref->yRecData[(coord_y+mv[1])*ref->width+(coord_x+mv[0])];
      }
    }

    /* Copy Chroma with boundary checking */
    for(y = ypos>>1; y < (ypos+width)>>1; y++)
    {
      for(x = xpos>>1; x < (xpos+width)>>1; x++)
      {
        coord_x = x;
        coord_y = y;
        overflow_neg_x = (x < 0)?1:0;
        overflow_neg_y = (y < 0)?1:0;

        overflow_pos_x = (x >= ref->width>>1 )?1:0;
        overflow_pos_y = (y >= ref->height>>1)?1:0;

        if(overflow_neg_x) {
          coord_x = 0;
        } else if(overflow_pos_x) {
          coord_x = (ref->width>>1)-1;
        }

        if(overflow_neg_y) {
          coord_y = 0;
        } else if(overflow_pos_y) {
          coord_y = (ref->height>>1)-1;
        }

        dst->uRecData[y*(dst->width>>1)+x] = ref->uRecData[(coord_y+(mv[1]>>1))*ref->width+(coord_x+(mv[0]>>1))];
        dst->vRecData[y*(dst->width>>1)+x] = ref->vRecData[(coord_y+(mv[1]>>1))*ref->width+(coord_x+(mv[0]>>1))];
      }
    }
  }
  else
  {
    /* Copy Luma */
    for(y = ypos; y < ypos+width; y++)
    {
      for(x = xpos; x < xpos+width; x++)
      {
        dst->yRecData[y*dst->width+x] = ref->yRecData[(y+mv[1])*ref->width+x+mv[0]];
      }
    }

    /* Copy Chroma */
    for(y = ypos>>1; y < (ypos+width)>>1; y++)
    {
      for(x = xpos>>1; x < (xpos+width)>>1; x++)
      {
        dst->uRecData[y*(dst->width>>1)+x] = ref->uRecData[(y+(mv[1]>>1))*(ref->width>>1)+x+(mv[0]>>1)];
        dst->vRecData[y*(dst->width>>1)+x] = ref->vRecData[(y+(mv[1]>>1))*(ref->width>>1)+x+(mv[0]>>1)];
      }
    }
  }

}