/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "inter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"


/*!
 \brief Set block mode (and init typedata)
 \param pic picture to use
 \param xCtb x CU position (smallest CU)
 \param yCtb y CU position (smallest CU)
 \param depth current CU depth
 \param mode mode to set
 \returns Void
*/
void inter_setBlockMode(picture* pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, CU_info* cur_cu)
{
  uint32_t x,y,d;
  /* Width in smallest CU */
  int width_in_SCU = pic->width_in_lcu<<MAX_DEPTH;
  int block_SCU_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  for(y = y_cu; y < y_cu+block_SCU_width; y++)
  {
    int CUpos = y*width_in_SCU;
    for(x = x_cu; x < x_cu+block_SCU_width; x++)
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
  if (overflow_neg_x || overflow_neg_y || overflow_pos_x || overflow_pos_y) {
    /* Copy Luma with boundary checking */
    for (y = ypos; y < ypos+width; y++) {
      for (x = xpos; x < xpos+width; x++) {
        coord_x = x;
        coord_y = y;
        overflow_neg_x = (x+mv[0] < 0)?1:0;
        overflow_neg_y = (y+mv[1] < 0)?1:0;

        overflow_pos_x = (x+mv[0] >= ref->width )?1:0;
        overflow_pos_y = (y+mv[1] >= ref->height)?1:0;

        if(overflow_neg_x) {
          coord_x = -mv[0];
        } else if(overflow_pos_x) {
          coord_x = ref->width-1-mv[0];
        }

        if(overflow_neg_y) {
          coord_y = -mv[1];
        } else if(overflow_pos_y) {
          coord_y = ref->height-1-mv[1];
        }

        dst->y_recdata[y*dst->width+x] = ref->y_recdata[(coord_y+mv[1])*ref->width+(coord_x+mv[0])];
      }
    }

    /* Copy Chroma with boundary checking */
    for (y = ypos>>1; y < (ypos+width)>>1; y++) {
      for (x = xpos>>1; x < (xpos+width)>>1; x++) {
        coord_x = x;
        coord_y = y;
        overflow_neg_x = (x+(mv[0]>>1) < 0)?1:0;
        overflow_neg_y = (y+(mv[1]>>1) < 0)?1:0;

        overflow_pos_x = (x+(mv[0]>>1) >= ref->width>>1 )?1:0;
        overflow_pos_y = (y+(mv[1]>>1) >= ref->height>>1)?1:0;

        if(overflow_neg_x) {
          coord_x = -(mv[0]>>1);
        } else if(overflow_pos_x) {
          coord_x = ((ref->width-mv[0])>>1)-1;
        }

        if(overflow_neg_y) {
          coord_y = -(mv[1]>>1);
        } else if(overflow_pos_y) {
          coord_y = ((ref->height-mv[1])>>1)-1;
        }

        dst->u_recdata[y*(dst->width>>1)+x] = ref->u_recdata[(coord_y+(mv[1]>>1))*(ref->width>>1)+(coord_x+(mv[0]>>1))];
        dst->v_recdata[y*(dst->width>>1)+x] = ref->v_recdata[(coord_y+(mv[1]>>1))*(ref->width>>1)+(coord_x+(mv[0]>>1))];
      }
    }
  } else {
    /* Copy Luma */
    for (y = ypos; y < ypos+width; y++) {
      for (x = xpos; x < xpos+width; x++) {
        dst->y_recdata[y*dst->width+x] = ref->y_recdata[(y+mv[1])*ref->width+x+mv[0]];
      }
    }

    /* Copy Chroma */
    for (y = ypos>>1; y < (ypos+width)>>1; y++) {
      for (x = xpos>>1; x < (xpos+width)>>1; x++) {
        dst->u_recdata[y*(dst->width>>1)+x] = ref->u_recdata[(y+(mv[1]>>1))*(ref->width>>1)+x+(mv[0]>>1)];
        dst->v_recdata[y*(dst->width>>1)+x] = ref->v_recdata[(y+(mv[1]>>1))*(ref->width>>1)+x+(mv[0]>>1)];
      }
    }
  }
}

/*!
 \brief Get MV prediction for current block
 \param encoder encoder control struct to use
 \param xCtb block x position in SCU
 \param yCtb block y position in SCU
 \param depth current block depth
 \param mv_pred[2][2] 2x motion vector prediction
 \returns Void
*/
void inter_get_mv_cand(encoder_control *encoder,int32_t xCtb, int32_t yCtb,int8_t depth, int16_t mv_cand[2][2])
{
  uint8_t cur_block_in_scu = (LCU_WIDTH>>depth) / CU_MIN_SIZE_PIXELS;
  uint8_t candidates = 0;

  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  CU_info *b0, *b1, *b2, *a0, *a1;

  b0 = b1 = b2 = a0 = a1 = NULL;

  if (xCtb != 0) {    
    a1 = &encoder->in.cur_pic->CU[depth][xCtb-1+(yCtb+cur_block_in_scu-1)*(encoder->in.width_in_lcu<<MAX_DEPTH)];
    if(!a1->coded) a1 = NULL;

    if (yCtb+cur_block_in_scu < encoder->in.height_in_lcu<<MAX_DEPTH) {
      a0 = &encoder->in.cur_pic->CU[depth][xCtb-1+(yCtb+cur_block_in_scu)*(encoder->in.width_in_lcu<<MAX_DEPTH)];
      if(!a0->coded) a0 = NULL;
    }
  }

  if (yCtb != 0) {
    b0 = &encoder->in.cur_pic->CU[depth][xCtb+cur_block_in_scu+(yCtb-1)*(encoder->in.width_in_lcu<<MAX_DEPTH)];
    if (!b0->coded) b0 = NULL;
    b1 = &encoder->in.cur_pic->CU[depth][xCtb+cur_block_in_scu-1+(yCtb-1)*(encoder->in.width_in_lcu<<MAX_DEPTH)];
    if (!b1->coded) b1 = NULL;

    if (xCtb != 0) {
      b2 = &encoder->in.cur_pic->CU[depth][xCtb-1+(yCtb-1)*(encoder->in.width_in_lcu<<MAX_DEPTH)];
      if(!b2->coded) b2 = NULL;
    }
  }

  /* Left predictors */
  if (a0 && a0->type == CU_INTER) {
    mv_cand[candidates][0] = a0->inter.mv[0];
    mv_cand[candidates][1] = a0->inter.mv[1];
    candidates++;
  } else if (a1 && a1->type == CU_INTER) {
    mv_cand[candidates][0] = a1->inter.mv[0];
    mv_cand[candidates][1] = a1->inter.mv[1];
    candidates++;
  }

  /* Top predictors */
  if (b0 && b0->type == CU_INTER) {
    mv_cand[candidates][0] = b0->inter.mv[0];
    mv_cand[candidates][1] = b0->inter.mv[1];
    candidates++;
  } else if (b1 && b1->type == CU_INTER) {
    mv_cand[candidates][0] = b1->inter.mv[0];
    mv_cand[candidates][1] = b1->inter.mv[1];
    candidates++;
  } else if(b2 && b2->type == CU_INTER) {
    mv_cand[candidates][0] = b2->inter.mv[0];
    mv_cand[candidates][1] = b2->inter.mv[1];
    candidates++;
  }

  /* Remove identical candidate */
  if(candidates == 2 && mv_cand[0][0] == mv_cand[1][0] && mv_cand[0][1] == mv_cand[1][1]) {
    candidates = 1;
  }

#if ENABLE_TEMPORAL_MVP
  if(candidates < 2) {
    //TODO: add temporal mv predictor
  }
#endif

  while (candidates < 2) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    candidates++;
  }

}