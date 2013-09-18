/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file search.c
    \brief searching
    \author Marko Viitanen
    \date 2013-04
    
    Search related functions
*/

#include "search.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "debug.h"


// Temporarily for debugging.
#define USE_INTRA_IN_P 0
#define RENDER_CU 1

/**
 * 
 *
 * pic: 
 * pic_data: picture color data starting from the block MV is being searched for.
 * ref_data: picture color data starting from the beginning of reference pic.
 * cur_cu: 
 */
void search_motion_vector(picture *pic, uint8_t *pic_data, uint8_t *ref_data, CU_info *cur_cu, unsigned step, int orig_x, int orig_y, int x, int y, unsigned depth)
{
  // TODO: Inter: Handle non-square blocks.
  unsigned block_width = CU_WIDTH_FROM_DEPTH(depth);
  unsigned block_height = block_width;
  unsigned cost;

  // TODO: Inter: Calculating error outside picture borders.
  // This prevents choosing vectors that need interpolating of borders to work.
  if (orig_x + x < 0 || orig_y + y < 0 || orig_x + x > pic->width - block_width || orig_y + y > pic->height - block_height) return;

  cost = SAD(pic_data, &ref_data[(orig_y + y) * pic->width + (orig_x + x)], block_width, block_height, pic->width) + 1;
  if (cost < cur_cu->inter.cost) {
    cur_cu->inter.cost = cost;
    cur_cu->inter.mv[0] = x << 2;
    cur_cu->inter.mv[1] = y << 2;
  }

  step /= 2;
  if (step > 0) {
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y, x, y - step, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y, x - step, y, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y, x + step, y, depth);
    search_motion_vector(pic, pic_data, ref_data, cur_cu, step, orig_x, orig_y, x, y + step, depth);
  }
}

void search_buildReferenceBorder(picture* pic, int32_t xCtb, int32_t yCtb,int16_t outwidth, int16_t* dst, int32_t dststride, int8_t chroma)
{
  int32_t leftColumn;  /*!< left column iterator */
  int16_t val;         /*!< variable to store extrapolated value */
  int32_t i;           /*!< index iterator */
  int16_t dcVal        = 1<<(g_bitDepth-1); /*!< default predictor value */
  int32_t topRow;      /*!< top row iterator */
  int32_t srcWidth     = (pic->width>>(chroma?1:0)); /*!< source picture width */
  int32_t srcHeight    = (pic->height>>(chroma?1:0));/*!< source picture height */
  uint8_t* srcPic      = (!chroma)?pic->yData: ((chroma==1)?pic->uData: pic->vData); /*!< input picture pointer */  
  int16_t SCU_width    = LCU_WIDTH>>(MAX_DEPTH+(chroma?1:0)); /*!< Smallest Coding Unit width */
  uint8_t* srcShifted  = &srcPic[xCtb*SCU_width+(yCtb*SCU_width)*srcWidth];  /*!< input picture pointer shifted to start from the left-top corner of the current block */
  int32_t width_in_SCU = pic->width_in_LCU<<MAX_DEPTH;     /*!< picture width in SCU */

  /* Fill left column */
  if(xCtb)
  {
    /* Loop SCU's */
    for(leftColumn = 1; leftColumn < outwidth/SCU_width; leftColumn++)
    {
      /* If over the picture height or block not yet searched, stop */
      if((yCtb+leftColumn)*SCU_width >= srcHeight || pic->CU[0][xCtb-1+(yCtb+leftColumn)*width_in_SCU].type == CU_NOTSET)
      {
        break;
      }
    }
    /* Copy the pixels to output */
    for(i = 0; i < leftColumn*SCU_width-1; i ++)
    {
      dst[(i+1)*dststride] = srcShifted[i*srcWidth-1];
    }

    /* if the loop was not completed, extrapolate the last pixel pushed to output */
    if(leftColumn != outwidth/SCU_width)
    {
      val = srcShifted[(leftColumn*SCU_width-1)*srcWidth-1];
      for(i = (leftColumn*SCU_width); i < outwidth; i++)
      {
        dst[i*dststride] = val;
      }
    }    
  }
  /* If left column not available, copy from toprow or use the default predictor */
  else
  {
    val = yCtb?srcShifted[-srcWidth]:dcVal;
    for(i = 0; i < outwidth; i++)
    {
      dst[i*dststride] = val;
    }
  }

  if(yCtb)
  {
    /* Loop top SCU's */
    for(topRow = 1; topRow < outwidth/SCU_width; topRow++)
    {
      if((xCtb+topRow)*SCU_width >= srcWidth || pic->CU[0][xCtb+topRow+(yCtb-1)*width_in_SCU].type == CU_NOTSET)
      {
        break;
      }
    }

    for(i = 0; i < topRow*SCU_width-1; i ++)
    {
      dst[i+1] = srcShifted[i-srcWidth];
    }

    if(topRow != outwidth/SCU_width)
    {
      val = srcShifted[(topRow*SCU_width)-srcWidth-1];
      for(i = (topRow*SCU_width); i < outwidth; i++)
      {
        dst[i] = val;
      }
    }
  }
  else
  {
    val = xCtb?srcShifted[-1]:dcVal;
    for(i = 1; i < outwidth; i++)
    {
      dst[i] = val;
    }
  }
  /* Topleft corner */
  dst[0] = (xCtb&&yCtb)?srcShifted[-srcWidth-1]:dst[dststride];

}

void search_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth)
{   
  uint8_t border_x = ((encoder->in.width)<( xCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_y = ((encoder->in.height)<( yCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_split_x = ((encoder->in.width)  < ( (xCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;
  uint8_t border_split_y = ((encoder->in.height) < ( (yCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */
  CU_info *cur_CU = &encoder->in.cur_pic->CU[depth][xCtb+yCtb*(encoder->in.width_in_LCU<<MAX_DEPTH)];

  cur_CU->intra.cost = 0xffffffff;
  cur_CU->inter.cost = 0xffffffff;

  /* Force split on border */
  if(depth != MAX_DEPTH)
  {
    if(border)
    {
      /* Split blocks and remember to change x and y block positions */
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      SET_SPLITDATA(cur_CU,1);
      search_tree(encoder,xCtb,yCtb,depth+1);
      if(!border_x || border_split_x)
      {
        search_tree(encoder,xCtb+change,yCtb,depth+1);
      }
      if(!border_y || border_split_y)
      {
        search_tree(encoder,xCtb,yCtb+change,depth+1);      
      }
      if(!border || (border_split_x && border_split_y) )
      {
        search_tree(encoder,xCtb+change,yCtb+change,depth+1);
      }
      /* We don't need to do anything else here */
      return;
    }
  }
  
  /* INTER SEARCH */
  if(depth >= MIN_INTER_SEARCH_DEPTH
      && depth <= MAX_INTER_SEARCH_DEPTH 
      && encoder->in.cur_pic->slicetype != SLICE_I) {
    /* Motion estimation on P-frame */
    if(encoder->in.cur_pic->slicetype != SLICE_B)
    {

    }

    {
      unsigned mv[2] = { 0, 0 }; // TODO: Take initial MV from adjacent blocks.
      picture *cur_pic = encoder->in.cur_pic;
        
      picture *ref_pic = encoder->ref->pics[0];

      int x = xCtb * CU_MIN_SIZE_PIXELS;
      int y = yCtb * CU_MIN_SIZE_PIXELS;
      uint8_t *cur_data = &cur_pic->yData[(y * cur_pic->width) + x];
      search_motion_vector(cur_pic, cur_data, ref_pic->yData, cur_CU, 8, x, y, 0, 0, depth);
    }

    cur_CU->type = CU_INTER;
    cur_CU->inter.mv_dir = 1;
    inter_setBlockMode(encoder->in.cur_pic,xCtb,yCtb,depth,cur_CU);
  }

  /* INTRA SEARCH */
  if (depth >= MIN_INTRA_SEARCH_DEPTH 
      && depth <= MAX_INTRA_SEARCH_DEPTH 
      && (encoder->in.cur_pic->slicetype == SLICE_I || USE_INTRA_IN_P)) {
    int x = 0,y = 0;
    uint8_t *base  = &encoder->in.cur_pic->yData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
    uint32_t width = LCU_WIDTH>>depth;

    /* INTRAPREDICTION */
    int16_t pred[LCU_WIDTH*LCU_WIDTH+1];
    int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    int16_t *recShift = &rec[(LCU_WIDTH>>(depth))*2+8+1];

    //int16_t *pred = (int16_t*)malloc(LCU_WIDTH*LCU_WIDTH*sizeof(int16_t));
    //int16_t *rec = (int16_t*)malloc((LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)*sizeof(int16_t));

    /* Build reconstructed block to use in prediction with extrapolated borders */
    search_buildReferenceBorder(encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth))*2+8, rec, (LCU_WIDTH>>(depth))*2+8, 0);
    cur_CU->intra.mode = (uint8_t)intra_prediction(encoder->in.cur_pic->yData,encoder->in.width,recShift,(LCU_WIDTH>>(depth))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH)),yCtb*(LCU_WIDTH>>(MAX_DEPTH)),width,pred,width,&cur_CU->intra.cost);
    //free(pred);
    //free(rec);
  }

  /* Split and search to max_depth */
  if(depth < MAX_INTRA_SEARCH_DEPTH && depth < MAX_INTER_SEARCH_DEPTH)
  {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);
    search_tree(encoder,xCtb,yCtb,depth+1);
    search_tree(encoder,xCtb+change,yCtb,depth+1);
    search_tree(encoder,xCtb,yCtb+change,depth+1);      
    search_tree(encoder,xCtb+change,yCtb+change,depth+1);
  }
}

uint32_t search_best_mode(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth)
{
  CU_info *cur_CU = &encoder->in.cur_pic->CU[depth][xCtb+yCtb*(encoder->in.width_in_LCU<<MAX_DEPTH)];
  uint32_t bestIntraCost = cur_CU->intra.cost;
  uint32_t bestInterCost = cur_CU->inter.cost;
  uint32_t bestCost = 0;
  uint32_t cost = 0;
  uint32_t lambdaCost = (4 * g_lambda_cost[encoder->QP]) << 4; //<<5; //TODO: Correct cost calculation
  
  /* Split and search to max_depth */
  if (depth != MAX_INTRA_SEARCH_DEPTH) {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);
    cost = search_best_mode(encoder,xCtb,yCtb,depth+1);
    cost += search_best_mode(encoder,xCtb+change,yCtb,depth+1);
    cost += search_best_mode(encoder,xCtb,yCtb+change,depth+1);
    cost += search_best_mode(encoder,xCtb+change,yCtb+change,depth+1);

    /* We split if the cost is better (0 cost -> not checked) */
    if(cost != 0 &&   (bestIntraCost != 0 && cost+lambdaCost < bestIntraCost) && 
                      (bestInterCost != 0 && cost+lambdaCost < bestInterCost && encoder->in.cur_pic->slicetype != SLICE_I))
    {
      /* Set split to 1 */
      picture_setBlockSplit(encoder->in.cur_pic,xCtb,yCtb,depth,1);
      bestCost = cost+lambdaCost;
    }
    /* Else, check if inter cost is smaller or the same as intra */
    else if(bestInterCost != 0 && (bestInterCost <= bestIntraCost || bestIntraCost == 0) && encoder->in.cur_pic->slicetype != SLICE_I)
    {
      /* Set split to 0 and mode to inter.mode */
      picture_setBlockSplit(encoder->in.cur_pic,xCtb,yCtb,depth,0);
      inter_setBlockMode(encoder->in.cur_pic,xCtb,yCtb,depth,cur_CU);
      bestCost = bestInterCost;
    }
    /* Else, dont split and recursively set block mode */
    else
    {
      /* Set split to 0 and mode to intra.mode */
      picture_setBlockSplit(encoder->in.cur_pic,xCtb,yCtb,depth,0);
      intra_setBlockMode(encoder->in.cur_pic,xCtb,yCtb,depth,cur_CU->intra.mode);
      bestCost = bestIntraCost;
    }
  }
  else if(bestInterCost != 0 && (bestInterCost <= bestIntraCost || bestIntraCost == 0) && encoder->in.cur_pic->slicetype != SLICE_I)
  {
    /* Set split to 0 and mode to inter.mode */
    picture_setBlockSplit(encoder->in.cur_pic,xCtb,yCtb,depth,0);
    inter_setBlockMode(encoder->in.cur_pic,xCtb,yCtb,depth,cur_CU);
    bestCost = bestInterCost;
  }
  else
  {
    /* Set split to 0 and mode to intra.mode */
    picture_setBlockSplit(encoder->in.cur_pic,xCtb,yCtb,depth,0);
    intra_setBlockMode(encoder->in.cur_pic,xCtb,yCtb,depth,cur_CU->intra.mode);
    bestCost = bestIntraCost;
  }

  return bestCost;
}

void search_slice_data(encoder_control* encoder)
{
  int16_t xCtb,yCtb;
  FILE *fp = 0, *fp2 = 0;

  if (RENDER_CU && encoder->frame == 1) {
    fp = open_cu_file("cu_search.html");
    fp2 = open_cu_file("cu_best.html");
  }


  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_LCU; yCtb++)
  {
    for(xCtb = 0; xCtb < encoder->in.width_in_LCU; xCtb++)
    {
      uint8_t depth = 0;
      /* Recursive function for looping through all the sub-blocks */
      search_tree(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);
      if (RENDER_CU && encoder->frame == 1) render_cu_file(encoder, depth, xCtb<<MAX_DEPTH, yCtb<<MAX_DEPTH, fp);

      /* Decide actual coding modes */
      search_best_mode(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);
      if (RENDER_CU && encoder->frame == 1) render_cu_file(encoder, depth, xCtb<<MAX_DEPTH, yCtb<<MAX_DEPTH, fp2);
    }
  }

  if (RENDER_CU && fp) {
    close_cu_file(fp);
    fp = 0;
  }
  if (RENDER_CU && fp2) {
    close_cu_file(fp2);
    fp2 = 0;
  }
}
