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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "intra.h"
#include "encoder.h"
#include "filter.h"
#include "search.h"



void search_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth)
{   
  uint8_t border_x = ((encoder->in.width)<( xCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_y = ((encoder->in.height)<( yCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */
  CU_info *cur_CU = &encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)];

  cur_CU->intra.cost = (uint32_t)-1;

  /* Force split on border */
  if(depth != MAX_DEPTH)
  {
    if(border)
    {
      /* Split blocks and remember to change x and y block positions */
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      SET_SPLITDATA(cur_CU,1);
      search_tree(encoder,xCtb,yCtb,depth+1);
      if(!border_x)
      {
        search_tree(encoder,xCtb+change,yCtb,depth+1);
      }
      if(!border_y)
      {
        search_tree(encoder,xCtb,yCtb+change,depth+1);      
      }
      if(!border)
      {
        search_tree(encoder,xCtb+change,yCtb+change,depth+1);
      }
      /* We don't need to do anything else here */
      return;
    }
  }

  
  if(encoder->in.cur_pic.slicetype != SLICE_I)
  {


  }

  /* INTRA SEARCH */
  if(depth > 0)
  {
    uint8_t *base  = &encoder->in.cur_pic.yData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
    uint8_t *baseU = &encoder->in.cur_pic.uData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
    uint8_t *baseV = &encoder->in.cur_pic.vData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
    uint32_t width = LCU_WIDTH>>depth;

    /* INTRAPREDICTION */
    /* ToDo: split to a function */
    int16_t pred[LCU_WIDTH*LCU_WIDTH];
    int16_t predU[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t predV[LCU_WIDTH*LCU_WIDTH>>2];

    int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    int16_t *recShift  = &rec[(LCU_WIDTH>>(depth))*2+8+1];
    int16_t *recShiftU = &rec[(LCU_WIDTH>>(depth+1))*2+8+1];
    uint8_t *recbase   = &encoder->in.cur_pic.yRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
    uint8_t *recbaseU  = &encoder->in.cur_pic.uRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
    uint8_t *recbaseV  = &encoder->in.cur_pic.vRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];

    /* Build reconstructed block to use in prediction with extrapolated borders */      
    intra_buildReferenceBorder(&encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth))*2+8, rec, (LCU_WIDTH>>(depth))*2+8, 0);
    cur_CU->intra.mode = (uint8_t)intra_prediction(encoder->in.cur_pic.yData,encoder->in.width,recShift,(LCU_WIDTH>>(depth))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH)),yCtb*(LCU_WIDTH>>(MAX_DEPTH)),width,pred,width,&cur_CU->intra.cost);
  }

  /* Split and search to max_depth */
  if(depth != 2)
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
  CU_info *cur_CU = &encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)];
  uint32_t bestCost = cur_CU->intra.cost;
  int8_t bestMode = cur_CU->type;
  uint32_t cost = 0;

  /* Split and search to max_depth */
  if(depth != MAX_DEPTH)
  {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);  
    cost = 4*g_lambda_cost[encoder->QP];
    cost += search_best_mode(encoder,xCtb,yCtb,depth+1);
    cost += search_best_mode(encoder,xCtb+change,yCtb,depth+1);
    cost += search_best_mode(encoder,xCtb,yCtb+change,depth+1);      
    cost += search_best_mode(encoder,xCtb+change,yCtb+change,depth+1);

    if(cost != 0 && cost < bestCost)
    {
      cur_CU->split = 1;
      bestCost = cost;
    }
    else
    {
      cur_CU->split = 0;
    }
  }

  return bestCost;
}

void search_slice_data(encoder_control* encoder)
{
  uint16_t xCtb,yCtb;

  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_LCU; yCtb++)
  {
    uint8_t lastCUy = (yCtb == (encoder->in.height_in_LCU-1))?1:0;
    for(xCtb = 0; xCtb < encoder->in.width_in_LCU; xCtb++)
    {
      uint8_t lastCUx = (xCtb == (encoder->in.width_in_LCU-1))?1:0;
      uint8_t depth = 0;

      /* Recursive function for looping through all the sub-blocks */
      search_tree(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);

      /* Decide actual coding modes */
      search_best_mode(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);
    }
  }
}