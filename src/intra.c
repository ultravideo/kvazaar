/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file intra.c
    \brief Intra functions
    \author Marko Viitanen
    \date 2013-03
    
    Intra functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "picture.h"
#include "intra.h"

/*!
 \brief Set intrablock mode (and init typedata)
 \param pic picture to use
 \param xCtb x CU position (smallest CU)
 \param yCtb y CU position (smallest CU)
 \param depth current CU depth
 \param mode mode to set
 \returns Void
*/
void intra_setBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, uint8_t mode)
{
  int x,y,d;
  //Width in smallest CU
  int width_in_SCU = pic->width/(LCU_WIDTH>>MAX_DEPTH);
  int block_SCU_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  for(y = yCtb; y < yCtb+block_SCU_width; y++)
  {
    for(x = xCtb; x < xCtb+block_SCU_width; x++)
    {
      int CUpos = y*width_in_SCU+x;
      for(d = 0; d < MAX_DEPTH; d++)
      {
        pic->CU[d][CUpos].type = CU_INTRA;
        if(pic->CU[d][CUpos].typedata == NULL)
        {
          pic->CU[d][CUpos].typedata = (void *)malloc(sizeof(CU_info_intra));
          memset(pic->CU[d][CUpos].typedata, 0, sizeof(CU_info_intra));
        }
        ((CU_info_intra*)pic->CU[d][CUpos].typedata)->mode = mode;
      }
    }
  }
}

/*!
 \brief get intrablock mode
 \param pic picture to use
 \param xCtb x CU position (smallest CU)
 \param yCtb y CU position (smallest CU)
 \param depth current CU depth
 \returns mode if it's present, otherwise -1
*/
int8_t intra_getBlockMode(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth)
{
  //Width in smallest CU
  int width_in_SCU = pic->width/(LCU_WIDTH>>MAX_DEPTH);
  int CUpos = yCtb*width_in_SCU+xCtb;
  if(pic->CU[depth][CUpos].typedata != NULL)
  {
    return ((CU_info_intra*)pic->CU[depth][CUpos].typedata)->mode;
  }  
  return -1;
}

/*!
\brief get intrablock mode
\param pic picture to use
\param xCtb x CU position (smallest CU)
\param yCtb y CU position (smallest CU)
\param depth current CU depth
\returns DC prediction or -1 if not available
*/
int16_t intra_getDCPred(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth)
{
  int32_t i, iSum = 0;
  int16_t pDcVal = 1<<(g_uiBitDepth-1);
  int32_t width = (LCU_WIDTH>>depth);
  int32_t SCUwidth = (LCU_WIDTH>>MAX_DEPTH);
  int32_t xpos = SCUwidth*xCtb;
  int32_t ypos = SCUwidth*yCtb;
  int32_t LCU_xpos = xpos - (xpos%LCU_WIDTH);
  int32_t LCU_ypos = ypos - (ypos%LCU_WIDTH);

  int8_t bAbove = (yCtb>0)?1:0;//*/(yCtb*SCUwidth>(LCU_WIDTH-1))?1:0;
  int8_t bLeft  = (xCtb>0)?1:0;//*/(xCtb*SCUwidth>(LCU_WIDTH-1))?1:0;
  int32_t add;


  if (bAbove)
  {
    add = (ypos-1)*pic->width;
    for (i = xpos+add; i < xpos+width+add ; i++)
    {
      iSum += pic->yRecData[i];
    }
  }  
  else if (bLeft)
  {
    iSum += pic->yRecData[ypos*pic->width+xpos-1]*width;
  }


  if (bLeft)
  {
    add = xpos-1;
    for (i = ypos*pic->width+add ; i < (ypos+width)*pic->width+add ; i+=pic->width)
    {
      iSum += pic->yRecData[i];
    }
  }  
  else if (bAbove)
  {
    iSum += pic->yRecData[(ypos-1)*pic->width+xpos]*width;
  }



  //if (bAbove && bLeft)
  if (bAbove || bLeft)
  {
    pDcVal = (iSum + width) / (width + width);
  }
  /*
  else if (bAbove)
  {
    pDcVal = (iSum + width/2) / width;
  }
  else if (bLeft)
  {
    pDcVal = (iSum + width/2) / width;
  }
  */
  
  return pDcVal;
}

/*! \brief Function for deriving intra luma predictions
  \param pic picture to use
  \param xCtb x CU position (smallest CU)
  \param yCtb y CU position (smallest CU)
  \param depth current CU depth
  \param preds output buffer for 3 predictions 
  \returns (predictions are found)?1:0
*/
int8_t intra_getDirLumaPredictor(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t* preds)
{
  int32_t iLeftIntraDir  = 1; //DC_IDX
  int32_t iAboveIntraDir = 1; //DC_IDX
  int32_t width_in_SCU = pic->width/(LCU_WIDTH>>MAX_DEPTH);
  int32_t CUpos = yCtb*width_in_SCU+xCtb;
  
  // Left PU predictor
  if(xCtb != 0 && pic->CU[depth][CUpos-1].type == CU_INTRA)
  {
    iLeftIntraDir = ((CU_info_intra*)pic->CU[depth][CUpos-1].typedata)->mode;
  }

  // Top PU predictor
  if(yCtb != 0 && pic->CU[depth][CUpos-width_in_SCU].type == CU_INTRA)
  {
    iAboveIntraDir = ((CU_info_intra*)pic->CU[depth][CUpos-width_in_SCU].typedata)->mode;
  }

  if(iLeftIntraDir == iAboveIntraDir)
  {
  
    if (iLeftIntraDir > 1) // angular modes
    {
      preds[0] = iLeftIntraDir;
      preds[1] = ((iLeftIntraDir + 29) % 32) + 2;
      preds[2] = ((iLeftIntraDir - 1 ) % 32) + 2;
    }
    else //non-angular
    {
      preds[0] = 0;//PLANAR_IDX;
      preds[1] = 1;//DC_IDX;
      preds[2] = 26;//VER_IDX; 
    }
  }
  else
  {
    preds[0] = iLeftIntraDir;
    preds[1] = iAboveIntraDir;
    
    if (iLeftIntraDir && iAboveIntraDir ) //both modes are non-planar
    {
      preds[2] = 0;//PLANAR_IDX;
    }
    else
    {
      preds[2] =  (iLeftIntraDir+iAboveIntraDir)<2? 26 : 1;
    }
  } 

  return 1;
}

/*! \brief Function for deriving planar intra prediction.
  \param pic picture to use
  \param xCtb x CU position (smallest CU)
  \param yCtb y CU position (smallest CU)
  \param depth current CU depth
  \param dst destination buffer for prediction
 
  This function derives the prediction samples for planar mode (intra coding).
*/
/*
void intra_getPlanarPred(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int16_t* dst)
{
  int32_t k, l, bottomLeft, topRight;
  int32_t horPred;
  int32_t width = (LCU_WIDTH>>depth);
  int16_t leftColumn[LCU_WIDTH], topRow[LCU_WIDTH], bottomRow[LCU_WIDTH], rightColumn[LCU_WIDTH];
  uint32_t blkSize = width;
  uint32_t offset2D = width;
  uint32_t shift1D = g_aucConvertToBit[ width ] + 2;
  uint32_t shift2D = shift1D + 1;

  // Get left and above reference column and row
  for(k=0;k<blkSize+1;k++)
  {
    topRow[k] = pSrc[k-srcStride];
    leftColumn[k] = pSrc[k*srcStride-1];
  }

  // Prepare intermediate variables used in interpolation
  bottomLeft = leftColumn[blkSize];
  topRight   = topRow[blkSize];
  for (k=0;k<blkSize;k++)
  {
    bottomRow[k]   = bottomLeft - topRow[k];
    rightColumn[k] = topRight   - leftColumn[k];
    topRow[k]      <<= shift1D;
    leftColumn[k]  <<= shift1D;
  }

  // Generate prediction signal
  for (k=0;k<blkSize;k++)
  {
    horPred = leftColumn[k] + offset2D;
    for (l=0;l<blkSize;l++)
    {
      horPred += rightColumn[k];
      topRow[l] += bottomRow[l];
      rpDst[k*dstStride+l] = ( (horPred + topRow[l]) >> shift2D );
    }
  }
}
*/