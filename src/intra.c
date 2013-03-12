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
        pic->CU[d][CUpos].intra.mode = mode;
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
  if(pic->CU[depth][CUpos].type == CU_INTRA)
  {
    return pic->CU[depth][CUpos].intra.mode;
  }  
  return -1;
}

/*!
\brief get intrablock mode
\param pic picture to use
\param xpos x-position
\param ypos y-position
\param width block width
\returns DC prediction or -1 if not available
*/
int16_t intra_getDCPred(uint8_t* pic, uint16_t picwidth,uint32_t xpos, uint32_t ypos, uint8_t width)
{
  int32_t i, iSum = 0;
  int16_t pDcVal = 1<<(g_uiBitDepth-1);  
  int8_t bAbove = xpos?1:0;
  int8_t bLeft  = ypos?1:0;
  int32_t add;


  if (bAbove)
  {
    add = (ypos-1)*picwidth;
    for (i = xpos+add; i < xpos+width+add ; i++)
    {
      iSum += pic[i];
    }
  }  
  else if (bLeft)
  {
    iSum += pic[ypos*picwidth+xpos-1]*width;
  }


  if (bLeft)
  {
    add = xpos-1;
    for (i = ypos*picwidth+add ; i < (ypos+width)*picwidth+add ; i+=picwidth)
    {
      iSum += pic[i];
    }
  }  
  else if (bAbove)
  {
    iSum += pic[(ypos-1)*picwidth+xpos]*width;
  }

  if (bAbove || bLeft)
  {
    pDcVal = (iSum + width) / (width + width);
  }
  
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
  if(xCtb && pic->CU[depth][CUpos-1].type == CU_INTRA)
  {
    iLeftIntraDir = pic->CU[depth][CUpos-1].intra.mode;
  }

  // Top PU predictor
  if(yCtb && ((yCtb*(LCU_WIDTH>>MAX_DEPTH))%LCU_WIDTH)!=0 && pic->CU[depth][CUpos-width_in_SCU].type == CU_INTRA)
  {
    iAboveIntraDir = pic->CU[depth][CUpos-width_in_SCU].intra.mode;
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

/*! \brief Function to test best intra prediction
  \param orig original picture data
  \param origstride original picture stride
  \param rec reconstructed picture data
  \param recstride reconstructed picture stride
  \param xpos source x-position
  \param ypos source y-position
  \param width block size to predict
  \param dst destination buffer for best prediction
  \param dststride destination width
  \param sad sad value of best mode
  \returns best intra mode
 
  This function derives the prediction samples for planar mode (intra coding).
*/
int16_t intra_prediction(uint8_t* orig,uint32_t origstride,uint8_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, uint32_t *sad)
{
  uint32_t bestSAD = 0xffffffff;
  uint32_t SAD = 0;
  int16_t bestMode = 1;
  int32_t x,y,i;
  uint32_t (*calcSAD)(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
  int16_t pred[LCU_WIDTH*LCU_WIDTH>>2];
  int16_t origBlock[LCU_WIDTH*LCU_WIDTH>>2];
  uint8_t *origShift = &orig[xpos+ypos*origstride];
  switch(width)
  {
    case 16:
      calcSAD = &SAD16x16;
      break;
    case 32:
      calcSAD = &SAD32x32;
      break;
    default:
  }
  i = 0;
  for(y = 0; y < width; y++)
  {
    for(x = 0; x < width; x++)
    {
      origBlock[i++] = origShift[x+y*origstride];
    }
  }

  /* Test planar */  
  intra_getPlanarPred(rec, recstride, xpos, ypos, width, pred, width);
  SAD = calcSAD(pred,width,origBlock,width);
  bestSAD = SAD; /* First to be tested */
  bestMode = 1;

  /* Test DC */  
  x = intra_getDCPred(rec, recstride, xpos, ypos, width);
  for(i = width*width-1; i >= 0; i--)
  {
    pred[i] = x;
  }
  SAD = calcSAD(pred,width,origBlock,width);
  if(SAD < bestSAD)
  {
    bestSAD = SAD;
    bestMode = 0;
  }

  /* Test directional predictions */
  
}

void intra_DCPredFiltering(uint8_t* pSrc, int32_t iSrcStride, uint8_t* rpDst, int32_t iDstStride, int32_t iWidth, int32_t iHeight )
{
  uint8_t* pDst = rpDst;
  int32_t x, y, iDstStride2, iSrcStride2;

  // boundary pixels processing
  pDst[0] = ((pSrc[-iSrcStride] + pSrc[-1] + 2 * pDst[0] + 2) >> 2);

  for ( x = 1; x < iWidth; x++ )
  {
    pDst[x] = ((pSrc[x - iSrcStride] +  3 * pDst[x] + 2) >> 2);
  }
  for ( y = 1, iDstStride2 = iDstStride, iSrcStride2 = iSrcStride-1; y < iHeight; y++, iDstStride2+=iDstStride, iSrcStride2+=iSrcStride )
  {
    pDst[iDstStride2] = ((pSrc[iSrcStride2] + 3 * pDst[iDstStride2] + 2) >> 2);
  }

  return;
}

/*! \brief Function for deriving planar intra prediction.
  \param src source pixel array
  \param srcstride source width
  \param xpos source x-position
  \param ypos source y-position
  \param width block size to predict
  \param dst destination buffer for prediction
  \param dststride destination width
 
  This function derives the prediction samples for planar mode (intra coding).
*/
void intra_getPlanarPred(uint8_t* src,uint32_t srcstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride)
{
  int8_t bAbove = ypos?1:0;
  int8_t bLeft  = xpos?1:0;
  int16_t pDcVal = 1<<(g_uiBitDepth-1);
  uint8_t* srcShifted = &src[xpos + ypos*srcstride];
  int32_t k, l, bottomLeft, topRight;
  int32_t horPred;
  int16_t leftColumn[LCU_WIDTH], topRow[LCU_WIDTH], bottomRow[LCU_WIDTH], rightColumn[LCU_WIDTH];
  uint32_t blkSize = width;
  uint32_t offset2D = width;
  uint32_t shift1D = g_aucConvertToBit[ width ] + 2;
  uint32_t shift2D = shift1D + 1;

  if(bAbove)
  {    
    for(k=0;k<blkSize;k++)
    {
      topRow[k] = srcShifted[k-srcstride];
    }
  }
  else
  {
    int16_t prediction = bLeft?srcShifted[-1]:pDcVal;
    for(k=0;k<blkSize;k++)
    {
      topRow[k] = prediction;
    }
  }

  if(bLeft)
  {
    for(k=0;k<blkSize;k++)
    {
      leftColumn[k] = srcShifted[k*srcstride-1];
    }
  }
  else
  {
    int16_t prediction = bAbove?srcShifted[-srcstride]:pDcVal;
    for(k=0;k<blkSize;k++)
    {
      leftColumn[k] = prediction;
    }
  }
  leftColumn[blkSize] = leftColumn[blkSize-1];
  topRow[blkSize] = topRow[blkSize-1];

  // Get left and above reference column and row
  

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
      dst[k*dststride+l] = ( (horPred + topRow[l]) >> shift2D );
    }
  }
}
