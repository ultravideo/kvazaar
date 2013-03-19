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


const uint8_t intraHorVerDistThres[4] = {0,7,1,0};

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
  uint32_t x,y,d;
  //Width in smallest CU
  int width_in_SCU = pic->width/(LCU_WIDTH>>MAX_DEPTH);
  int block_SCU_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  for(y = yCtb; y < yCtb+block_SCU_width; y++)
  {
    int CUpos = y*width_in_SCU;
    for(x = xCtb; x < xCtb+block_SCU_width; x++)
    {      
      for(d = 0; d < MAX_DEPTH; d++)
      {
        pic->CU[d][CUpos+x].type = CU_INTRA;
        pic->CU[d][CUpos+x].intra.mode = mode;
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
\returns DC prediction
*/
int16_t intra_getDCPred(int16_t* pic, uint16_t picwidth,uint32_t xpos, uint32_t ypos, uint8_t width)
{
  int32_t i, iSum = 0;
  int16_t pDcVal = 1<<(g_uiBitDepth-1);  
  
  /* Average of pixels on top and left */
  for (i = -picwidth; i < width-picwidth ; i++)
  {
    iSum += pic[i];
  }

  for (i = -1 ; i < width*picwidth-1 ; i+=picwidth)
  {
    iSum += pic[i];
  }

  pDcVal = (iSum + width) / (width + width);
  
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


void intra_filter(int16_t* ref, uint32_t stride,uint32_t width, int8_t mode)
{
  #define FWIDTH (LCU_WIDTH+1)
  int16_t filtered[FWIDTH*FWIDTH];
  int16_t* filteredShift = &filtered[FWIDTH+1];
  int x,y;
  if(!mode)
  {
    //pF[ -1 ][ -1 ] = ( p[ -1 ][ 0 ] + 2 * p[ -1 ][ -1 ] + p[ 0 ][ -1 ] + 2 )  >>  2	(8 35)
    filteredShift[-FWIDTH-1] = (ref[-1]+2*ref[-(int32_t)stride-1]+ref[-(int32_t)stride]+2) >> 2;

    //pF[ -1 ][ y ] = ( p[ -1 ][ y + 1 ] + 2 * p[ -1 ][ y ] + p[ -1 ][ y - 1 ] + 2 )  >>  2 for y = 0..nTbS * 2 - 2	(8 36)
    for(y = 0; y < (int32_t)width*2-1; y++)
    {
      filteredShift[y*FWIDTH-1] = (ref[(y+1)*stride-1] + 2*ref[y*stride-1] + ref[(y-1)*stride-1]+2) >> 2;
    }
    
    //pF[ -1 ][ nTbS * 2 - 1 ] = p[ -1 ][ nTbS * 2 - 1 ]		(8 37)
    filteredShift[(width*2-1)*FWIDTH-1] = ref[(width*2-1)*stride-1];

    //pF[ x ][ -1 ] = ( p[ x - 1 ][ -1 ] + 2 * p[ x ][ -1 ] + p[ x + 1 ][ -1 ] + 2 )  >>  2 for x = 0..nTbS * 2 - 2	(8 38)
    for(x = 0; x < (int32_t)width*2-1; x++)
    {
      filteredShift[x-FWIDTH] = (ref[x-1-stride] + 2*ref[x-stride] + ref[x+1-stride]+2) >> 2;
    }

    //pF[ nTbS * 2 - 1 ][ -1 ] = p[ nTbS * 2 - 1 ][ -1 ]	
    filteredShift[(width*2-1)-FWIDTH] = ref[(width*2-1)-stride];


    /* Copy filtered samples to the input array */
    for(x = -1; x < (int32_t)width*2; x++)
    {
      ref[x-stride] = filtered[x+1];
    }
    for(y = 0; y < (int32_t)width*2; y++)
    {
      ref[y*stride-1] = filtered[(y+1)*FWIDTH];
    }
  }
  else
  {
    printf("UNHANDLED: %s: %d\r\n", __FILE__, __LINE__);
    exit(1);
  }
  #undef FWIDTH
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
int16_t intra_prediction(uint8_t* orig,uint32_t origstride,int16_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, uint32_t *sad)
{
  typedef uint32_t (*SADfunction)(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
  uint32_t bestSAD = 0xffffffff;
  uint32_t SAD = 0;
  int16_t bestMode = 1;
  int32_t x,y,i;
  uint32_t (*calcSAD)(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
  int16_t pred[LCU_WIDTH*LCU_WIDTH>>2];
  int16_t origBlock[LCU_WIDTH*LCU_WIDTH>>2];
  uint8_t *origShift = &orig[xpos+ypos*origstride];  
  SADfunction SADarray[4] = {&SAD4x4,&SAD8x8,&SAD16x16,&SAD32x32};
  uint8_t threshold = intraHorVerDistThres[g_toBits[width]]; /*!< Intra filtering threshold */
  #define COPY_PRED_TO_DST() for(y = 0; y < (int32_t)width; y++)  {   for(x = 0; x < (int32_t)width; x++)  {  dst[x+y*dststride] = pred[x+y*width];  }   }
  #define CHECK_FOR_BEST(mode)  SAD = calcSAD(pred,width,origBlock,width); \
                                if(SAD < bestSAD)\
                                {\
                                  bestSAD = SAD;\
                                  bestMode = mode;\
                                  COPY_PRED_TO_DST();\
                                }

  /* Choose SAD function according to width */
  calcSAD = SADarray[g_toBits[width]];

  /* Store original block for SAD computation */
  i = 0;
  for(y = 0; y < (int32_t)width; y++)
  {
    for(x = 0; x < (int32_t)width; x++)
    {
      origBlock[i++] = origShift[x+y*origstride];
    }
  }


  /* Test DC */  
  x = intra_getDCPred(rec, recstride, xpos, ypos, width);
  for(i = 0; i < (int32_t)(width*width); i++)
  {
    pred[i] = x;
  }
  CHECK_FOR_BEST(1);

  /* Check angular that don't require filtering */
  for(i = 2; i < 35; i++)
  {
    if(MIN(abs(i-26),abs(i-10)) <= threshold)
    {
      intra_getAngularPred(rec,recstride,pred, width,width,width,i, xpos?1:0, ypos?1:0, 0);
      CHECK_FOR_BEST(i);
    }
  }
  
  /*Apply filter*/
  intra_filter(rec,recstride,width,0);

  /* Test planar */  
  intra_getPlanarPred(rec, recstride, xpos, ypos, width, pred, width);
  CHECK_FOR_BEST(0);
  

  /* Test directional predictions */
  /* ToDo: add conditions to skip some modes on borders */
  
  //chroma can use only 26 and 10
  /* Test angular predictions which require filtered samples */
  for(i = 2; i < 35; i++)
  {
    if(MIN(abs(i-26),abs(i-10)) > threshold)
    {
      intra_getAngularPred(rec,recstride,pred, width,width,width,i, xpos?1:0, ypos?1:0, 0);
      CHECK_FOR_BEST(i);
    }
  }

  *sad = bestSAD;
  #undef COPY_PRED_TO_DST
  #undef CHECK_FOR_BEST

  return bestMode;
}

void intra_recon(int16_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, int8_t mode, int8_t chroma)
{
  int32_t x,y,i;
  int16_t pred[LCU_WIDTH*LCU_WIDTH>>2];
  //int16_t* recShift = &rec[xpos+ypos*recstride];
  #define COPY_PRED_TO_DST() for(y = 0; y < (int32_t)width; y++)  {   for(x = 0; x < (int32_t)width; x++)  {  dst[x+y*dststride] = pred[x+y*width];  }   }

  /* Filtering apply if luma and not DC */
  if(!chroma && mode != 1/*&& width > 4*/)
  {
    uint8_t threshold = intraHorVerDistThres[g_toBits[width]];
    if(MIN(abs(mode-26),abs(mode-10)) > threshold)
    {
      intra_filter(rec,recstride,width,0);
    }
  }

  /* planar */  
  if(mode == 0)
  {
    intra_getPlanarPred(rec, recstride, xpos, ypos, width, pred, width); 
  }
  /* DC */
  else if(mode == 1)
  {
    i = intra_getDCPred(rec, recstride, xpos, ypos, width);
    for(y = 0; y < (int32_t)width; y++)
    {
      for(x = 0; x < (int32_t)width; x++)
      {
        dst[x+y*dststride] = i;
      }
    }
    return;
  }
  /* directional predictions */
  else
  {
    intra_getAngularPred(rec, recstride,pred, width,width,width,mode, xpos?1:0, ypos?1:0, 0);
  }

  COPY_PRED_TO_DST();

  #undef COPY_PRED_TO_DST
}


/*! \brief this functions build a reference block (only borders) used for intra predictions
    \param pic picture to use as a source, should contain full CU-data
    \param outwidth width of the prediction block
    \param chroma signaling if chroma is used, 0 = luma, 1 = U and 2 = V
    

*/
void intra_buildReferenceBorder(picture* pic, int32_t xCtb, int32_t yCtb,int8_t outwidth, int16_t* dst, int32_t dststride, int8_t chroma)
{
  int32_t leftColumn;  /*!< left column iterator */
  int16_t val;         /*!< variable to store extrapolated value */
  int32_t i;           /*!< index iterator */
  int16_t dcVal        = 1<<(g_uiBitDepth-1); /*!< default predictor value */
  int32_t topRow;      /*!< top row iterator */
  int32_t srcWidth     = (pic->width>>(chroma?1:0)); /*!< source picture width */
  int32_t srcHeight    = (pic->height>>(chroma?1:0));/*!< source picture height */
  uint8_t* srcPic      = (!chroma)?pic->yRecData: ((chroma==1)?pic->uRecData: pic->vRecData); /*!< input picture pointer */  
  int16_t SCU_width    = LCU_WIDTH>>(MAX_DEPTH+(chroma?1:0)); /*!< Smallest Coding Unit width */
  uint8_t* srcShifted  = &srcPic[xCtb*SCU_width+(yCtb*SCU_width)*srcWidth];  /*!< input picture pointer shifted to start from the left-top corner of the current block */
  int32_t width_in_SCU = srcWidth/SCU_width;     /*!< picture width in SCU */

  //memset(dst,127,outwidth*outwidth*sizeof(int16_t));

  /* Fill left column */
  if(xCtb)
  {
    /* Loop SCU's */
    for(leftColumn = 1; leftColumn < outwidth/SCU_width; leftColumn++)
    {
      /* If over the picture height or block not yet coded, stop */
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
  /*
  {
    FILE* test = fopen("blockout.yuv","wb");
    int x,y;
    uint8_t outvalue;
    for(y = 0; y < outwidth; y++)
    {
      for(x = 0; x < outwidth; x++)
      {
        outvalue = dst[x+y*outwidth];
        fwrite(&outvalue,1,1,test);
      }
    }
    fclose(test);
  }
  */
}

void intra_getAngularPred(int16_t* pSrc, int32_t srcStride, int16_t* rpDst, int32_t dstStride, int32_t width, int32_t height, int32_t dirMode, int8_t leftAvail,int8_t topAvail, int8_t filter)
{
  int32_t k,l;
  int32_t blkSize        = width;
  int16_t* pDst          = rpDst;

  // Map the mode index to main prediction direction and angle
  int8_t modeHor       = dirMode < 18;
  int8_t modeVer       = !modeHor;
  int32_t intraPredAngle = modeVer ? (int32_t)dirMode - 26 : modeHor ? -((int32_t)dirMode - 10) : 0;
  int32_t absAng         = abs(intraPredAngle);
  int32_t signAng        = intraPredAngle < 0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32_t angTable[9]    = {0,    2,    5,   9,  13,  17,  21,  26,  32};
  int32_t invAngTable[9] = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / Angle
  int32_t invAngle       = invAngTable[absAng];

  // Do angular predictions
  int16_t* refMain;
  int16_t* refSide;
  int16_t  refAbove[2*LCU_WIDTH+1];
  int16_t  refLeft[2*LCU_WIDTH+1];

  absAng             = angTable[absAng];
  intraPredAngle     = signAng * absAng;

  // Initialise the Main and Left reference array.
  if (intraPredAngle < 0)
  {
    int32_t invAngleSum = 128;       // rounding for (shift by 8)
    for (k=0;k<blkSize+1;k++)
    {
      refAbove[k+blkSize-1] = pSrc[k-srcStride-1];
    }
    for (k=0;k<blkSize+1;k++)
    {
      refLeft[k+blkSize-1] = pSrc[(k-1)*srcStride-1];
    }
    refMain = (modeVer ? refAbove : refLeft) + (blkSize-1);
    refSide = (modeVer ? refLeft : refAbove) + (blkSize-1);

    // Extend the Main reference to the left.    
    for (k=-1; k>blkSize*intraPredAngle>>5; k--)
    {
      invAngleSum += invAngle;
      refMain[k] = refSide[invAngleSum>>8];
    }
  }
  else
  {
    for (k=0;k<2*blkSize+1;k++)
    {
      refAbove[k] = pSrc[k-srcStride-1];
    }
    for (k=0;k<2*blkSize+1;k++)
    {
      refLeft[k] = pSrc[(k-1)*srcStride-1];
    }
    refMain = modeVer ? refAbove : refLeft;
    refSide = modeVer ? refLeft  : refAbove;
  }

  if (intraPredAngle == 0)
  {
    for (k=0;k<blkSize;k++)
    {
      for (l=0;l<blkSize;l++)
      {
        pDst[k*dstStride+l] = refMain[l+1];
      }
    }

    if(filter)
    {
      for (k=0;k<blkSize;k++)
      {
        pDst[k*dstStride] = CLIP(0, (1<<g_uiBitDepth)-1, pDst[k*dstStride] + (( refSide[k+1] - refSide[0] ) >> 1) );
      }
    }
  }
  else
  {
    int32_t deltaPos=0;
    int32_t deltaInt;
    int32_t deltaFract;
    int32_t refMainIndex;

    for (k=0;k<blkSize;k++)
    {
      deltaPos += intraPredAngle;
      deltaInt   = deltaPos >> 5;
      deltaFract = deltaPos & (32 - 1);

      if (deltaFract)
      {
        // Do linear filtering
        for (l=0;l<blkSize;l++)
        {
          refMainIndex        = l+deltaInt+1;
          pDst[k*dstStride+l] = (int16_t) ( ((32-deltaFract)*refMain[refMainIndex]+deltaFract*refMain[refMainIndex+1]+16) >> 5 );
        }
      }
      else
      {
        // Just copy the integer samples
        for (l=0;l<blkSize;l++)
        {
          pDst[k*dstStride+l] = refMain[l+deltaInt+1];
        }
      }
    }
  }

  // Flip the block if this is the horizontal mode
  if (modeHor)
  {
    int16_t tmp;
    for (k=0;k<blkSize-1;k++)
    {
      for (l=k+1;l<blkSize;l++)
      {
        tmp                 = pDst[k*dstStride+l];
        pDst[k*dstStride+l] = pDst[l*dstStride+k];
        pDst[l*dstStride+k] = tmp;
      }
    }
  }

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
void intra_getPlanarPred(int16_t* src,int32_t srcstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride)
{
  int16_t pDcVal = 1<<(g_uiBitDepth-1);
  int32_t k, l, bottomLeft, topRight;
  int32_t horPred;
  int32_t leftColumn[LCU_WIDTH], topRow[LCU_WIDTH], bottomRow[LCU_WIDTH], rightColumn[LCU_WIDTH];
  uint32_t blkSize = width;
  uint32_t offset2D = width;
  uint32_t shift1D = g_aucConvertToBit[ width ] + 2;
  uint32_t shift2D = shift1D + 1;


  for(k = 0; k < (int32_t)blkSize+1; k++)
  {
    topRow[k] = src[k-srcstride];
    leftColumn[k] = src[k*srcstride-1];
  }

  // Get left and above reference column and row
  

  // Prepare intermediate variables used in interpolation
  bottomLeft = leftColumn[blkSize];
  topRight   = topRow[blkSize];
  for (k = 0; k < (int32_t)blkSize; k++)
  {
    bottomRow[k]   = bottomLeft - topRow[k];
    rightColumn[k] = topRight   - leftColumn[k];
    topRow[k]      <<= shift1D;
    leftColumn[k]  <<= shift1D;
  }

  // Generate prediction signal
  for (k = 0; k < (int32_t)blkSize; k++)
  {
    horPred = leftColumn[k] + offset2D;
    for (l = 0; l < (int32_t)blkSize; l++)
    {
      horPred += rightColumn[k];
      topRow[l] += bottomRow[l];
      dst[k*dststride+l] = ( (horPred + topRow[l]) >> shift2D );
    }
  }
}
