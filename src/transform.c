/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file transform.c
    \brief Transform functions
    \author Marko Viitanen
    \date 2012-09
    
    Transform functions
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "nal.h"
#include "transform.h"


const int16_t g_aiT4[4][4] =
{
  { 64, 64, 64, 64},
  { 83, 36,-36,-83},
  { 64,-64,-64, 64},
  { 36,-83, 83,-36}
};

const int16_t g_aiT8[8][8] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64},
  { 89, 75, 50, 18,-18,-50,-75,-89},
  { 83, 36,-36,-83,-83,-36, 36, 83},
  { 75,-18,-89,-50, 50, 89, 18,-75},
  { 64,-64,-64, 64, 64,-64,-64, 64},
  { 50,-89, 18, 75,-75,-18, 89,-50},
  { 36,-83, 83,-36,-36, 83,-83, 36},
  { 18,-50, 75,-89, 89,-75, 50,-18}
};

const int16_t g_aiT16[16][16] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

const int16_t g_aiT32[32][32] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13,  4, -4,-13,-22,-31,-38,-46,-54,-61,-67,-73,-78,-82,-85,-88,-90,-90},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90,-90,-87,-80,-70,-57,-43,-25, -9,  9, 25, 43, 57, 70, 80, 87, 90},
  { 90, 82, 67, 46, 22, -4,-31,-54,-73,-85,-90,-88,-78,-61,-38,-13, 13, 38, 61, 78, 88, 90, 85, 73, 54, 31,  4,-22,-46,-67,-82,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89, 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 88, 67, 31,-13,-54,-82,-90,-78,-46, -4, 38, 73, 90, 85, 61, 22,-22,-61,-85,-90,-73,-38,  4, 46, 78, 90, 82, 54, 13,-31,-67,-88},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87,-87,-57, -9, 43, 80, 90, 70, 25,-25,-70,-90,-80,-43,  9, 57, 87},
  { 85, 46,-13,-67,-90,-73,-22, 38, 82, 88, 54, -4,-61,-90,-78,-31, 31, 78, 90, 61,  4,-54,-88,-82,-38, 22, 73, 90, 67, 13,-46,-85},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 82, 22,-54,-90,-61, 13, 78, 85, 31,-46,-90,-67,  4, 73, 88, 38,-38,-88,-73, -4, 67, 90, 46,-31,-85,-78,-13, 61, 90, 54,-22,-82},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80,-80, -9, 70, 87, 25,-57,-90,-43, 43, 90, 57,-25,-87,-70,  9, 80},
  { 78, -4,-82,-73, 13, 85, 67,-22,-88,-61, 31, 90, 54,-38,-90,-46, 46, 90, 38,-54,-90,-31, 61, 88, 22,-67,-85,-13, 73, 82,  4,-78},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75, 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 73,-31,-90,-22, 78, 67,-38,-90,-13, 82, 61,-46,-88, -4, 85, 54,-54,-85,  4, 88, 46,-61,-82, 13, 90, 38,-67,-78, 22, 90, 31,-73},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70,-70, 43, 87, -9,-90,-25, 80, 57,-57,-80, 25, 90,  9,-87,-43, 70},
  { 67,-54,-78, 38, 85,-22,-90,  4, 90, 13,-88,-31, 82, 46,-73,-61, 61, 73,-46,-82, 31, 88,-13,-90, -4, 90, 22,-85,-38, 78, 54,-67},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 61,-73,-46, 82, 31,-88,-13, 90, -4,-90, 22, 85,-38,-78, 54, 67,-67,-54, 78, 38,-85,-22, 90,  4,-90, 13, 88,-31,-82, 46, 73,-61},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57,-57, 80, 25,-90,  9, 87,-43,-70, 70, 43,-87, -9, 90,-25,-80, 57},
  { 54,-85, -4, 88,-46,-61, 82, 13,-90, 38, 67,-78,-22, 90,-31,-73, 73, 31,-90, 22, 78,-67,-38, 90,-13,-82, 61, 46,-88,  4, 85,-54},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50, 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 46,-90, 38, 54,-90, 31, 61,-88, 22, 67,-85, 13, 73,-82,  4, 78,-78, -4, 82,-73,-13, 85,-67,-22, 88,-61,-31, 90,-54,-38, 90,-46},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43,-43, 90,-57,-25, 87,-70, -9, 80,-80,  9, 70,-87, 25, 57,-90, 43},
  { 38,-88, 73, -4,-67, 90,-46,-31, 85,-78, 13, 61,-90, 54, 22,-82, 82,-22,-54, 90,-61,-13, 78,-85, 31, 46,-90, 67,  4,-73, 88,-38},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 31,-78, 90,-61,  4, 54,-88, 82,-38,-22, 73,-90, 67,-13,-46, 85,-85, 46, 13,-67, 90,-73, 22, 38,-82, 88,-54, -4, 61,-90, 78,-31},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25,-25, 70,-90, 80,-43, -9, 57,-87, 87,-57,  9, 43,-80, 90,-70, 25},
  { 22,-61, 85,-90, 73,-38, -4, 46,-78, 90,-82, 54,-13,-31, 67,-88, 88,-67, 31, 13,-54, 82,-90, 78,-46,  4, 38,-73, 90,-85, 61,-22},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18, 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  { 13,-38, 61,-78, 88,-90, 85,-73, 54,-31,  4, 22,-46, 67,-82, 90,-90, 82,-67, 46,-22, -4, 31,-54, 73,-85, 90,-88, 78,-61, 38,-13},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9, -9, 25,-43, 57,-70, 80,-87, 90,-90, 87,-80, 70,-57, 43,-25,  9},
  {  4,-13, 22,-31, 38,-46, 54,-61, 67,-73, 78,-82, 85,-88, 90,-90, 90,-90, 88,-85, 82,-78, 73,-67, 61,-54, 46,-38, 31,-22, 13, -4}
};

const int32_t g_quantTSDefault4x4[16] =
{
  16,16,16,16,
  16,16,16,16,
  16,16,16,16,
  16,16,16,16
};

const int32_t g_quantIntraDefault8x8[64] =
{
  16,16,16,16,17,18,21,24,
  16,16,16,16,17,19,22,25,
  16,16,17,18,20,22,25,29,
  16,16,18,21,24,27,31,36,
  17,17,20,24,30,35,41,47,
  18,19,22,27,35,44,54,65,
  21,22,25,31,41,54,70,88,
  24,25,29,36,47,65,88,115
};

const int32_t g_quantInterDefault8x8[64] =
{
  16,16,16,16,17,18,20,24,
  16,16,16,17,18,20,24,25,
  16,16,17,18,20,24,25,28,
  16,17,18,20,24,25,28,33,
  17,18,20,24,25,28,33,41,
  18,20,24,25,28,33,41,54,
  20,24,25,28,33,41,54,71,
  24,25,28,33,41,54,71,91
};

const uint8_t g_aucChromaScale[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};

int32_t* g_quant_coeff[4][6][6];
int32_t* g_de_quant_coeff[4][6][6];

const  uint8_t g_scalingListNum[4]={6,6,6,2};
const  uint16_t g_scalingListSize[4] = {16,64,256,1024}; 
const  uint8_t g_scalingListSizeX[4] = { 4, 8, 16,  32};
const int16_t g_quantScales[6] = { 26214,23302,20560,18396,16384,14564 };
const int16_t g_invQuantScales[6] = { 40,45,51,57,64,72 };
//static int32_t  m_scalingListDC[4][6];

void scalinglist_init()
{
  uint32_t sizeId,listId,qp;
  for(sizeId = 0; sizeId < 4; sizeId++)
  {
    for(listId = 0; listId < g_scalingListNum[sizeId]; listId++)
    {
      for(qp = 0; qp < 6; qp++)
      {
        if(!(sizeId == 3 && listId == 3))
        {
          g_quant_coeff   [sizeId][listId][qp] = (int32_t*)malloc(sizeof(int32_t)*g_scalingListSize[sizeId]);
          memset(g_quant_coeff[sizeId][listId][qp],0,sizeof(int32_t)*g_scalingListSize[sizeId]);
          g_de_quant_coeff   [sizeId][listId][qp] = (int32_t*)malloc(sizeof(int32_t)*g_scalingListSize[sizeId]);
          memset(g_de_quant_coeff[sizeId][listId][qp],0,sizeof(int32_t)*g_scalingListSize[sizeId]);
        
        //m_dequantCoef [sizeId][listId][qp][SCALING_LIST_SQT] = new Int [g_scalingListSize[sizeId]];
        //m_errScale    [sizeId][listId][qp][SCALING_LIST_SQT] = new double [g_scalingListSize[sizeId]];
        }

      }
    }
  }
  //Alias
  for(qp = 0; qp < 6; qp++)
  {
    g_quant_coeff[3][3][qp] = g_quant_coeff[3][1][qp];
    g_de_quant_coeff[3][3][qp] = g_de_quant_coeff[3][1][qp];
  }
}

void scalinglist_destroy()
{
  uint32_t sizeId,listId,qp;
  for(sizeId = 0; sizeId < 4; sizeId++)
  {
    for(listId = 0; listId < g_scalingListNum[sizeId]; listId++)
    {
      for(qp = 0; qp < 6; qp++)
      {
         free(g_quant_coeff[sizeId][listId][qp]);
         free(g_de_quant_coeff[sizeId][listId][qp]);
      }
    }
  }
}


void scalinglist_process()
{
  uint32_t size,list,qp;
  for(size=0;size</*SCALING_LIST_SIZE_NUM*/4;size++)
  {
    int32_t* list_ptr = (int32_t *)g_quantIntraDefault8x8; /* Default to "8x8" intra */

    for(list = 0; list < g_scalingListNum[size]; list++)
    {
      switch(size)
      {
        case 0: /* 4x4 */
          list_ptr = (int32_t *)g_quantTSDefault4x4;
          break;
        case 1: /* 8x8 */
        case 2: /* 16x16 */
          if(list > 2) list_ptr = (int32_t*)g_quantInterDefault8x8;
          break;
        case 3: /* 32x32 */
          if(list > 0) list_ptr = (int32_t*)g_quantInterDefault8x8;
          break;
      }

      for(qp=0;qp</*SCALING_LIST_REM_NUM*/6;qp++)
      {
        scalinglist_set(list_ptr,list,size,qp);
      }
    }
  }
}

void scalinglist_processEnc( int32_t *coeff, int32_t *quantcoeff, int32_t quantScales, uint32_t height,uint32_t width, uint32_t ratio, int32_t sizuNum, uint32_t dc, uint8_t flat)
{
  uint32_t j,i;
  int32_t nsqth = (height < width) ? 4: 1; //height ratio for NSQT
  int32_t nsqtw = (width < height) ? 4: 1; //width ratio for NSQT 
  if(flat)
  {
    for(j=0;j<height*width;j++)
    {
      *quantcoeff++ = quantScales>>4;
    }
  }
  else
  {
    for(j=0;j<height;j++)
    {
      for(i=0;i<width;i++)
      {
        uint32_t coeffpos  = sizuNum * (j * nsqth / ratio) + i * nsqtw /ratio;
        quantcoeff[j*width + i] = quantScales / ((coeffpos>63)?1:coeff[coeffpos]);
      }
    }
    if(ratio > 1)
    {
      quantcoeff[0] = quantScales / dc;
    }
  }
}

void scalinglist_processDec( int32_t *coeff, int32_t *dequantcoeff, int32_t invQuantScales, uint32_t height,uint32_t width, uint32_t ratio, int32_t sizuNum, uint32_t dc, uint8_t flat)
{
  uint32_t j,i;
  for(j=0;j<height;j++)
  {
    for(i=0;i<width;i++)
    {
      dequantcoeff[j*width + i] = invQuantScales * coeff[sizuNum * (j / ratio) + i / ratio];
    }
  }
  if(ratio > 1)
  {
    dequantcoeff[0] = invQuantScales * dc;
  }
}

void scalinglist_set(int32_t *coeff, uint32_t listId, uint32_t sizeId, uint32_t qp)
{
  uint32_t width  = g_scalingListSizeX[sizeId];
  uint32_t height = g_scalingListSizeX[sizeId];
  uint32_t ratio  = g_scalingListSizeX[sizeId]/MIN(8,g_scalingListSizeX[sizeId]);
  int32_t *quantcoeff   = g_quant_coeff[sizeId][listId][qp];
  int32_t *dequantcoeff = g_de_quant_coeff[sizeId][listId][qp];

  scalinglist_processEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*SCALING_LIST_DC*/16, 0);

  scalinglist_processDec(coeff,dequantcoeff,g_invQuantScales[qp],height,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*SCALING_LIST_DC*/16, 0);

  //ToDo: support NSQT
  //if(sizeId == /*SCALING_LIST_32x32*/3 || sizeId == /*SCALING_LIST_16x16*/2) //for NSQT
  //{
  //  quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_VER*/1];
  //  scalinglist_processEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width>>2,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);

  //  quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_HOR*/2];
  //  scalinglist_processEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height>>2,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);
  //}
}


void partialButterfly4(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j;  
  int32_t E[2],O[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {    
    /* E and O */
    E[0] = src[0] + src[3];
    O[0] = src[0] - src[3];
    E[1] = src[1] + src[2];
    O[1] = src[1] - src[2];

    dst[0] = (g_aiT4[0][0]*E[0] + g_aiT4[0][1]*E[1] + add)>>shift;
    dst[2*line] = (g_aiT4[2][0]*E[0] + g_aiT4[2][1]*E[1] + add)>>shift;
    dst[line] = (g_aiT4[1][0]*O[0] + g_aiT4[1][1]*O[1] + add)>>shift;
    dst[3*line] = (g_aiT4[3][0]*O[0] + g_aiT4[3][1]*O[1] + add)>>shift;

    src += 4;
    dst ++;
  }
}

void partialButterflyInverse4(short *src,short *dst,int shift, int line)
{
  int j;
  int E[2],O[2];
  int add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */    
    O[0] = g_aiT4[1][0]*src[line] + g_aiT4[3][0]*src[3*line];
    O[1] = g_aiT4[1][1]*src[line] + g_aiT4[3][1]*src[3*line];
    E[0] = g_aiT4[0][0]*src[0] + g_aiT4[2][0]*src[2*line];
    E[1] = g_aiT4[0][1]*src[0] + g_aiT4[2][1]*src[2*line];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    dst[0] = CLIP( -32768, 32767, (E[0] + O[0] + add)>>shift );
    dst[1] = CLIP( -32768, 32767, (E[1] + O[1] + add)>>shift );
    dst[2] = CLIP( -32768, 32767, (E[1] - O[1] + add)>>shift );
    dst[3] = CLIP( -32768, 32767, (E[0] - O[0] + add)>>shift );
            
    src   ++;
    dst += 4;
  }
}

// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm 
// give identical results
void fastForwardDst(short *block,short *coeff,int32_t shift)  // input block, output coeff
{
  int32_t i, c[4];
  int32_t rnd_factor = 1<<(shift-1);
  for (i=0; i<4; i++)
  {
    // int32_termediate Variables
    c[0] = block[4*i+0] + block[4*i+3];
    c[1] = block[4*i+1] + block[4*i+3];
    c[2] = block[4*i+0] - block[4*i+1];
    c[3] = 74* block[4*i+2];

    coeff[   i] =  ( 29 * c[0] + 55 * c[1]         + c[3]               + rnd_factor ) >> shift;
    coeff[ 4+i] =  ( 74 * (block[4*i+0]+ block[4*i+1] - block[4*i+3])   + rnd_factor ) >> shift;
    coeff[ 8+i] =  ( 29 * c[2] + 55 * c[0]         - c[3]               + rnd_factor ) >> shift;
    coeff[12+i] =  ( 55 * c[2] - 29 * c[1]         + c[3]               + rnd_factor ) >> shift;
  }
}

void fastInverseDst(short *tmp,short *block,int shift)  // input tmp, output block
{
  int i, c[4];
  int rnd_factor = 1<<(shift-1);
  for (i=0; i<4; i++)
  {  
    // Intermediate Variables
    c[0] = tmp[  i] + tmp[ 8+i];
    c[1] = tmp[8+i] + tmp[12+i];
    c[2] = tmp[  i] - tmp[12+i];
    c[3] = 74* tmp[4+i];

    block[4*i+0] = CLIP( -32768, 32767, ( 29 * c[0] + 55 * c[1]     + c[3]               + rnd_factor ) >> shift );
    block[4*i+1] = CLIP( -32768, 32767, ( 55 * c[2] - 29 * c[1]     + c[3]               + rnd_factor ) >> shift );
    block[4*i+2] = CLIP( -32768, 32767, ( 74 * (tmp[i] - tmp[8+i]  + tmp[12+i])      + rnd_factor ) >> shift );
    block[4*i+3] = CLIP( -32768, 32767, ( 55 * c[0] + 29 * c[2]     - c[3]               + rnd_factor ) >> shift );
  }
}


void partialButterfly8(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;  
  int32_t E[4],O[4];
  int32_t EE[2],EO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {  
    /* E and O*/
    for (k=0;k<4;k++)
    {
      E[k] = src[k] + src[7-k];
      O[k] = src[k] - src[7-k];
    }    
    /* EE and EO */
    EE[0] = E[0] + E[3];    
    EO[0] = E[0] - E[3];
    EE[1] = E[1] + E[2];
    EO[1] = E[1] - E[2];

    dst[0] = (g_aiT8[0][0]*EE[0] + g_aiT8[0][1]*EE[1] + add)>>shift;
    dst[4*line] = (g_aiT8[4][0]*EE[0] + g_aiT8[4][1]*EE[1] + add)>>shift; 
    dst[2*line] = (g_aiT8[2][0]*EO[0] + g_aiT8[2][1]*EO[1] + add)>>shift;
    dst[6*line] = (g_aiT8[6][0]*EO[0] + g_aiT8[6][1]*EO[1] + add)>>shift; 

    dst[line] = (g_aiT8[1][0]*O[0] + g_aiT8[1][1]*O[1] + g_aiT8[1][2]*O[2] + g_aiT8[1][3]*O[3] + add)>>shift;
    dst[3*line] = (g_aiT8[3][0]*O[0] + g_aiT8[3][1]*O[1] + g_aiT8[3][2]*O[2] + g_aiT8[3][3]*O[3] + add)>>shift;
    dst[5*line] = (g_aiT8[5][0]*O[0] + g_aiT8[5][1]*O[1] + g_aiT8[5][2]*O[2] + g_aiT8[5][3]*O[3] + add)>>shift;
    dst[7*line] = (g_aiT8[7][0]*O[0] + g_aiT8[7][1]*O[1] + g_aiT8[7][2]*O[2] + g_aiT8[7][3]*O[3] + add)>>shift;

    src += 8;
    dst ++;
  }
}

void partialButterflyInverse8(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t E[4],O[4];
  int32_t EE[2],EO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++) 
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<4;k++)
    {
      O[k] = g_aiT8[ 1][k]*src[line] + g_aiT8[ 3][k]*src[3*line] + g_aiT8[ 5][k]*src[5*line] + g_aiT8[ 7][k]*src[7*line];
    }

    EO[0] = g_aiT8[2][0]*src[ 2*line ] + g_aiT8[6][0]*src[ 6*line ];
    EO[1] = g_aiT8[2][1]*src[ 2*line ] + g_aiT8[6][1]*src[ 6*line ];
    EE[0] = g_aiT8[0][0]*src[ 0      ] + g_aiT8[4][0]*src[ 4*line ];
    EE[1] = g_aiT8[0][1]*src[ 0      ] + g_aiT8[4][1]*src[ 4*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    E[0] = EE[0] + EO[0];
    E[3] = EE[0] - EO[0];
    E[1] = EE[1] + EO[1];
    E[2] = EE[1] - EO[1];
    for (k=0;k<4;k++)
    {
      dst[ k   ] = MAX( -32768, MIN(32767, (E[k] + O[k] + add)>>shift ));
      dst[ k+4 ] = MAX( -32768, MIN(32767, (E[3-k] - O[3-k] + add)>>shift ));
    }   
    src ++;
    dst += 8;
  }
}



void partialButterfly16(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t E[8],O[8];
  int32_t EE[4],EO[4];
  int32_t EEE[2],EEO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++) 
  {    
    /* E and O*/
    for (k=0;k<8;k++)
    {
      E[k] = src[k] + src[15-k];
      O[k] = src[k] - src[15-k];
    } 
    /* EE and EO */
    for (k=0;k<4;k++)
    {
      EE[k] = E[k] + E[7-k];
      EO[k] = E[k] - E[7-k];
    }
    /* EEE and EEO */
    EEE[0] = EE[0] + EE[3];    
    EEO[0] = EE[0] - EE[3];
    EEE[1] = EE[1] + EE[2];
    EEO[1] = EE[1] - EE[2];

    dst[ 0      ] = (g_aiT16[ 0][0]*EEE[0] + g_aiT16[ 0][1]*EEE[1] + add)>>shift;        
    dst[ 8*line ] = (g_aiT16[ 8][0]*EEE[0] + g_aiT16[ 8][1]*EEE[1] + add)>>shift;    
    dst[ 4*line ] = (g_aiT16[ 4][0]*EEO[0] + g_aiT16[ 4][1]*EEO[1] + add)>>shift;        
    dst[ 12*line] = (g_aiT16[12][0]*EEO[0] + g_aiT16[12][1]*EEO[1] + add)>>shift;

    for (k=2;k<16;k+=4)
    {
      dst[ k*line ] = (g_aiT16[k][0]*EO[0] + g_aiT16[k][1]*EO[1] + g_aiT16[k][2]*EO[2] + g_aiT16[k][3]*EO[3] + add)>>shift;      
    }

    for (k=1;k<16;k+=2)
    {
      dst[ k*line ] = (g_aiT16[k][0]*O[0] + g_aiT16[k][1]*O[1] + g_aiT16[k][2]*O[2] + g_aiT16[k][3]*O[3] + 
        g_aiT16[k][4]*O[4] + g_aiT16[k][5]*O[5] + g_aiT16[k][6]*O[6] + g_aiT16[k][7]*O[7] + add)>>shift;
    }

    src += 16;
    dst ++; 

  }
}


void partialButterflyInverse16(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t E[8],O[8];
  int32_t EE[4],EO[4];
  int32_t EEE[2],EEO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {    
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<8;k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ line] + g_aiT16[ 3][k]*src[ 3*line] + g_aiT16[ 5][k]*src[ 5*line] + g_aiT16[ 7][k]*src[ 7*line] + 
        g_aiT16[ 9][k]*src[ 9*line] + g_aiT16[11][k]*src[11*line] + g_aiT16[13][k]*src[13*line] + g_aiT16[15][k]*src[15*line];
    }
    for (k=0;k<4;k++)
    {
      EO[k] = g_aiT16[ 2][k]*src[ 2*line] + g_aiT16[ 6][k]*src[ 6*line] + g_aiT16[10][k]*src[10*line] + g_aiT16[14][k]*src[14*line];
    }
    EEO[0] = g_aiT16[4][0]*src[ 4*line ] + g_aiT16[12][0]*src[ 12*line ];
    EEE[0] = g_aiT16[0][0]*src[ 0      ] + g_aiT16[ 8][0]*src[ 8*line  ];
    EEO[1] = g_aiT16[4][1]*src[ 4*line ] + g_aiT16[12][1]*src[ 12*line ];
    EEE[1] = g_aiT16[0][1]*src[ 0      ] + g_aiT16[ 8][1]*src[ 8*line  ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */ 
    for (k=0;k<2;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }    
    for (k=0;k<4;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      
      dst[k]   = MAX( -32768, MIN(32767, (E[k] + O[k] + add)>>shift));
      dst[k+8] = MAX( -32768, MIN(32767, (E[7-k] - O[7-k] + add)>>shift));
    }   
    src ++; 
    dst += 16;
  }
}



void partialButterfly32(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t E[16],O[16];
  int32_t EE[8],EO[8];
  int32_t EEE[4],EEO[4];
  int32_t EEEE[2],EEEO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {    
    /* E and O*/
    for (k=0;k<16;k++)
    {
      E[k] = src[k] + src[31-k];
      O[k] = src[k] - src[31-k];
    } 
    /* EE and EO */
    for (k=0;k<8;k++)
    {
      EE[k] = E[k] + E[15-k];
      EO[k] = E[k] - E[15-k];
    }
    /* EEE and EEO */
    for (k=0;k<4;k++)
    {
      EEE[k] = EE[k] + EE[7-k];
      EEO[k] = EE[k] - EE[7-k];
    }
    /* EEEE and EEEO */
    EEEE[0] = EEE[0] + EEE[3];    
    EEEO[0] = EEE[0] - EEE[3];
    EEEE[1] = EEE[1] + EEE[2];
    EEEO[1] = EEE[1] - EEE[2];

    dst[ 0       ] = (g_aiT32[ 0][0]*EEEE[0] + g_aiT32[ 0][1]*EEEE[1] + add)>>shift;
    dst[ 16*line ] = (g_aiT32[16][0]*EEEE[0] + g_aiT32[16][1]*EEEE[1] + add)>>shift;
    dst[ 8*line  ] = (g_aiT32[ 8][0]*EEEO[0] + g_aiT32[ 8][1]*EEEO[1] + add)>>shift; 
    dst[ 24*line ] = (g_aiT32[24][0]*EEEO[0] + g_aiT32[24][1]*EEEO[1] + add)>>shift;
    for (k=4;k<32;k+=8)
    {
      dst[ k*line ] = (g_aiT32[k][0]*EEO[0] + g_aiT32[k][1]*EEO[1] + g_aiT32[k][2]*EEO[2] + g_aiT32[k][3]*EEO[3] + add)>>shift;
    }       
    for (k=2;k<32;k+=4)
    {
      dst[ k*line ] = (g_aiT32[k][0]*EO[0] + g_aiT32[k][1]*EO[1] + g_aiT32[k][2]*EO[2] + g_aiT32[k][3]*EO[3] + 
        g_aiT32[k][4]*EO[4] + g_aiT32[k][5]*EO[5] + g_aiT32[k][6]*EO[6] + g_aiT32[k][7]*EO[7] + add)>>shift;
    }       
    for (k=1;k<32;k+=2)
    {
      dst[ k*line ] = (g_aiT32[k][ 0]*O[ 0] + g_aiT32[k][ 1]*O[ 1] + g_aiT32[k][ 2]*O[ 2] + g_aiT32[k][ 3]*O[ 3] + 
        g_aiT32[k][ 4]*O[ 4] + g_aiT32[k][ 5]*O[ 5] + g_aiT32[k][ 6]*O[ 6] + g_aiT32[k][ 7]*O[ 7] +
        g_aiT32[k][ 8]*O[ 8] + g_aiT32[k][ 9]*O[ 9] + g_aiT32[k][10]*O[10] + g_aiT32[k][11]*O[11] + 
        g_aiT32[k][12]*O[12] + g_aiT32[k][13]*O[13] + g_aiT32[k][14]*O[14] + g_aiT32[k][15]*O[15] + add)>>shift;
    }
    src += 32;
    dst ++;
  }
}


void partialButterflyInverse32(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t E[16],O[16];
  int32_t EE[8],EO[8];
  int32_t EEE[4],EEO[4];
  int32_t EEEE[2],EEEO[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (k=0;k<16;k++)
    {
      O[k] = g_aiT32[ 1][k]*src[ line  ] + g_aiT32[ 3][k]*src[ 3*line  ] + g_aiT32[ 5][k]*src[ 5*line  ] + g_aiT32[ 7][k]*src[ 7*line  ] + 
        g_aiT32[ 9][k]*src[ 9*line  ] + g_aiT32[11][k]*src[ 11*line ] + g_aiT32[13][k]*src[ 13*line ] + g_aiT32[15][k]*src[ 15*line ] + 
        g_aiT32[17][k]*src[ 17*line ] + g_aiT32[19][k]*src[ 19*line ] + g_aiT32[21][k]*src[ 21*line ] + g_aiT32[23][k]*src[ 23*line ] + 
        g_aiT32[25][k]*src[ 25*line ] + g_aiT32[27][k]*src[ 27*line ] + g_aiT32[29][k]*src[ 29*line ] + g_aiT32[31][k]*src[ 31*line ];
    }
    for (k=0;k<8;k++)
    {
      EO[k] = g_aiT32[ 2][k]*src[ 2*line  ] + g_aiT32[ 6][k]*src[ 6*line  ] + g_aiT32[10][k]*src[ 10*line ] + g_aiT32[14][k]*src[ 14*line ] + 
        g_aiT32[18][k]*src[ 18*line ] + g_aiT32[22][k]*src[ 22*line ] + g_aiT32[26][k]*src[ 26*line ] + g_aiT32[30][k]*src[ 30*line ];
    }
    for (k=0;k<4;k++)
    {
      EEO[k] = g_aiT32[4][k]*src[ 4*line ] + g_aiT32[12][k]*src[ 12*line ] + g_aiT32[20][k]*src[ 20*line ] + g_aiT32[28][k]*src[ 28*line ];
    }
    EEEO[0] = g_aiT32[8][0]*src[ 8*line ] + g_aiT32[24][0]*src[ 24*line ];
    EEEO[1] = g_aiT32[8][1]*src[ 8*line ] + g_aiT32[24][1]*src[ 24*line ];
    EEEE[0] = g_aiT32[0][0]*src[ 0      ] + g_aiT32[16][0]*src[ 16*line ];    
    EEEE[1] = g_aiT32[0][1]*src[ 0      ] + g_aiT32[16][1]*src[ 16*line ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    EEE[0] = EEEE[0] + EEEO[0];
    EEE[3] = EEEE[0] - EEEO[0];
    EEE[1] = EEEE[1] + EEEO[1];
    EEE[2] = EEEE[1] - EEEO[1];    
    for (k=0;k<4;k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+4] = EEE[3-k] - EEO[3-k];
    }    
    for (k=0;k<8;k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+8] = EE[7-k] - EO[7-k];
    }    
    for (k=0;k<16;k++)
    {
      dst[k]    = MAX( -32768, MIN(32767, (E[k] + O[k] + add)>>shift ));
      dst[k+16] = MAX( -32768, MIN(32767, (E[15-k] - O[15-k] + add)>>shift ));
    }
    src ++;
    dst += 32;
  }
}


/** forward transform (2D)
*  \param block input residual
*  \param coeff transform coefficients
*  \param blockSize width of transform
*/
void transform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int32_t uiMode)
{

  int32_t shift_1st = g_aucConvertToBit[blockSize]  + 1 + g_uiBitIncrement; // log2(iWidth) - 1 + g_uiBitIncrement
  int32_t shift_2nd = g_aucConvertToBit[blockSize]  + 8;                   // log2(iHeight) + 6

  int16_t tmp[LCU_WIDTH*LCU_WIDTH];
  
  if(blockSize== 4)
  {
    if (uiMode != 65535)
    {
      fastForwardDst(block,tmp,shift_1st); // Forward DST BY FAST ALGORITHM, block input, tmp output
      fastForwardDst(tmp,coeff,shift_2nd); // Forward DST BY FAST ALGORITHM, tmp input, coeff output
    }
    else
    {
      partialButterfly4(block, tmp, shift_1st, blockSize);
      partialButterfly4(tmp, coeff, shift_2nd, blockSize);
    }

  }
  else
  {
    switch(blockSize)
    {
      case 8:
      {
        partialButterfly8( block, tmp, shift_1st, blockSize );
        partialButterfly8( tmp, coeff, shift_2nd, blockSize );
        break;
      }
      case 16:
      {
        partialButterfly16( block, tmp, shift_1st, blockSize );
        partialButterfly16( tmp, coeff, shift_2nd, blockSize );
        break;
      }
      case 32:
      {
        partialButterfly32( block, tmp, shift_1st, blockSize );
        partialButterfly32( tmp, coeff, shift_2nd, blockSize );
        break;
      }
    }
  }
}

/*! \brief NxN inverse transform (2D)
 \param coeff input data (transform coefficients)
 \param block output data (residual)
 \param blockSize input data (width of transform)
 \param uiMode
*/
void itransform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int32_t uiMode)
  //(Int bitDepth, Short *coeff,Short *block, Int iWidth, Int iHeight, UInt uiMode)
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (g_bitDepth-8);
  int16_t tmp[LCU_WIDTH*LCU_WIDTH];

  if( blockSize == 4)
  {
    if (uiMode != 65535)
    {
      fastInverseDst(coeff,tmp,shift_1st);    // Inverse DST by FAST Algorithm, coeff input, tmp output
      fastInverseDst(tmp,block,shift_2nd); // Inverse DST by FAST Algorithm, tmp input, coeff output
    }
    else
    {
      partialButterflyInverse4(coeff,tmp,shift_1st,blockSize);
      partialButterflyInverse4(tmp,block,shift_2nd,blockSize);
    }
  }
  else if( blockSize == 8)
  {
    partialButterflyInverse8(coeff,tmp,shift_1st,blockSize);
    partialButterflyInverse8(tmp,block,shift_2nd,blockSize);
  }
  else if( blockSize == 16)
  {
    partialButterflyInverse16(coeff,tmp,shift_1st,blockSize);
    partialButterflyInverse16(tmp,block,shift_2nd,blockSize);
  }
  else if( blockSize == 32)
  {
    partialButterflyInverse32(coeff,tmp,shift_1st,blockSize);
    partialButterflyInverse32(tmp,block,shift_2nd,blockSize);
  }
}


#define QUANT_SHIFT 14

void quant(encoder_control* encoder, int16_t* pSrc, int16_t* pDes, int32_t iWidth,
           int32_t iHeight, uint32_t *uiAcSum, int8_t eTType, int8_t scanIdx )
{
  int16_t*   piCoef    = pSrc;
  int16_t*   piQCoef   = pDes;
  uint32_t*  scan;
 
  int8_t useRDOQForTransformSkip = 0;
  uint32_t log2BlockSize = g_aucConvertToBit[ iWidth ] + 2;

  //uint32_t scanIdx = SCAN_DIAG;

  scan = g_auiSigLastScan[ scanIdx ][ log2BlockSize - 1 ];
  {
  int32_t deltaU[LCU_WIDTH*LCU_WIDTH] ;
  int32_t iQpBase = encoder->QP;

  int32_t qpScaled;
  int32_t qpBDOffset = 0;//(eTType == 0)? pcCU->getSlice()->getSPS()->getQpBDOffsetY() : pcCU->getSlice()->getSPS()->getQpBDOffsetC();

  if(eTType == 0)
  {
    qpScaled = iQpBase + qpBDOffset;
  }
  else
  {
    qpScaled = MAX( -qpBDOffset, MIN(57, iQpBase));
    if(qpScaled < 0)
    {
      qpScaled = qpScaled +  qpBDOffset;
    }
    else
    {
      qpScaled = g_aucChromaScale[ qpScaled ] + qpBDOffset;
    }
  }
  
  //New block for variable definitions
  {
  int32_t n;
  uint32_t dir = 0;//SCALING_LIST_SQT;
    
  uint32_t uiLog2TrSize = g_aucConvertToBit[ iWidth ] + 2;
  int32_t scalingListType = (/*pcCU->isint32_tra(uiAbsPartIdx)*/0 ? 0 : 3) + (int8_t)("\0\3\1\2"[eTType]);
  
  int32_t *piQuantCoeff = g_quant_coeff[uiLog2TrSize-2][scalingListType][/*m_cQP.m_iRem*/qpScaled%6];

  uint32_t uiBitDepth = g_bitDepth;

  int32_t iTransformShift = /*MAX_TR_DYNAMIC_RANGE*/15 - uiBitDepth - uiLog2TrSize; // Represents scaling through forward transform
  int32_t iQBits = QUANT_SHIFT + /*cQpBase.m_iPer +*/qpScaled/6 + iTransformShift;
  int32_t iAdd = ((encoder->in.cur_pic.slicetype == SLICE_I) ? 171 : 85) << (iQBits-9);

  int32_t qBits8 = iQBits-8;
  for( n = 0; n < iWidth*iHeight; n++ )
  {
    int32_t iLevel;
    int32_t  iSign;
    int64_t tmpLevel;
    iLevel  = piCoef[n];
    iSign   = (iLevel < 0 ? -1: 1);

    tmpLevel  = (int64_t)abs(iLevel) * piQuantCoeff[n];
    iLevel    = (int32_t)((tmpLevel + iAdd ) >> iQBits);
    deltaU[n] = (int32_t)((tmpLevel - (iLevel<<iQBits) )>> qBits8);

    iLevel *= iSign;        
    piQCoef[n] = CLIP( -32768, 32767, iLevel);
  } // for n
  }
  }

}


void dequant(encoder_control* encoder, int16_t* piQCoef, int16_t* piCoef, int32_t iWidth, int32_t iHeight,int8_t eTType)
{
  int32_t iShift,iAdd,iCoeffQ;
  uint32_t uiLog2TrSize = g_aucConvertToBit[ iWidth ] + 2;
  int16_t clipQCoef;
  int32_t n;
  int32_t iTransformShift = 15 - g_bitDepth - (g_aucConvertToBit[ iWidth ] + 2);
  int32_t qpScaled;
  int32_t iQpBase = encoder->QP;
  int32_t scalingListType = (/*pcCU->isintra(uiAbsPartIdx)*/1 ? 0 : 3) + (int8_t)("\0\3\1\2"[eTType]);
  int32_t *piDequantCoef;

  if(eTType == 0)
  {
    qpScaled = iQpBase;
  }
  else
  {
    qpScaled = CLIP( 0, 57, iQpBase);
    if(qpScaled < 0)
    {
      qpScaled = qpScaled;
    }
    else
    {
      qpScaled = g_aucChromaScale[ qpScaled ];
    }
  }
  piDequantCoef = g_de_quant_coeff[uiLog2TrSize-2][scalingListType][qpScaled%6];

  iShift = 20 - 14 - iTransformShift;

  iShift += 4;
  

  if(iShift >qpScaled/6)
  {
    iAdd = 1 << (iShift - qpScaled/6 - 1);
      
    for(n = 0; n < iWidth*iHeight; n++ )
    {
      clipQCoef = CLIP( -32768, 32767, piQCoef[n] );
      iCoeffQ = ((clipQCoef * piDequantCoef[n]) + iAdd ) >> (iShift -  qpScaled/6);
      piCoef[n] = CLIP(-32768,32767,iCoeffQ);
    }
  }
  else
  {
    for(n = 0; n < iWidth*iHeight; n++ )
    {
      clipQCoef = CLIP( -32768, 32767, piQCoef[n] );
      iCoeffQ   = CLIP( -32768, 32767, clipQCoef * piDequantCoef[n] ); // Clip to avoid possible overflow in following shift left operation
      piCoef[n] = CLIP( -32768, 32767, iCoeffQ << ( qpScaled/6 - iShift ) );
    }
  }
  /*
  iAdd = 1 << (iShift-1);
  scale = g_invQuantScales[encoder->QP%6] << (encoder->QP/6);

  for(n = 0; n < iWidth*iHeight; n++ )
  {
    clipQCoef = CLIP( -32768, 32767, piQCoef[n]);
    iCoeffQ   = ( clipQCoef * scale + iAdd ) >> iShift;
    piCoef[n] = CLIP( -32768, 32767, iCoeffQ);
  }
  */
}
