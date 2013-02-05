/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
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

int32_t* g_quant_coeff[4][6][6][3];

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

uint8_t g_scalingListNum[4]={6,6,6,2};
uint16_t g_scalingListSize[4] = {16,64,256,1024}; 
uint8_t g_scalingListSizeX[4] = { 4, 8, 16,  32};
int16_t g_quantScales[6] = { 26214,23302,20560,18396,16384,14564 };    

void scalinglist_init()
{
  uint32_t sizeId,listId,qp,dir;
  for(sizeId = 0; sizeId < 4; sizeId++)
  {
    for(listId = 0; listId < g_scalingListNum[sizeId]; listId++)
    {
      for(qp = 0; qp < 6; qp++)
      {
         g_quant_coeff   [sizeId][listId][qp][0] = (int32_t*)malloc(g_scalingListSize[sizeId]);
        //m_dequantCoef [sizeId][listId][qp][SCALING_LIST_SQT] = new Int [g_scalingListSize[sizeId]];
        //m_errScale    [sizeId][listId][qp][SCALING_LIST_SQT] = new double [g_scalingListSize[sizeId]];
        
        if(sizeId == /*SCALING_LIST_8x8*/1 || (sizeId == /*SCALING_LIST_16x16*/2 && listId < 2))
        {
          for(dir = /*SCALING_LIST_VER*/1; dir < /*SCALING_LIST_DIR_NUM*/3; dir++)
          {
            g_quant_coeff   [sizeId][listId][qp][dir] = (int32_t*)malloc(g_scalingListSize[sizeId]);
            //m_dequantCoef [sizeId][listId][qp][dir] = new Int [g_scalingListSize[sizeId]];
            //m_errScale    [sizeId][listId][qp][dir] = new double [g_scalingListSize[sizeId]];
          }
        }
      }
    }
  }
}

void scalinglist_destroy()
{
  uint32_t sizeId,listId,qp,dir;
  for(sizeId = 0; sizeId < 4; sizeId++)
  {
    for(listId = 0; listId < g_scalingListNum[sizeId]; listId++)
    {
      for(qp = 0; qp < 6; qp++)
      {
         free(g_quant_coeff[sizeId][listId][qp][0]);
        
        if(sizeId == /*SCALING_LIST_8x8*/1 || (sizeId == /*SCALING_LIST_16x16*/2 && listId < 2))
        {
          for(dir = /*SCALING_LIST_VER*/1; dir < /*SCALING_LIST_DIR_NUM*/3; dir++)
          {
            free(g_quant_coeff[sizeId][listId][qp][dir]);
          }
        }
      }
    }
  }
}

void processScalingListEnc( int32_t *coeff, int32_t *quantcoeff, int32_t quantScales, uint32_t height, uint32_t width, uint32_t ratio, int32_t sizuNum, uint32_t dc)
{
  uint32_t j,i;
  int32_t nsqth = (height < width) ? 4: 1; //height ratio for NSQT
  int32_t nsqtw = (width < height) ? 4: 1; //width ratio for NSQT
  for(j=0;j<height;j++)
  {
    for(i=0;i<width;i++)
    {
      quantcoeff[j*width + i] = quantScales / coeff[sizuNum * (j * nsqth / ratio) + i * nsqtw /ratio];
    }
  }
  if(ratio > 1)
  {
    quantcoeff[0] = quantScales / dc;
  }
}

void scalinglist_set(int32_t *coeff, uint32_t listId, uint32_t sizeId, uint32_t qp)
{
  uint32_t width = g_scalingListSizeX[sizeId];
  uint32_t height = g_scalingListSizeX[sizeId];
  uint32_t ratio = g_scalingListSizeX[sizeId]/MIN(8,g_scalingListSizeX[sizeId]);
  int32_t *quantcoeff;
  quantcoeff   = g_quant_coeff[listId][qp][sizeId][/*SCALING_LIST_SQT*/0];

  processScalingListEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);

  if(sizeId == /*SCALING_LIST_32x32*/3 || sizeId == /*SCALING_LIST_16x16*/2) //for NSQT
  {
    quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_VER*/1];
    processScalingListEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width>>2,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);

    quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_HOR*/2];
    processScalingListEnc(coeff,quantcoeff,g_quantScales[qp]<<4,height>>2,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);
  }
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


/** forward transform (2D)
*  \param block input residual
*  \param coeff transform coefficients
*  \param blockSize width of transform
*/
void transform2d(int16_t *block,int16_t *coeff, int8_t blockSize, int8_t uiMode)
{

  int32_t shift_1st = g_aucConvertToBit[blockSize]  + 1;// + g_uiBitIncrement; // log2(iWidth) - 1 + g_uiBitIncrement
  int32_t shift_2nd = g_aucConvertToBit[blockSize]  + 8;                   // log2(iHeight) + 6

  short tmp[ 64 * 64 ];
  /*
  if(blockSize== 4)
  {
    if (uiMode != REG_DCT)
    {
      fastForwardDst(block,tmp,shift_1st); // Forward DST BY FAST ALGORITHM, block input, tmp output
      fastForwardDst(tmp,coeff,shift_2nd); // Forward DST BY FAST ALGORITHM, tmp input, coeff output
    }
    else
    {
      partialButterfly4(block, tmp, shift_1st, iHeight);
      partialButterfly4(tmp, coeff, shift_2nd, iWidth);
    }

  }
  else*/ 
  if(blockSize == 8)
  {
    partialButterfly8( block, tmp, shift_1st, blockSize );
    partialButterfly8( tmp, coeff, shift_2nd, blockSize );
  }
  else if(blockSize == 16)
  {
    partialButterfly16( block, tmp, shift_1st, blockSize );
    partialButterfly16( tmp, coeff, shift_2nd, blockSize );
  }
  else if(blockSize == 32)
  {
    partialButterfly32( block, tmp, shift_1st, blockSize );
    partialButterfly32( tmp, coeff, shift_2nd, blockSize );
  }
}



#define QUANT_SHIFT 14

void quant(encoder_control* encoder, int16_t* pSrc, int16_t* pDes, /*int32_t** pArlDes,*/ int32_t iWidth,
           int32_t iHeight, uint32_t *uiAcSum, int8_t eTType/*, uint32_t uiAbsPartIdx*/ )
{
  int16_t*   piCoef    = pSrc;
  int16_t*   piQCoef   = pDes;
  //int32_t*   piArlCCoef = (*pArlDes);  
  uint32_t*  scan;
 
  int8_t useRDOQForTransformSkip = 0;
  //!(m_useTansformSkipFast && pcCU->getTransformSkip(uiAbsPartIdx,eTType));

  uint32_t log2BlockSize = g_aucConvertToBit[ iWidth ] + 2;

  uint32_t scanIdx = SCAN_DIAG;
  //pcCU->getCoefScanIdx(uiAbsPartIdx, iWidth, eTType==TEXT_LUMA, pcCU->isint32_tra(uiAbsPartIdx));
  if(scanIdx == SCAN_ZIGZAG)
  {
    scanIdx = SCAN_DIAG;
  }

  scan = g_auiSigLastScan[ scanIdx ][ log2BlockSize - 1 ];
  
  {
  int32_t deltaU[32*32] ;

  //QpParam cQpBase;
  int32_t iQpBase = encoder->QP;

  int32_t qpScaled;
  int32_t qpBDOffset = 0;//(eTType == 0)? pcCU->getSlice()->getSPS()->getQpBDOffsetY() : pcCU->getSlice()->getSPS()->getQpBDOffsetC();

  if(eTType == 1)
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
     //qpScaled = g_aucChromaScale[ qpScaled ] + qpBDOffset;
    }
  }
  
  //New block for variable definitions
  {
  int32_t n;
  uint32_t dir = 0;//SCALING_LIST_SQT;
    
  uint32_t uiLog2TrSize = g_aucConvertToBit[ iWidth ] + 2;
  int32_t scalingListType = (/*pcCU->isint32_tra(uiAbsPartIdx)*/0 ? 0 : 3) + ("\0\3\1\2"[eTType]);
  
  int32_t *piQuantCoeff = g_quant_coeff[scalingListType][/*m_cQP.m_iRem*/0][uiLog2TrSize-2][dir];

  uint32_t uiBitDepth = g_uiBitDepth + g_uiBitIncrement;

  int32_t iTransformShift = /*MAX_TR_DYNAMIC_RANGE*/15 - uiBitDepth - uiLog2TrSize; // Represents scaling through forward transform
  int32_t iQBits = QUANT_SHIFT + /*cQpBase.m_iPer +*/ iTransformShift;
  int32_t    iAdd = (encoder->in.cur_pic.type == NAL_IDR_SLICE ? 171 : 85) << (iQBits-9);

  int32_t iQBitsC = QUANT_SHIFT + /*cQpBase.m_iPer +*/ iTransformShift - /*ARL_C_PRECISION*/ 7;  
  int32_t iAddC   = 1 << (iQBitsC-1);

  int32_t qBits8 = iQBits-8;
  for( n = 0; n < iWidth*iHeight; n++ )
  {
    int32_t iLevel;
    int32_t  iSign;
    int64_t tmpLevel;
    uint32_t uiBlockPos = n;
    iLevel  = piCoef[uiBlockPos];
    iSign   = (iLevel < 0 ? -1: 1);      


    tmpLevel = (int64_t)abs(iLevel) * piQuantCoeff[uiBlockPos];
    /*
    if( m_bUseAdaptQpSelect )
    {
      piArlCCoef[uiBlockPos] = (int32_t)((tmpLevel + iAddC ) >> iQBitsC);
    }
    */
    iLevel = (int32_t)((tmpLevel + iAdd ) >> iQBits);
    deltaU[uiBlockPos] = (int32_t)((tmpLevel - (iLevel<<iQBits) )>> qBits8);

    uiAcSum += iLevel;
    iLevel *= iSign;        
    piQCoef[uiBlockPos] = MAX( -32768, MIN(32767, iLevel));
  } // for n
  }
  }

}