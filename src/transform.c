/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 * 
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see 
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 */

#include "transform.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "nal.h"

//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
// 
const int16_t g_t4[4][4] =
{
  { 64, 64, 64, 64},
  { 83, 36,-36,-83},
  { 64,-64,-64, 64},
  { 36,-83, 83,-36}
};

const int16_t g_t8[8][8] =
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

const int16_t g_t16[16][16] =
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

const int16_t g_t32[32][32] =
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

const int32_t g_quant_default_4x4[16] =
{
  16,16,16,16,
  16,16,16,16,
  16,16,16,16,
  16,16,16,16
};

const int32_t g_quant_intra_default_8x8[64] =
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

const int32_t g_quant_inter_default_8x8[64] =
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

const uint8_t g_chroma_scale[58]=
{
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
  17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,
  33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,
  45,46,47,48,49,50,51
};

int32_t *g_quant_coeff[4][6][6];
int32_t *g_de_quant_coeff[4][6][6];
double *g_error_scale[4][6][6];

const uint8_t g_scaling_list_num[4]    = { 6, 6, 6, 2};
const uint16_t g_scaling_list_size[4]  = {   16,  64, 256,1024}; 
const uint8_t g_scaling_list_size_x[4] = { 4, 8,16,32};
const int16_t g_quant_scales[6]        = { 26214,23302,20560,18396,16384,14564 };
const int16_t g_inv_quant_scales[6]    = { 40,45,51,57,64,72 };

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS
// 


/**
 * \brief Initialize scaling lists
 * 
 */
void scalinglist_init()
{
  uint32_t sizeId,listId,qp;

  for (sizeId = 0; sizeId < 4; sizeId++) {
    for (listId = 0; listId < g_scaling_list_num[sizeId]; listId++) {
      for (qp = 0; qp < 6; qp++) {
        if (!(sizeId == 3 && listId == 3)) {
          g_quant_coeff[sizeId][listId][qp]    = (int32_t*)calloc(g_scaling_list_size[sizeId], sizeof(int32_t));
          g_de_quant_coeff[sizeId][listId][qp] = (int32_t*)calloc(g_scaling_list_size[sizeId], sizeof(int32_t));
          g_error_scale[sizeId][listId][qp]    = (double*)calloc(g_scaling_list_size[sizeId], sizeof(double));
        }
      }
    }
  }
  // alias, assign pointer to an existing array
  for (qp = 0; qp < 6; qp++) {
    g_quant_coeff[3][3][qp]    = g_quant_coeff[3][1][qp];
    g_de_quant_coeff[3][3][qp] = g_de_quant_coeff[3][1][qp];
    g_error_scale[3][3][qp]    = g_error_scale[3][1][qp];
  }
}

/**
 * \brief Destroy scaling list allocated memory
 * 
 */
void scalinglist_destroy()
{
  uint32_t sizeId,listId,qp;

  for (sizeId = 0; sizeId < 4; sizeId++) {
    for (listId = 0; listId < g_scaling_list_num[sizeId]; listId++) {
      for (qp = 0; qp < 6; qp++) {
        if (!(sizeId == 3 && listId == 3)) {
          free(   g_quant_coeff[sizeId][listId][qp]);
          free(g_de_quant_coeff[sizeId][listId][qp]);
          free(   g_error_scale[sizeId][listId][qp]);
        }
      }
    }
  }
}


/**
 * \brief
 * 
 */
void scalinglist_process()
{
  #define SCALING_LIST_SIZE_NUM 4
  #define SCALING_LIST_REM_NUM 6
  uint32_t size,list,qp;

  for (size = 0; size < SCALING_LIST_SIZE_NUM; size++) {
    int32_t* list_ptr = (int32_t *)g_quant_intra_default_8x8; // Default to "8x8" intra

    for (list = 0; list < g_scaling_list_num[size]; list++) {
      switch(size) {
        case 0: // 4x4
          list_ptr = (int32_t *)g_quant_default_4x4;
          break;
        case 1: // 8x8
        case 2: // 16x16
          if (list > 2) list_ptr = (int32_t*)g_quant_inter_default_8x8;
          break;
        case 3: // 32x32
          if (list > 0) list_ptr = (int32_t*)g_quant_inter_default_8x8;
          break;
      }

      for (qp = 0; qp < SCALING_LIST_REM_NUM; qp++) {
        scalinglist_set(list_ptr, list, size, qp);
        scalinglist_set_err_scale(list, size, qp);
      }
    }
  }
  #undef SCALING_LIST_SIZE_NUM
  #undef SCALING_LIST_REM_NUM
}


/** set error scale coefficients
 * \param list List ID
 * \param uiSize Size
 * \param uiQP Quantization parameter
 */
#define MAX_TR_DYNAMIC_RANGE 15
void scalinglist_set_err_scale(uint32_t list,uint32_t size, uint32_t qp)
{
  uint32_t log2_tr_size   = g_convert_to_bit[ g_scaling_list_size_x[size] ] + 2;
  int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - g_bitdepth - log2_tr_size;  // Represents scaling through forward transform

  uint32_t i,max_num_coeff = g_scaling_list_size[size];
  int32_t *quantcoeff      = g_quant_coeff[size][list][qp];
  double *err_scale        = g_error_scale[size][list][qp];

  // Compensate for scaling of bitcount in Lagrange cost function
  double scale = (double)(1<<15);
  // Compensate for scaling through forward transform
  scale = scale*pow(2.0,-2.0*transform_shift);
  for(i=0;i<max_num_coeff;i++) {
    err_scale[i] = scale / quantcoeff[i] / quantcoeff[i] / (1<<(2*(g_bitdepth-8)));
  }
}

/**
 * \brief get scaling list for encoder
 * 
 */
void scalinglist_process_enc( int32_t *coeff, int32_t *quantcoeff, int32_t quant_scales, uint32_t height,uint32_t width, uint32_t ratio, int32_t size_num, uint32_t dc, uint8_t flat)
{
  uint32_t j,i;
  int32_t nsqth = (height < width) ? 4: 1; //!< height ratio for NSQT
  int32_t nsqtw = (width < height) ? 4: 1; //!< width ratio for NSQT

  // Flat scaling list
  if (flat) {
    for (j = 0; j < height * width; j++) {
      *quantcoeff++ = quant_scales>>4;
    }
  } else {
    for (j = 0; j < height; j++) {
      for (i = 0; i < width; i++) {
        uint32_t coeffpos  = size_num * (j * nsqth / ratio) + i * nsqtw / ratio;
        quantcoeff[j*width + i] = quant_scales / ((coeffpos > 63) ? 1 : coeff[coeffpos]);
      }
    }
    if (ratio > 1) {
      quantcoeff[0] = quant_scales / dc;
    }
  }
}

/**
 * \brief get scaling list for decoder
 * 
 */
void scalinglist_process_dec( int32_t *coeff, int32_t *dequantcoeff, int32_t inv_quant_scales, uint32_t height,uint32_t width, uint32_t ratio, int32_t size_num, uint32_t dc, uint8_t flat)
{
  uint32_t j,i;

  // Flat scaling list
  if (flat) {
    for (j = 0; j < height * width; j++) {
      *dequantcoeff++ = inv_quant_scales<<4;
    }
  } else {
    for (j = 0; j < height; j++) {
      for (i = 0; i < width; i++) {
        dequantcoeff[j*width + i] = inv_quant_scales * coeff[size_num * (j / ratio) + i / ratio];
      }
    }
    if (ratio > 1) {
      dequantcoeff[0] = inv_quant_scales * dc;
    }
  }
}

/**
 * \brief set scaling lists
 * 
 */
void scalinglist_set(int32_t *coeff, uint32_t listId, uint32_t sizeId, uint32_t qp)
{
  #define SCALING_LIST_DC 16
  uint32_t width  = g_scaling_list_size_x[sizeId];
  uint32_t height = g_scaling_list_size_x[sizeId];
  uint32_t ratio  = g_scaling_list_size_x[sizeId] / MIN(8, g_scaling_list_size_x[sizeId]);
  int32_t *quantcoeff   = g_quant_coeff[sizeId][listId][qp];
  int32_t *dequantcoeff = g_de_quant_coeff[sizeId][listId][qp];

  // Encoder list
  scalinglist_process_enc(coeff, quantcoeff, g_quant_scales[qp]<<4, height, width, ratio,
                          MIN(8, g_scaling_list_size_x[sizeId]), SCALING_LIST_DC, ENABLE_SCALING_LIST ? 0 : 1);
  // Decoder list
  scalinglist_process_dec(coeff, dequantcoeff, g_inv_quant_scales[qp], height, width, ratio,
                          MIN(8, g_scaling_list_size_x[sizeId]), SCALING_LIST_DC, ENABLE_SCALING_LIST ? 0 : 1);


  // TODO: support NSQT
  // if(sizeId == /*SCALING_LIST_32x32*/3 || sizeId == /*SCALING_LIST_16x16*/2) { //for NSQT
  //   quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_VER*/1];
  //   scalinglist_process_enc(coeff,quantcoeff,g_quantScales[qp]<<4,height,width>>2,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);

  //   quantcoeff   = g_quant_coeff[listId][qp][sizeId-1][/*SCALING_LIST_HOR*/2];
  //   scalinglist_process_enc(coeff,quantcoeff,g_quantScales[qp]<<4,height>>2,width,ratio,MIN(8,g_scalingListSizeX[sizeId]),/*scalingList->getScalingListDC(sizeId,listId)*/0);
  // }
  #undef SCALING_LIST_DC
}


void partial_butterfly_4(short *src, short *dst,int32_t shift, int32_t line)
{
  int32_t j;  
  int32_t e[2],o[2];
  int32_t add = 1<<(shift - 1);

  for (j = 0; j < line; j++) {    
    // E and O
    e[0] = src[0] + src[3];
    o[0] = src[0] - src[3];
    e[1] = src[1] + src[2];
    o[1] = src[1] - src[2];

    dst[0]      = (short)((g_t4[0][0]*e[0] + g_t4[0][1]*e[1] + add) >> shift);
    dst[2*line] = (short)((g_t4[2][0]*e[0] + g_t4[2][1]*e[1] + add) >> shift);
    dst[line]   = (short)((g_t4[1][0]*o[0] + g_t4[1][1]*o[1] + add) >> shift);
    dst[3*line] = (short)((g_t4[3][0]*o[0] + g_t4[3][1]*o[1] + add) >> shift);

    src += 4;
    dst ++;
  }
}

void partial_butterfly_inverse_4(short *src,short *dst,int shift, int line)
{
  int j;
  int e[2],o[2];
  int add = 1<<(shift - 1);

  for (j = 0; j < line; j++) {    
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    o[0] = g_t4[1][0]*src[line] + g_t4[3][0]*src[3*line];
    o[1] = g_t4[1][1]*src[line] + g_t4[3][1]*src[3*line];
    e[0] = g_t4[0][0]*src[0]    + g_t4[2][0]*src[2*line];
    e[1] = g_t4[0][1]*src[0]    + g_t4[2][1]*src[2*line];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    dst[0] = (short)CLIP(-32768, 32767, (e[0] + o[0] + add) >> shift);
    dst[1] = (short)CLIP(-32768, 32767, (e[1] + o[1] + add) >> shift);
    dst[2] = (short)CLIP(-32768, 32767, (e[1] - o[1] + add) >> shift);
    dst[3] = (short)CLIP(-32768, 32767, (e[0] - o[0] + add) >> shift);
            
    src++;
    dst += 4;
  }
}

// Fast DST Algorithm. Full matrix multiplication for DST and Fast DST algorithm 
// gives identical results
void fast_forward_dst(short *block, short *coeff, int32_t shift)  // input block, output coeff
{
  int32_t i, c[4];
  int32_t rnd_factor = 1<<(shift - 1);
  for (i = 0; i < 4; i++) {
    // int32_termediate Variables
    c[0] = block[4*i + 0] + block[4*i + 3];
    c[1] = block[4*i + 1] + block[4*i + 3];
    c[2] = block[4*i + 0] - block[4*i + 1];
    c[3] = 74* block[4*i + 2];

    coeff[   i] =  (short)(( 29*c[0] + 55*c[1]         + c[3]                     + rnd_factor ) >> shift);
    coeff[ 4+i] =  (short)(( 74*(block[4*i + 0]+ block[4*i + 1] - block[4*i + 3]) + rnd_factor ) >> shift);
    coeff[ 8+i] =  (short)(( 29*c[2] + 55*c[0]         - c[3]                     + rnd_factor ) >> shift);
    coeff[12+i] =  (short)(( 55*c[2] - 29*c[1]         + c[3]                     + rnd_factor ) >> shift);
  }
}

void fast_inverse_dst(short *tmp,short *block,int shift)  // input tmp, output block
{
  int i, c[4];
  int rnd_factor = 1<<(shift-1);
  for (i = 0; i < 4; i++) {  
    // Intermediate Variables
    c[0] = tmp[    i] + tmp[ 8 + i];
    c[1] = tmp[8 + i] + tmp[12 + i];
    c[2] = tmp[    i] - tmp[12 + i];
    c[3] = 74 * tmp[4 + i];

    block[4*i + 0] = (short)CLIP(-32768, 32767, ( 29*c[0] + 55*c[1]     + c[3]            + rnd_factor ) >> shift);
    block[4*i + 1] = (short)CLIP(-32768, 32767, ( 55*c[2] - 29*c[1]     + c[3]            + rnd_factor ) >> shift);
    block[4*i + 2] = (short)CLIP(-32768, 32767, ( 74*(tmp[i] - tmp[8 + i]  + tmp[12 + i]) + rnd_factor ) >> shift);
    block[4*i + 3] = (short)CLIP(-32768, 32767, ( 55*c[0] + 29*c[2]     - c[3]            + rnd_factor ) >> shift);
  }
}


void partial_butterfly_8(short *src, short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;  
  int32_t e[4],o[4];
  int32_t ee[2],eo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {  
    // E and O
    for (k = 0; k < 4; k++) {
      e[k] = src[k] + src[7 - k];
      o[k] = src[k] - src[7 - k];
    }    
    // EE and EO
    ee[0] = e[0] + e[3];
    eo[0] = e[0] - e[3];
    ee[1] = e[1] + e[2];
    eo[1] = e[1] - e[2];

    dst[0]      = (short)((g_t8[0][0]*ee[0] + g_t8[0][1]*ee[1] + add) >> shift);
    dst[4*line] = (short)((g_t8[4][0]*ee[0] + g_t8[4][1]*ee[1] + add) >> shift); 
    dst[2*line] = (short)((g_t8[2][0]*eo[0] + g_t8[2][1]*eo[1] + add) >> shift);
    dst[6*line] = (short)((g_t8[6][0]*eo[0] + g_t8[6][1]*eo[1] + add) >> shift); 

    dst[line]   = (short)((g_t8[1][0]*o[0] + g_t8[1][1]*o[1] + g_t8[1][2]*o[2] + g_t8[1][3]*o[3] + add) >> shift);
    dst[3*line] = (short)((g_t8[3][0]*o[0] + g_t8[3][1]*o[1] + g_t8[3][2]*o[2] + g_t8[3][3]*o[3] + add) >> shift);
    dst[5*line] = (short)((g_t8[5][0]*o[0] + g_t8[5][1]*o[1] + g_t8[5][2]*o[2] + g_t8[5][3]*o[3] + add) >> shift);
    dst[7*line] = (short)((g_t8[7][0]*o[0] + g_t8[7][1]*o[1] + g_t8[7][2]*o[2] + g_t8[7][3]*o[3] + add) >> shift);

    src += 8;
    dst++;
  }
}

void partial_butterfly_inverse_8(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[4],o[4];
  int32_t ee[2],eo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {

    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 4; k++) {
      o[k] = g_t8[ 1][k]*src[line] + g_t8[ 3][k]*src[3*line] + g_t8[ 5][k]*src[5*line] + g_t8[ 7][k]*src[7*line];
    }

    eo[0] = g_t8[2][0]*src[ 2*line ] + g_t8[6][0]*src[ 6*line ];
    eo[1] = g_t8[2][1]*src[ 2*line ] + g_t8[6][1]*src[ 6*line ];
    ee[0] = g_t8[0][0]*src[ 0      ] + g_t8[4][0]*src[ 4*line ];
    ee[1] = g_t8[0][1]*src[ 0      ] + g_t8[4][1]*src[ 4*line ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    e[0] = ee[0] + eo[0];
    e[3] = ee[0] - eo[0];
    e[1] = ee[1] + eo[1];
    e[2] = ee[1] - eo[1];
    for (k = 0; k < 4; k++) {
      dst[ k   ] = (int16_t)MAX(-32768, MIN(32767, (e[k] + o[k] + add)>>shift));
      dst[ k+4 ] = (int16_t)MAX(-32768, MIN(32767, (e[3-k] - o[3-k] + add)>>shift));
    }
    src++;
    dst += 8;
  }
}


void partial_butterfly_16(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[8],o[8];
  int32_t ee[4],eo[4];
  int32_t eee[2],eeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {    
    // E and O
    for (k = 0; k < 8; k++) {
      e[k] = src[k] + src[15 - k];
      o[k] = src[k] - src[15 - k];
    } 
    // EE and EO
    for (k = 0; k < 4; k++) {
      ee[k] = e[k] + e[7 - k];
      eo[k] = e[k] - e[7 - k];
    }
    // EEE and EEO
    eee[0] = ee[0] + ee[3];
    eeo[0] = ee[0] - ee[3];
    eee[1] = ee[1] + ee[2];
    eeo[1] = ee[1] - ee[2];

    dst[0      ] = (short)((g_t16[ 0][0]*eee[0] + g_t16[ 0][1]*eee[1] + add) >> shift);
    dst[8*line ] = (short)((g_t16[ 8][0]*eee[0] + g_t16[ 8][1]*eee[1] + add) >> shift);
    dst[4*line ] = (short)((g_t16[ 4][0]*eeo[0] + g_t16[ 4][1]*eeo[1] + add) >> shift);
    dst[12*line] = (short)((g_t16[12][0]*eeo[0] + g_t16[12][1]*eeo[1] + add) >> shift);

    for (k = 2; k < 16; k += 4) {
      dst[k*line] = (short)((g_t16[k][0]*eo[0] + g_t16[k][1]*eo[1] + g_t16[k][2]*eo[2] + g_t16[k][3]*eo[3] + add) >> shift);
    }

    for (k = 1; k < 16; k += 2) {
      dst[k*line] = (short)((g_t16[k][0]*o[0] + g_t16[k][1]*o[1] + g_t16[k][2]*o[2] + g_t16[k][3]*o[3] + 
                     g_t16[k][4]*o[4] + g_t16[k][5]*o[5] + g_t16[k][6]*o[6] + g_t16[k][7]*o[7] + add) >> shift);
    }

    src += 16;
    dst++;
  }
}


void partial_butterfly_inverse_16(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[8],o[8];
  int32_t ee[4],eo[4];
  int32_t eee[2],eeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {    
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 8; k++)  {
      o[k] = g_t16[ 1][k]*src[  line] + g_t16[ 3][k]*src[ 3*line] + g_t16[ 5][k]*src[ 5*line] + g_t16[ 7][k]*src[ 7*line] + 
             g_t16[ 9][k]*src[9*line] + g_t16[11][k]*src[11*line] + g_t16[13][k]*src[13*line] + g_t16[15][k]*src[15*line];
    }
    for (k = 0; k < 4; k++) {
      eo[k] = g_t16[ 2][k]*src[ 2*line] + g_t16[ 6][k]*src[ 6*line] + g_t16[10][k]*src[10*line] + g_t16[14][k]*src[14*line];
    }
    eeo[0] = g_t16[4][0]*src[ 4*line ] + g_t16[12][0]*src[ 12*line ];
    eee[0] = g_t16[0][0]*src[ 0      ] + g_t16[ 8][0]*src[ 8*line  ];
    eeo[1] = g_t16[4][1]*src[ 4*line ] + g_t16[12][1]*src[ 12*line ];
    eee[1] = g_t16[0][1]*src[ 0      ] + g_t16[ 8][1]*src[ 8*line  ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    for (k = 0; k < 2; k++) {
      ee[k]   = eee[k]   + eeo[k];
      ee[k+2] = eee[1-k] - eeo[1-k];
    }    
    for (k = 0; k < 4; k++) {
      e[k]   = ee[k] + eo[k];
      e[k+4] = ee[3-k] - eo[3-k];
    }    
    for (k = 0; k < 8; k++) {
      dst[k]   = (short)MAX(-32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k+8] = (short)MAX(-32768, MIN(32767, (e[7-k] - o[7-k] + add) >> shift));
    }
    src++;
    dst += 16;
  }
}



void partial_butterfly_32(short *src,short *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[16],o[16];
  int32_t ee[8],eo[8];
  int32_t eee[4],eeo[4];
  int32_t eeee[2],eeeo[2];
  int32_t add = 1<<(shift-1);

  for (j = 0; j < line; j++) {    
    // E and O
    for (k = 0; k < 16; k++) {
      e[k] = src[k] + src[31-k];
      o[k] = src[k] - src[31-k];
    }
    // EE and EO
    for (k = 0; k < 8; k++) {
      ee[k] = e[k] + e[15-k];
      eo[k] = e[k] - e[15-k];
    }
    // EEE and EEO
    for (k = 0; k < 4; k++) {
      eee[k] = ee[k] + ee[7-k];
      eeo[k] = ee[k] - ee[7-k];
    }
    // EEEE and EEEO
    eeee[0] = eee[0] + eee[3];
    eeeo[0] = eee[0] - eee[3];
    eeee[1] = eee[1] + eee[2];
    eeeo[1] = eee[1] - eee[2];

    dst[0      ] = (short)((g_t32[ 0][0]*eeee[0] + g_t32[ 0][1]*eeee[1] + add) >> shift);
    dst[16*line] = (short)((g_t32[16][0]*eeee[0] + g_t32[16][1]*eeee[1] + add) >> shift);
    dst[ 8*line] = (short)((g_t32[ 8][0]*eeeo[0] + g_t32[ 8][1]*eeeo[1] + add) >> shift);
    dst[24*line] = (short)((g_t32[24][0]*eeeo[0] + g_t32[24][1]*eeeo[1] + add) >> shift);
    for (k = 4; k < 32; k += 8) {
      dst[k*line] = (short)((g_t32[k][0]*eeo[0] + g_t32[k][1]*eeo[1] + g_t32[k][2]*eeo[2] + g_t32[k][3]*eeo[3] + add) >> shift);
    }
    for (k = 2; k < 32; k += 4) {
      dst[k*line] = (short)((g_t32[k][0]*eo[0] + g_t32[k][1]*eo[1] + g_t32[k][2]*eo[2] + g_t32[k][3]*eo[3] +
                             g_t32[k][4]*eo[4] + g_t32[k][5]*eo[5] + g_t32[k][6]*eo[6] + g_t32[k][7]*eo[7] + add) >> shift);
    }
    for (k = 1; k < 32; k += 2) {
      dst[k*line] = (short)((g_t32[k][ 0]*o[ 0] + g_t32[k][ 1]*o[ 1] + g_t32[k][ 2]*o[ 2] + g_t32[k][ 3]*o[ 3] +
                             g_t32[k][ 4]*o[ 4] + g_t32[k][ 5]*o[ 5] + g_t32[k][ 6]*o[ 6] + g_t32[k][ 7]*o[ 7] +
                             g_t32[k][ 8]*o[ 8] + g_t32[k][ 9]*o[ 9] + g_t32[k][10]*o[10] + g_t32[k][11]*o[11] + 
                             g_t32[k][12]*o[12] + g_t32[k][13]*o[13] + g_t32[k][14]*o[14] + g_t32[k][15]*o[15] + add) >> shift);
    }
    src += 32;
    dst++;
  }
}


void partial_butterfly_inverse_32(int16_t *src,int16_t *dst,int32_t shift, int32_t line)
{
  int32_t j,k;
  int32_t e[16],o[16];
  int32_t ee[8],eo[8];
  int32_t eee[4],eeo[4];
  int32_t eeee[2],eeeo[2];
  int32_t add = 1<<(shift-1);

  for (j=0; j<line; j++) {
    // Utilizing symmetry properties to the maximum to minimize the number of multiplications
    for (k = 0; k < 16; k++) {
      o[k] = g_t32[ 1][k]*src[ line  ]   + g_t32[ 3][k]*src[ 3*line  ] + g_t32[ 5][k]*src[ 5*line  ] + g_t32[ 7][k]*src[ 7*line  ] + 
             g_t32[ 9][k]*src[ 9*line  ] + g_t32[11][k]*src[ 11*line ] + g_t32[13][k]*src[ 13*line ] + g_t32[15][k]*src[ 15*line ] + 
             g_t32[17][k]*src[ 17*line ] + g_t32[19][k]*src[ 19*line ] + g_t32[21][k]*src[ 21*line ] + g_t32[23][k]*src[ 23*line ] + 
             g_t32[25][k]*src[ 25*line ] + g_t32[27][k]*src[ 27*line ] + g_t32[29][k]*src[ 29*line ] + g_t32[31][k]*src[ 31*line ];
    }
    for (k = 0; k < 8; k++) {
      eo[k] = g_t32[ 2][k]*src[ 2*line  ] + g_t32[ 6][k]*src[ 6*line  ] + g_t32[10][k]*src[ 10*line ] + g_t32[14][k]*src[ 14*line ] + 
              g_t32[18][k]*src[ 18*line ] + g_t32[22][k]*src[ 22*line ] + g_t32[26][k]*src[ 26*line ] + g_t32[30][k]*src[ 30*line ];
    }
    for (k = 0; k < 4; k++) {
      eeo[k] = g_t32[4][k]*src[ 4*line ] + g_t32[12][k]*src[ 12*line ] + g_t32[20][k]*src[ 20*line ] + g_t32[28][k]*src[ 28*line ];
    }
    eeeo[0] = g_t32[8][0]*src[ 8*line ] + g_t32[24][0]*src[ 24*line ];
    eeeo[1] = g_t32[8][1]*src[ 8*line ] + g_t32[24][1]*src[ 24*line ];
    eeee[0] = g_t32[0][0]*src[ 0      ] + g_t32[16][0]*src[ 16*line ];    
    eeee[1] = g_t32[0][1]*src[ 0      ] + g_t32[16][1]*src[ 16*line ];

    // Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector
    eee[0] = eeee[0] + eeeo[0];
    eee[3] = eeee[0] - eeeo[0];
    eee[1] = eeee[1] + eeeo[1];
    eee[2] = eeee[1] - eeeo[1];    
    for (k = 0; k < 4; k++) {
      ee[k]   = eee[k]   + eeo[k];
      ee[k+4] = eee[3-k] - eeo[3-k];
    }
    for (k = 0; k < 8; k++) {
      e[k]   = ee[k]   + eo[k];
      e[k+8] = ee[7-k] - eo[7-k];
    }
    for (k=0;k<16;k++) {
      dst[k]    = (short)MAX( -32768, MIN(32767, (e[k] + o[k] + add) >> shift));
      dst[k+16] = (short)MAX( -32768, MIN(32767, (e[15-k] - o[15-k] + add) >> shift));
    }
    src++;
    dst += 32;
  }
}


/** 
 * \brief forward transform (2D)
 * \param block input residual
 * \param coeff transform coefficients
 * \param block_size width of transform
 */
void transform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode)
{

  int32_t shift_1st = g_convert_to_bit[block_size]  + 1 + g_bit_increment;
  int32_t shift_2nd = g_convert_to_bit[block_size]  + 8;

  int16_t tmp[LCU_WIDTH * LCU_WIDTH];
  
  if(block_size== 4) {
    if (mode != 65535) {
      // Forward DST BY FAST ALGORITHM
      fast_forward_dst(block,tmp,shift_1st);
      fast_forward_dst(tmp,coeff,shift_2nd);
    } else {
      partial_butterfly_4(block, tmp, shift_1st, block_size);
      partial_butterfly_4(tmp, coeff, shift_2nd, block_size);
    }
  } else {
    switch(block_size) {
      case 8:
        partial_butterfly_8( block, tmp, shift_1st, block_size );
        partial_butterfly_8( tmp, coeff, shift_2nd, block_size );
        break;
      case 16:
        partial_butterfly_16( block, tmp, shift_1st, block_size );
        partial_butterfly_16( tmp, coeff, shift_2nd, block_size );
        break;      
      case 32:
        partial_butterfly_32( block, tmp, shift_1st, block_size );
        partial_butterfly_32( tmp, coeff, shift_2nd, block_size );
        break;
    }
  }
}

/**
 * \brief NxN inverse transform (2D)
 * \param coeff input data (transform coefficients)
 * \param block output data (residual)
 * \param block_size input data (width of transform)
 * \param mode
 */
void itransform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode)  
{
  int32_t shift_1st = 7;
  int32_t shift_2nd = 12 - (g_bitdepth - 8);
  int16_t tmp[LCU_WIDTH*LCU_WIDTH];

  if( block_size == 4) {
    if (mode != 65535) {
      // Inverse DST by FAST Algorithm
      fast_inverse_dst(coeff, tmp, shift_1st);
      fast_inverse_dst(tmp, block, shift_2nd);
    } else {
      partial_butterfly_inverse_4(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_4(tmp, block, shift_2nd, block_size);
    }
  } else {
    switch(block_size) {
    case 8:
      partial_butterfly_inverse_8(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_8(tmp, block, shift_2nd, block_size);
      break;
    case 16:
      partial_butterfly_inverse_16(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_16(tmp, block, shift_2nd, block_size);
      break;      
    case 32:
      partial_butterfly_inverse_32(coeff, tmp, shift_1st, block_size);
      partial_butterfly_inverse_32(tmp, block, shift_2nd, block_size);
      break;
    }
  }
}


#define QUANT_SHIFT 14
#define MAX_TR_DYNAMIC_RANGE 15
/**
 * \brief quantize transformed coefficents
 * 
 */
void quant(encoder_control *encoder, int16_t *coef, int16_t *q_coef, int32_t width,
           int32_t height, uint32_t *ac_sum, int8_t type, int8_t scan_idx, int8_t block_type )
{
  uint32_t log2_block_size = g_convert_to_bit[ width ] + 2;
  uint32_t *scan = g_sig_last_scan[ scan_idx ][ log2_block_size - 1 ];

  #if ENABLE_SIGN_HIDING == 1
  int32_t delta_u[LCU_WIDTH*LCU_WIDTH>>2];
  #endif
  int32_t qp_base = encoder->QP;

  int32_t qp_scaled;
  int32_t qp_offset = 0;
  if(type == 0) {
    qp_scaled = qp_base + qp_offset;
  } else {
    qp_scaled = CLIP(-qp_offset, 57, qp_base);
    if(qp_scaled < 0) {
      qp_scaled = qp_scaled + qp_offset;
    } else {
      qp_scaled = g_chroma_scale[qp_scaled] + qp_offset;
    }
  }

  //New block for variable definitions
  {
  int32_t n; 
  uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
  int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);
  
  int32_t *quant_coeff = g_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];

  int32_t transform_shift = MAX_TR_DYNAMIC_RANGE - g_bitdepth - log2_tr_size; //!< Represents scaling through forward transform
  int32_t q_bits = QUANT_SHIFT + qp_scaled/6 + transform_shift;
  int32_t add = ((encoder->in.cur_pic->slicetype == SLICE_I) ? 171 : 85) << (q_bits - 9);

  int32_t q_bits8 = q_bits - 8;
  for (n = 0; n < width * height; n++) {
    int32_t level;
    int32_t  sign;

    level = coef[n];
    sign  = (level < 0 ? -1: 1);

    level = ((int64_t)abs(level) * quant_coeff[n] + add ) >> q_bits;

    #if ENABLE_SIGN_HIDING == 1
    delta_u[n] = (int32_t)( ((int64_t)abs(coef[n]) * quant_coeff[n] - (level<<q_bits) )>> q_bits8 );
    *ac_sum += level;
    #endif

    level *= sign;
    q_coef[n] = (int16_t)(CLIP( -32768, 32767, level));
  }

  #if ENABLE_SIGN_HIDING == 1
  if(*ac_sum >= 2) {
    #define SCAN_SET_SIZE 16
    #define LOG2_SCAN_SET_SIZE 4
    int32_t n,last_cg = -1, abssum = 0, subset, subpos;
    for(subset = (width*height - 1)>>LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
      int32_t first_nz_pos_in_cg = SCAN_SET_SIZE, last_nz_pos_in_cg=-1;
      subpos = subset<<LOG2_SCAN_SET_SIZE;
      abssum = 0;

      // Find last coeff pos
      for (n = SCAN_SET_SIZE - 1; n >= 0; n--)  {
        if (q_coef[scan[n + subpos]])  {
          last_nz_pos_in_cg = n;
          break;
        }
      }

      // First coeff pos
      for (n = 0; n <SCAN_SET_SIZE; n++) {
        if (q_coef[scan[n + subpos]]) {
          first_nz_pos_in_cg = n;
          break;
        }
      }

      // Sum all quant coeffs between first and last
      for(n = first_nz_pos_in_cg; n <= last_nz_pos_in_cg; n++) {
        abssum += q_coef[scan[n + subpos]];
      }

      if(last_nz_pos_in_cg >= 0 && last_cg == -1) {
        last_cg = 1;
      }

      if(last_nz_pos_in_cg - first_nz_pos_in_cg >= 4) {
        int32_t signbit = (q_coef[scan[subpos + first_nz_pos_in_cg]] > 0 ? 0 : 1) ;
        if(signbit != (abssum&0x1)) { // compare signbit with sum_parity
          int32_t min_cost_inc = 0x7fffffff,  min_pos =-1, cur_cost=0x7fffffff;
          int16_t final_change = 0, cur_change=0;
          for(n = (last_cg == 1 ? last_nz_pos_in_cg : SCAN_SET_SIZE - 1); n >= 0; n--) {
            uint32_t blkPos  = scan[n + subpos];
            if(q_coef[blkPos] != 0) {
              if(delta_u[blkPos] > 0) {
                cur_cost = -delta_u[blkPos];
                cur_change=1;
              } else if(n == first_nz_pos_in_cg && abs(q_coef[blkPos]) == 1) {
                cur_cost=0x7fffffff;
              } else {
                cur_cost = delta_u[blkPos];
                cur_change =-1;
              }
            } else if(n < first_nz_pos_in_cg && ((coef[blkPos] >= 0)?0:1) != signbit) {
              cur_cost = 0x7fffffff;
            } else {
              cur_cost   = -delta_u[blkPos];
              cur_change = 1;
            }

            if(cur_cost < min_cost_inc) {
              min_cost_inc = cur_cost;
              final_change = cur_change;
              min_pos      = blkPos;
            }
          } // CG loop

          if(q_coef[min_pos] == 32767 || q_coef[min_pos] == -32768) {
            final_change = -1;
          }

          if(coef[min_pos] >= 0) q_coef[min_pos] += final_change;
          else q_coef[min_pos] -= final_change;

        } // Hide
      }
      if (last_cg == 1) last_cg=0;      
    }

    #undef SCAN_SET_SIZE
    #undef LOG2_SCAN_SET_SIZE
  }
  #endif
  }
}

/**
 * \brief inverse quantize transformed and quantized coefficents
 * 
 */
void dequant(encoder_control *encoder, int16_t *q_coef, int16_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type)
{
  int32_t shift,add,coeff_q,clip_q_coef;
  int32_t n;
  int32_t transform_shift = 15 - g_bitdepth - (g_convert_to_bit[ width ] + 2);
  int32_t qp_scaled;
  int32_t qp_base = encoder->QP;

  if (type == 0) {
    qp_scaled = qp_base;
  } else {
    qp_scaled = CLIP( 0, 57, qp_base);
    if (qp_scaled < 0) {
      qp_scaled = qp_scaled;
    } else {
      qp_scaled = g_chroma_scale[qp_scaled];
    }
  }
  
  shift = 20 - QUANT_SHIFT - transform_shift;

  UNREFERENCED_PARAMETER(block_type);
  #if ENABLE_SCALING_LIST == 1
  {
    uint32_t log2_tr_size = g_convert_to_bit[ width ] + 2;
    int32_t scalinglist_type = (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

    dequant_coef = g_de_quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
    shift += 4;

    if (shift >qp_scaled / 6) {
      add = 1 << (shift - qp_scaled/6 - 1);
    
      for (n = 0; n < width * height; n++) {
        clip_q_coef = CLIP(-32768, 32767, q_coef[n]);
        coeff_q = ((clip_q_coef * dequant_coef[n]) + add ) >> (shift -  qp_scaled/6);
        coef[n] = (int16_t)CLIP(-32768,32767,coeff_q);
      }
    } else {
      for (n = 0; n < width * height; n++) {
        // Clip to avoid possible overflow in following shift left operation
        clip_q_coef = CLIP(-32768, 32767, q_coef[n]);
        coeff_q   = CLIP(-32768, 32767, clip_q_coef * dequant_coef[n]); 
        coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q << (qp_scaled/6 - shift));
      }
    }
  }
  #else
  {
    int32_t scale = g_inv_quant_scales[qp_scaled%6] << (qp_scaled/6);
    add = 1 << (shift-1);

    for (n = 0; n < width*height; n++) {
      clip_q_coef = CLIP(-32768, 32767, q_coef[n]);
      coeff_q   = (clip_q_coef * scale + add) >> shift;
      coef[n] = (int16_t)CLIP(-32768, 32767, coeff_q);
    }
  }
  #endif
}
