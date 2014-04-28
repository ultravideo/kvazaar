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

#include "picture.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "sao.h"


#define PSNRMAX (255.0 * 255.0)


yuv_t * alloc_yuv_t(int luma_size)
{
  // Get buffers with separate mallocs in order to take advantage of
  // automatic buffer overrun checks.
  yuv_t *yuv = (yuv_t *)malloc(sizeof(*yuv));
  yuv->y = (pixel *)malloc(luma_size * sizeof(*yuv->y));
  yuv->u = (pixel *)malloc(luma_size / 2 * sizeof(*yuv->u));
  yuv->v = (pixel *)malloc(luma_size / 2 * sizeof(*yuv->v));
  yuv->size = luma_size;

  return yuv;
}

void dealloc_yuv_t(yuv_t * yuv)
{
  free(yuv->y);
  free(yuv->u);
  free(yuv->v);
  free(yuv);
}


/**
 * \brief BLock Image Transfer from one buffer to another.
 *
 * It's a stupidly simple loop that copies pixels.
 *
 * \param orig  Start of the originating buffer.
 * \param dst  Start of the destination buffer.
 * \param width  Width of the copied region.
 * \param height  Height of the copied region.
 * \param orig_stride  Width of a row in the originating buffer.
 * \param dst_stride  Width of a row in the destination buffer.
 *
 * This should be inlined, but it's defined here for now to see if Visual
 * Studios LTCG will inline it.
 */
void picture_blit_pixels(const pixel *orig, pixel *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride)
{
  unsigned y, x;
  //Fix problem with reading unitialized memory
  if (width > orig_stride) width = orig_stride;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      dst[x] = orig[x];
    }
    // Move pointers to the next row.
    orig += orig_stride;
    dst += dst_stride;
  }
}

void picture_blit_coeffs(const coefficient *orig, coefficient *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride)
{
  unsigned y, x;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      dst[x] = orig[x];
    }
    // Move pointers to the next row.
    orig += orig_stride;
    dst += dst_stride;
  }
}


/**
 * \brief Allocate memory for picture_list
 * \param size  initial array size
 * \return picture_list pointer, NULL on failure
 */
picture_list * picture_list_init(int size)
{
  picture_list *list = (picture_list *)malloc(sizeof(picture_list));
  list->size = size;
  if (size > 0) {
    list->pics = (picture**)malloc(sizeof(picture*) * size);
  }

  list->used_size = 0;

  return list;
}

/**
 * \brief Resize picture_list array
 * \param list  picture_list pointer
 * \param size  new array size
 * \return 1 on success, 0 on failure
 */
int picture_list_resize(picture_list *list, unsigned size)
{
  unsigned int i;
  picture** old_pics = NULL;

  // No need to do anything when resizing to same size
  if (size == list->size) {
    return 1;
  }

  // Save the old list
  if (list->used_size > 0) {
    old_pics = list->pics;
  }

  // allocate space for the new list
  list->pics = (picture**)malloc(sizeof(picture*)*size);

  // Copy everything from the old list to the new if needed.
  if (old_pics != NULL) {
    for (i = 0; i < list->used_size; ++i) {
      list->pics[i] = old_pics[i];
    }

    free(old_pics);
  }

  return 1;
}

/**
 * \brief Free memory allocated to the picture_list
 * \param list picture_list pointer
 * \return 1 on success, 0 on failure
 */
int picture_list_destroy(picture_list *list)
{
  unsigned int i;
  if (list->used_size > 0) {
    for (i = 0; i < list->used_size; ++i) {
      picture_destroy(list->pics[i]);
      FREE_POINTER(list->pics[i]);
    }
  }

  if (list->size > 0) {
    free(list->pics);
  }
  free(list);
  return 1;
}

/**
 * \brief Add picture to the front of the picturelist
 * \param pic picture pointer to add
 * \param picture_list list to use
 * \return 1 on success
 */
int picture_list_add(picture_list *list,picture* pic)
{
  int i = 0;
  if (list->size == list->used_size) {
    if (!picture_list_resize(list, list->size*2)) return 0;
  }

  for (i = list->used_size; i > 0; i--) {
    list->pics[i] = list->pics[i - 1];
  }

  list->pics[0] = pic;
  list->used_size++;
  return 1;
}

/**
 * \brief Add picture to picturelist
 * \param pic picture pointer to add
 * \param picture_list list to use
 * \return 1 on success
 */
int picture_list_rem(picture_list *list, unsigned n, int8_t destroy)
{
  // Must be within list boundaries
  if (n >= list->used_size)
  {
    return 0;
  }

  if (destroy) {
    picture_destroy(list->pics[n]);
    FREE_POINTER(list->pics[n]);
  }

  // The last item is easy to remove
  if (n == list->used_size - 1) {
    list->pics[n] = NULL;
    list->used_size--;
  } else {
    // Shift all following pics one backward in the list
    for (; n < list->used_size - 1; ++n) {
      list->pics[n] = list->pics[n + 1];
    }
    list->pics[list->used_size - 1] = NULL;
    list->used_size--;
  }

  return 1;
}

/**
 * \brief Allocate new picture
 * \param pic picture pointer
 * \return picture pointer
 */
picture *picture_init(int32_t width, int32_t height,
                      int32_t width_in_lcu, int32_t height_in_lcu)
{
  picture *pic = (picture *)malloc(sizeof(picture));
  unsigned int luma_size = width * height;
  unsigned int chroma_size = luma_size / 4;

  if (!pic) return 0;

  memset(pic, 0, sizeof(picture));

  pic->width  = width;
  pic->height = height;
  pic->width_in_lcu  = width_in_lcu;
  pic->height_in_lcu = height_in_lcu;
  pic->referenced = 0;
  // Allocate buffers
  pic->y_data = MALLOC(pixel, luma_size);
  pic->u_data = MALLOC(pixel, chroma_size);
  pic->v_data = MALLOC(pixel, chroma_size);
  pic->data[COLOR_Y] = pic->y_data;
  pic->data[COLOR_U] = pic->u_data;
  pic->data[COLOR_V] = pic->v_data;

  // Reconstruction buffers
  pic->y_recdata = MALLOC(pixel, luma_size);
  pic->u_recdata = MALLOC(pixel, chroma_size);
  pic->v_recdata = MALLOC(pixel, chroma_size);
  pic->recdata[COLOR_Y] = pic->y_recdata;
  pic->recdata[COLOR_U] = pic->u_recdata;
  pic->recdata[COLOR_V] = pic->v_recdata;

  memset(pic->u_recdata, 128, (chroma_size));
  memset(pic->v_recdata, 128, (chroma_size));

  {
    // Allocate height_in_scu x width_in_scu x sizeof(CU_info)
    unsigned height_in_scu = height_in_lcu << MAX_DEPTH;
    unsigned width_in_scu = width_in_lcu << MAX_DEPTH;
    unsigned cu_array_size = height_in_scu * width_in_scu;
    pic->cu_array = (cu_info*)malloc(sizeof(cu_info) * cu_array_size);
    memset(pic->cu_array, 0, sizeof(cu_info) * cu_array_size);
  }
  
  pic->coeff_y = NULL; pic->coeff_u = NULL; pic->coeff_v = NULL;
  pic->pred_y = NULL; pic->pred_u = NULL; pic->pred_v = NULL;

  pic->slice_sao_luma_flag = 1;
  pic->slice_sao_chroma_flag = 1;
  pic->sao_luma = MALLOC(sao_info, width_in_lcu * height_in_lcu);
  pic->sao_chroma = MALLOC(sao_info, width_in_lcu * height_in_lcu);

  return pic;
}

/**
 * \brief Free memory allocated to picture
 * \param pic picture pointer
 * \return 1 on success, 0 on failure
 */
int picture_destroy(picture *pic)
{
  free(pic->u_data);
  free(pic->v_data);
  free(pic->y_data);
  pic->y_data = pic->u_data = pic->v_data = NULL;
  pic->data[COLOR_Y] = pic->data[COLOR_U] = pic->data[COLOR_V] = NULL;

  free(pic->y_recdata);
  free(pic->u_recdata);
  free(pic->v_recdata);
  pic->y_recdata = pic->u_recdata = pic->v_recdata = NULL;
  pic->recdata[COLOR_Y] = pic->recdata[COLOR_U] = pic->recdata[COLOR_V] = NULL;

  FREE_POINTER(pic->cu_array);

  FREE_POINTER(pic->coeff_y);
  FREE_POINTER(pic->coeff_u);
  FREE_POINTER(pic->coeff_v);

  FREE_POINTER(pic->pred_y);
  FREE_POINTER(pic->pred_u);
  FREE_POINTER(pic->pred_v);

  FREE_POINTER(pic->sao_luma);
  FREE_POINTER(pic->sao_chroma);

  return 1;
}

/**
 * \brief Calculates image PSNR value
 */
double image_psnr(pixel *frame1, pixel *frame2, int32_t x, int32_t y)
{
  uint64_t error_sum = 0;
  int32_t error = 0;
  int32_t pixels = x * y;
  int32_t i;

  for (i = 0; i < pixels; ++i) {
    error = frame1[i] - frame2[i];
    error_sum += error * error;
  }

  // Avoid division by zero
  if (error_sum == 0) return 99.0;

  return 10 * log10((pixels * PSNRMAX) / ((double)error_sum));
}

/**
 * \brief  Calculate SATD between two 4x4 blocks inside bigger arrays.
 * From HM 13.0
 */
static unsigned satd_16bit_4x4(pixel *piOrg,pixel *piCur)
{
  int32_t k, satd = 0, diff[16], m[16], d[16];
  int32_t iStrideOrg = 4, iStrideCur = 4;
  for( k = 0; k < 16; k+=4 ) {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  /*===== hadamard transform =====*/
  m[ 0] = diff[ 0] + diff[12];
  m[ 1] = diff[ 1] + diff[13];
  m[ 2] = diff[ 2] + diff[14];
  m[ 3] = diff[ 3] + diff[15];
  m[ 4] = diff[ 4] + diff[ 8];
  m[ 5] = diff[ 5] + diff[ 9];
  m[ 6] = diff[ 6] + diff[10];
  m[ 7] = diff[ 7] + diff[11];
  m[ 8] = diff[ 4] - diff[ 8];
  m[ 9] = diff[ 5] - diff[ 9];
  m[10] = diff[ 6] - diff[10];
  m[11] = diff[ 7] - diff[11];
  m[12] = diff[ 0] - diff[12];
  m[13] = diff[ 1] - diff[13];
  m[14] = diff[ 2] - diff[14];
  m[15] = diff[ 3] - diff[15];

  d[ 0] = m[ 0] + m[ 4];
  d[ 1] = m[ 1] + m[ 5];
  d[ 2] = m[ 2] + m[ 6];
  d[ 3] = m[ 3] + m[ 7];
  d[ 4] = m[ 8] + m[12];
  d[ 5] = m[ 9] + m[13];
  d[ 6] = m[10] + m[14];
  d[ 7] = m[11] + m[15];
  d[ 8] = m[ 0] - m[ 4];
  d[ 9] = m[ 1] - m[ 5];
  d[10] = m[ 2] - m[ 6];
  d[11] = m[ 3] - m[ 7];
  d[12] = m[12] - m[ 8];
  d[13] = m[13] - m[ 9];
  d[14] = m[14] - m[10];
  d[15] = m[15] - m[11];

  m[ 0] = d[ 0] + d[ 3];
  m[ 1] = d[ 1] + d[ 2];
  m[ 2] = d[ 1] - d[ 2];
  m[ 3] = d[ 0] - d[ 3];
  m[ 4] = d[ 4] + d[ 7];
  m[ 5] = d[ 5] + d[ 6];
  m[ 6] = d[ 5] - d[ 6];
  m[ 7] = d[ 4] - d[ 7];
  m[ 8] = d[ 8] + d[11];
  m[ 9] = d[ 9] + d[10];
  m[10] = d[ 9] - d[10];
  m[11] = d[ 8] - d[11];
  m[12] = d[12] + d[15];
  m[13] = d[13] + d[14];
  m[14] = d[13] - d[14];
  m[15] = d[12] - d[15];

  d[ 0] = m[ 0] + m[ 1];
  d[ 1] = m[ 0] - m[ 1];
  d[ 2] = m[ 2] + m[ 3];
  d[ 3] = m[ 3] - m[ 2];
  d[ 4] = m[ 4] + m[ 5];
  d[ 5] = m[ 4] - m[ 5];
  d[ 6] = m[ 6] + m[ 7];
  d[ 7] = m[ 7] - m[ 6];
  d[ 8] = m[ 8] + m[ 9];
  d[ 9] = m[ 8] - m[ 9];
  d[10] = m[10] + m[11];
  d[11] = m[11] - m[10];
  d[12] = m[12] + m[13];
  d[13] = m[12] - m[13];
  d[14] = m[14] + m[15];
  d[15] = m[15] - m[14];

  for (k=0; k<16; ++k) {
    satd += abs(d[k]);
  }
  satd = ((satd+1)>>1);

  return satd;
}

/**
 * \brief  Calculate SATD between two 8x8 blocks inside bigger arrays.
 */
static unsigned satd_16bit_8x8_general(pixel *piOrg, int32_t iStrideOrg,
                                       pixel *piCur, int32_t iStrideCur)
{
  int32_t k, i, j, jj, sad=0;
  int32_t diff[64], m1[8][8], m2[8][8], m3[8][8];

  for (k = 0; k < 64; k += 8) {
    diff[k+0] = piOrg[0] - piCur[0];
    diff[k+1] = piOrg[1] - piCur[1];
    diff[k+2] = piOrg[2] - piCur[2];
    diff[k+3] = piOrg[3] - piCur[3];
    diff[k+4] = piOrg[4] - piCur[4];
    diff[k+5] = piOrg[5] - piCur[5];
    diff[k+6] = piOrg[6] - piCur[6];
    diff[k+7] = piOrg[7] - piCur[7];

    piCur += iStrideCur;
    piOrg += iStrideOrg;
  }

  // horizontal
  for (j = 0; j < 8; ++j) {
    jj = j << 3;
    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  // vertical
  for (i = 0; i < 8; ++i) {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j) {
      sad += abs(m2[i][j]);
    }
  }

  sad = (sad + 2) >> 2;

  return sad;
}

// Function macro for defining hadamard calculating functions
// for fixed size blocks. They calculate hadamard for integer
// multiples of 8x8 with the 8x8 hadamard function.
#define SATD_NXN(n, pixel_type, suffix) \
  static unsigned satd_ ## suffix ## _ ## n ## x ## n( \
                  pixel_type *block1, pixel_type *block2) \
  { \
    unsigned y, x; \
    unsigned sum = 0; \
    for (y = 0; y < (n); y += 8) { \
      unsigned row = y * (n); \
      for (x = 0; x < (n); x += 8) { \
        sum += satd_16bit_8x8_general(&block1[row + x], (n), &block2[row + x], (n)); \
      } \
    } \
    return sum; \
    }

// These macros define sadt_16bit_NxN for N = 8, 16, 32, 64
SATD_NXN(8, pixel, 16bit)
SATD_NXN(16, pixel, 16bit)
SATD_NXN(32, pixel, 16bit)
SATD_NXN(64, pixel, 16bit)

// Function macro for defining SAD calculating functions
// for fixed size blocks.
#define SAD_NXN(n, pixel_type, suffix) \
  static unsigned sad_ ## suffix ## _ ##  n ## x ## n( \
                  pixel_type *block1, pixel_type *block2) \
  { \
    unsigned x, y, row; \
    unsigned sum = 0; \
    for(y = 0; y < (n); y++) { \
      row = y * (n); \
      for (x = 0; x < (n); ++x) { \
        sum += abs(block1[row + x] - block2[row + x]); \
      } \
    } \
    return sum; \
  }

// These macros define sad_16bit_nxn functions for n = 4, 8, 16, 32, 64
// with function signatures of cost_16bit_nxn_func.
// They are used through get_sad_16bit_nxn_func.
SAD_NXN(4, pixel, 16bit)
SAD_NXN(8, pixel, 16bit)
SAD_NXN(16, pixel, 16bit)
SAD_NXN(32, pixel, 16bit)
SAD_NXN(64, pixel, 16bit)

/**
 * \brief  Get a function that calculates SATD for NxN block.
 *
 * \param n  Width of the region for which SATD is calculated.
 *
 * \returns  Pointer to cost_16bit_nxn_func.
 */
cost_16bit_nxn_func get_satd_16bit_nxn_func(unsigned n)
{
  switch (n) {
  case 4:
    return &satd_16bit_4x4;
  case 8:
    return &satd_16bit_8x8;
  case 16:
    return &satd_16bit_16x16;
  case 32:
    return &satd_16bit_32x32;
  case 64:
    return &satd_16bit_64x64;
  default:
    return NULL;
    }
  }

/**
 * \brief  Get a function that calculates SAD for NxN block.
 *
 * \param n  Width of the region for which SAD is calculated.
 *
 * \returns  Pointer to cost_16bit_nxn_func.
 */
cost_16bit_nxn_func get_sad_16bit_nxn_func(unsigned n)
  {
  switch (n) {
  case 4:
    return &sad_16bit_4x4;
  case 8:
    return &sad_16bit_8x8;
  case 16:
    return &sad_16bit_16x16;
  case 32:
    return &sad_16bit_32x32;
  case 64:
    return &sad_16bit_64x64;
  default:
    return NULL;
  }
}

#if 0
/**
 * \brief  Calculate SATD for NxN block of size N.
 *
 * \param block1  Start of the first block.
 * \param block2  Start of the second block.
 * \param n       Width of the region for which SAD is calculated.
 *
 * \returns       Sum of Absolute Transformed Differences (SATD)
 */
static unsigned satd_nxn_16bit(pixel *block1, pixel *block2, unsigned n)
{
  cost_16bit_nxn_func sad_func = get_satd_16bit_nxn_func(n);
  return sad_func(block1, block2);
  }

/**
 * \brief Calculate SAD for NxN block of size N.
 *
 * \param block1  Start of the first block.
 * \param block2  Start of the second block.
 * \param n       Width of the region for which SAD is calculated.
 *
 * \returns       Sum of Absolute Differences
 */
static unsigned sad_nxn_16bit(pixel *block1, pixel *block2, unsigned n)
{
  cost_16bit_nxn_func sad_func = get_sad_16bit_nxn_func(n);
  if (sad_func) {
    return sad_func(block1, block2);
  } else {
    unsigned row, x;
    unsigned sum = 0;
    for (row = 0; row < n; row += n) {
      for (x = 0; x < n; ++x) {
        sum += abs(block1[row + x] - block2[row + x]);
  }
    }
  return sum;
}
}
#endif

/**
 * \brief Diagonally interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned cor_sad(const pixel *pic_data, const pixel *ref_data,
                        int block_width, int block_height, unsigned pic_width)
{
  pixel ref = *ref_data;
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_width + x] - ref);
    }
  }

  return sad;
}

/**
 * \brief Vertically interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned ver_sad(const pixel *pic_data, const pixel *ref_data,
                        int block_width, int block_height, unsigned pic_width)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_width + x] - ref_data[x]);
    }
  }

  return sad;
}

/**
 * \brief Horizontally interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width   Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned hor_sad(const pixel *pic_data, const pixel *ref_data,
                        int block_width, int block_height, unsigned pic_width, unsigned ref_width)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_width + x] - ref_data[y * ref_width]);
    }
  }

  return sad;
}


#if defined(__SSE2__)
#include "picture-sse2.c"
#elif defined(__ALTIVEC__)
#include "picture-altivec.c"
#else
//Generic implementations
/**
 * \brief Calculate Sum of Absolute Differences (SAD)
 *
 * Calculate Sum of Absolute Differences (SAD) between two rectangular regions
 * located in arbitrary points in the picture.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param stride  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned reg_sad(const pixel * const data1, const pixel * const data2,
                        const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  int y, x;
  unsigned sad = 0;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }

  return sad;
}
#endif

/**
 * \brief  Handle special cases of comparing blocks that are not completely
 *         inside the frame.
 *
 * \param pic  First frame.
 * \param ref  Second frame.
 * \param pic_x  X coordinate of the first block.
 * \param pic_y  Y coordinate of the first block.
 * \param ref_x  X coordinate of the second block.
 * \param ref_y  Y coordinate of the second block.
 * \param block_width  Width of the blocks.
 * \param block_height  Height of the blocks.
 */
static unsigned interpolated_sad(const picture *pic, const picture *ref,
                                 int pic_x, int pic_y, int ref_x, int ref_y,
                                 int block_width, int block_height)
{
  pixel *pic_data, *ref_data;

  int left, right, top, bottom;
  int result = 0;

  // Change the movement vector to point right next to the frame. This doesn't
  // affect the result but removes some special cases.
  if (ref_x > ref->width)            ref_x = ref->width;
  if (ref_y > ref->height)           ref_y = ref->height;
  if (ref_x + block_width < 0)  ref_x = -block_width;
  if (ref_y + block_height < 0) ref_y = -block_height;

  // These are the number of pixels by how far the movement vector points
  // outside the frame. They are always >= 0. If all of them are 0, the
  // movement vector doesn't point outside the frame.
  left   = (ref_x < 0) ? -ref_x : 0;
  top    = (ref_y < 0) ? -ref_y : 0;
  right  = (ref_x + block_width  > ref->width)  ? ref_x + block_width  - ref->width  : 0;
  bottom = (ref_y + block_height > ref->height) ? ref_y + block_height - ref->height : 0;

  // Center picture to the current block and reference to the point where
  // movement vector is pointing to. That point might be outside the buffer,
  // but that is ok because we project the movement vector to the buffer
  // before dereferencing the pointer.
  pic_data = &pic->y_data[pic_y * pic->width + pic_x];
  ref_data = &ref->y_data[ref_y * ref->width + ref_x];

  // The handling of movement vectors that point outside the picture is done
  // in the following way.
  // - Correct the index of ref_data so that it points to the top-left
  //   of the area we want to compare against.
  // - Correct the index of pic_data to point inside the current block, so
  //   that we compare the right part of the block to the ref_data.
  // - Reduce block_width and block_height so that the the size of the area
  //   being compared is correct.
  if (top && left) {
    result += cor_sad(pic_data,
                      &ref_data[top * ref->width + left],
                      left, top, pic->width);
    result += ver_sad(&pic_data[left],
                      &ref_data[top * ref->width + left],
                      block_width - left, top, pic->width);
    result += hor_sad(&pic_data[top * pic->width],
                      &ref_data[top * ref->width + left],
                      left, block_height - top, pic->width, ref->width);
    result += reg_sad(&pic_data[top * pic->width + left],
                      &ref_data[top * ref->width + left],
                      block_width - left, block_height - top, pic->width, ref->width);
  } else if (top && right) {
    result += ver_sad(pic_data,
                      &ref_data[top * ref->width],
                      block_width - right, top, pic->width);
    result += cor_sad(&pic_data[block_width - right],
                      &ref_data[top * ref->width + (block_width - right - 1)],
                      right, top, pic->width);
    result += reg_sad(&pic_data[top * pic->width],
                      &ref_data[top * ref->width],
                      block_width - right, block_height - top, pic->width, ref->width);
    result += hor_sad(&pic_data[top * pic->width + (block_width - right)],
                      &ref_data[top * ref->width + (block_width - right - 1)],
                      right, block_height - top, pic->width, ref->width);
  } else if (bottom && left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height - bottom, pic->width, ref->width);
    result += reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height - bottom, pic->width, ref->width);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->width],
                      &ref_data[(block_height - bottom - 1) * ref->width + left],
                      left, bottom, pic->width);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->width + left],
                      &ref_data[(block_height - bottom - 1) * ref->width + left],
                      block_width - left, bottom, pic->width);
  } else if (bottom && right) {
    result += reg_sad(pic_data,
                      ref_data,
                      block_width - right, block_height - bottom, pic->width, ref->width);
    result += hor_sad(&pic_data[block_width - right],
                      &ref_data[block_width - right - 1],
                      right, block_height - bottom, pic->width, ref->width);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->width],
                      &ref_data[(block_height - bottom - 1) * ref->width],
                      block_width - right, bottom, pic->width);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->width + block_width - right],
                      &ref_data[(block_height - bottom - 1) * ref->width + block_width - right - 1],
                      right, bottom, pic->width);
  } else if (top) {
    result += ver_sad(pic_data,
                      &ref_data[top * ref->width],
                      block_width, top, pic->width);
    result += reg_sad(&pic_data[top * pic->width],
                      &ref_data[top * ref->width],
                      block_width, block_height - top, pic->width, ref->width);
  } else if (bottom) {
    result += reg_sad(pic_data,
                      ref_data,
                      block_width, block_height - bottom, pic->width, ref->width);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->width],
                      &ref_data[(block_height - bottom - 1) * ref->width],
                      block_width, bottom, pic->width);
  } else if (left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height, pic->width, ref->width);
    result += reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height, pic->width, ref->width);
  } else if (right) {
    result += reg_sad(pic_data,
                      ref_data,
                      block_width - right, block_height, pic->width, ref->width);
    result += hor_sad(&pic_data[block_width - right],
                      &ref_data[block_width - right - 1],
                      right, block_height, pic->width, ref->width);
  } else {
    result += reg_sad(pic_data, ref_data, block_width, block_height, pic->width, ref->width);
  }

  return result;
}

/**
 * \brief  Get Sum of Absolute Differences (SAD) between two blocks in two
 *         different frames.
 *
 * \param pic  First frame.
 * \param ref  Second frame.
 * \param pic_x  X coordinate of the first block.
 * \param pic_y  Y coordinate of the first block.
 * \param ref_x  X coordinate of the second block.
 * \param ref_y  Y coordinate of the second block.
 * \param block_width  Width of the blocks.
 * \param block_height  Height of the blocks.
 */
unsigned calc_sad(const picture *pic, const picture *ref,
                  int pic_x, int pic_y, int ref_x, int ref_y,
                  int block_width, int block_height)
{
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);
  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const pixel *pic_data = &pic->y_data[pic_y * pic->width + pic_x];
    const pixel *ref_data = &ref->y_data[ref_y * ref->width + ref_x];
    return reg_sad(pic_data, ref_data, block_width, block_height, pic->width, ref->width);
  } else {
    // Call a routine that knows how to interpolate pixels outside the frame.
    return interpolated_sad(pic, ref, pic_x, pic_y, ref_x, ref_y, block_width, block_height);
  }
}
