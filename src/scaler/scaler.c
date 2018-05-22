/*****************************************************************************
* This file is part of Kvazaar HEVC encoder.
*
* Copyright (C) 2013-2015 Tampere University of Technology and others (see
* COPYING file).
*
* Kvazaar is free software: you can redistribute it and/or modify it under
* the terms of the GNU Lesser General Public License as published by the
* Free Software Foundation; either version 2.1 of the License, or (at your
* option) any later version.
*
* Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along
* with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

#include "scaler.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

#define SHIFT(x,y) (((y) < 0) ? ((x)>>(-(y))) : ((x)<<(y)))

//Define filters for scaling operations
//Values from SHM
static const int filter16[8][16][12] = {
  // ratio <= 20/19
  {
    {0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 2, -6, 127, 7, -2, 0, 0, 0, 0},
    {0, 0, 0, 3, -12, 125, 16, -5, 1, 0, 0, 0},
    {0, 0, 0, 4, -16, 120, 26, -7, 1, 0, 0, 0},
    {0, 0, 0, 5, -18, 114, 36, -10, 1, 0, 0, 0},
    {0, 0, 0, 5, -20, 107, 46, -12, 2, 0, 0, 0},
    {0, 0, 0, 5, -21, 99, 57, -15, 3, 0, 0, 0},
    {0, 0, 0, 5, -20, 89, 68, -18, 4, 0, 0, 0},
    {0, 0, 0, 4, -19, 79, 79, -19, 4, 0, 0, 0},
    {0, 0, 0, 4, -18, 68, 89, -20, 5, 0, 0, 0},
    {0, 0, 0, 3, -15, 57, 99, -21, 5, 0, 0, 0},
    {0, 0, 0, 2, -12, 46, 107, -20, 5, 0, 0, 0},
    {0, 0, 0, 1, -10, 36, 114, -18, 5, 0, 0, 0},
    {0, 0, 0, 1, -7, 26, 120, -16, 4, 0, 0, 0},
    {0, 0, 0, 1, -5, 16, 125, -12, 3, 0, 0, 0},
    {0, 0, 0, 0, -2, 7, 127, -6, 2, 0, 0, 0}
  },
  // 20/19 < ratio <= 5/4
  {
    {0, 2, 0, -14, 33, 86, 33, -14, 0, 2, 0, 0},
    {0, 1, 1, -14, 29, 85, 38, -13, -1, 2, 0, 0},
    {0, 1, 2, -14, 24, 84, 43, -12, -2, 2, 0, 0},
    {0, 1, 2, -13, 19, 83, 48, -11, -3, 2, 0, 0},
    {0, 0, 3, -13, 15, 81, 53, -10, -4, 3, 0, 0},
    {0, 0, 3, -12, 11, 79, 57, -8, -5, 3, 0, 0},
    {0, 0, 3, -11, 7, 76, 62, -5, -7, 3, 0, 0},
    {0, 0, 3, -10, 3, 73, 65, -2, -7, 3, 0, 0},
    {0, 0, 3, -9, 0, 70, 70, 0, -9, 3, 0, 0},
    {0, 0, 3, -7, -2, 65, 73, 3, -10, 3, 0, 0},
    {0, 0, 3, -7, -5, 62, 76, 7, -11, 3, 0, 0},
    {0, 0, 3, -5, -8, 57, 79, 11, -12, 3, 0, 0},
    {0, 0, 3, -4, -10, 53, 81, 15, -13, 3, 0, 0},
    {0, 0, 2, -3, -11, 48, 83, 19, -13, 2, 1, 0},
    {0, 0, 2, -2, -12, 43, 84, 24, -14, 2, 1, 0},
    {0, 0, 2, -1, -13, 38, 85, 29, -14, 1, 1, 0}
  },
  // 5/4 < ratio <= 5/3
  {
    {0, 5, -6, -10, 37, 76, 37, -10, -6, 5, 0, 0},
    {0, 5, -4, -11, 33, 76, 40, -9, -7, 5, 0, 0},
    {-1, 5, -3, -12, 29, 75, 45, -7, -8, 5, 0, 0},
    {-1, 4, -2, -13, 25, 75, 48, -5, -9, 5, 1, 0},
    {-1, 4, -1, -13, 22, 73, 52, -3, -10, 4, 1, 0},
    {-1, 4, 0, -13, 18, 72, 55, -1, -11, 4, 2, -1},
    {-1, 4, 1, -13, 14, 70, 59, 2, -12, 3, 2, -1},
    {-1, 3, 1, -13, 11, 68, 62, 5, -12, 3, 2, -1},
    {-1, 3, 2, -13, 8, 65, 65, 8, -13, 2, 3, -1},
    {-1, 2, 3, -12, 5, 62, 68, 11, -13, 1, 3, -1},
    {-1, 2, 3, -12, 2, 59, 70, 14, -13, 1, 4, -1},
    {-1, 2, 4, -11, -1, 55, 72, 18, -13, 0, 4, -1},
    {0, 1, 4, -10, -3, 52, 73, 22, -13, -1, 4, -1},
    {0, 1, 5, -9, -5, 48, 75, 25, -13, -2, 4, -1},
    {0, 0, 5, -8, -7, 45, 75, 29, -12, -3, 5, -1},
    {0, 0, 5, -7, -9, 40, 76, 33, -11, -4, 5, 0},
  },
  // 5/3 < ratio <= 2
  {
    {2, -3, -9, 6, 39, 58, 39, 6, -9, -3, 2, 0},
    {2, -3, -9, 4, 38, 58, 43, 7, -9, -4, 1, 0},
    {2, -2, -9, 2, 35, 58, 44, 9, -8, -4, 1, 0},
    {1, -2, -9, 1, 34, 58, 46, 11, -8, -5, 1, 0},
    {1, -1, -8, -1, 31, 57, 47, 13, -7, -5, 1, 0},
    {1, -1, -8, -2, 29, 56, 49, 15, -7, -6, 1, 1},
    {1, 0, -8, -3, 26, 55, 51, 17, -7, -6, 1, 1},
    {1, 0, -7, -4, 24, 54, 52, 19, -6, -7, 1, 1},
    {1, 0, -7, -5, 22, 53, 53, 22, -5, -7, 0, 1},
    {1, 1, -7, -6, 19, 52, 54, 24, -4, -7, 0, 1},
    {1, 1, -6, -7, 17, 51, 55, 26, -3, -8, 0, 1},
    {1, 1, -6, -7, 15, 49, 56, 29, -2, -8, -1, 1},
    {0, 1, -5, -7, 13, 47, 57, 31, -1, -8, -1, 1},
    {0, 1, -5, -8, 11, 46, 58, 34, 1, -9, -2, 1},
    {0, 1, -4, -8, 9, 44, 58, 35, 2, -9, -2, 2},
    {0, 1, -4, -9, 7, 43, 58, 38, 4, -9, -3, 2},
  },
  // 2 < ratio <= 5/2
  {
    {-2, -7, 0, 17, 35, 43, 35, 17, 0, -7, -5, 2},
    {-2, -7, -1, 16, 34, 43, 36, 18, 1, -7, -5, 2},
    {-1, -7, -1, 14, 33, 43, 36, 19, 1, -6, -5, 2},
    {-1, -7, -2, 13, 32, 42, 37, 20, 3, -6, -5, 2},
    {0, -7, -3, 12, 31, 42, 38, 21, 3, -6, -5, 2},
    {0, -7, -3, 11, 30, 42, 39, 23, 4, -6, -6, 1},
    {0, -7, -4, 10, 29, 42, 40, 24, 5, -6, -6, 1},
    {1, -7, -4, 9, 27, 41, 40, 25, 6, -5, -6, 1},
    {1, -6, -5, 7, 26, 41, 41, 26, 7, -5, -6, 1},
    {1, -6, -5, 6, 25, 40, 41, 27, 9, -4, -7, 1},
    {1, -6, -6, 5, 24, 40, 42, 29, 10, -4, -7, 0},
    {1, -6, -6, 4, 23, 39, 42, 30, 11, -3, -7, 0},
    {2, -5, -6, 3, 21, 38, 42, 31, 12, -3, -7, 0},
    {2, -5, -6, 3, 20, 37, 42, 32, 13, -2, -7, -1},
    {2, -5, -6, 1, 19, 36, 43, 33, 14, -1, -7, -1},
    {2, -5, -7, 1, 18, 36, 43, 34, 16, -1, -7, -2}
  },
  // 5/2 < ratio <= 20/7
  {
    {-6, -3, 5, 19, 31, 36, 31, 19, 5, -3, -6, 0},
    {-6, -4, 4, 18, 31, 37, 32, 20, 6, -3, -6, -1},
    {-6, -4, 4, 17, 30, 36, 33, 21, 7, -3, -6, -1},
    {-5, -5, 3, 16, 30, 36, 33, 22, 8, -2, -6, -2},
    {-5, -5, 2, 15, 29, 36, 34, 23, 9, -2, -6, -2},
    {-5, -5, 2, 15, 28, 36, 34, 24, 10, -2, -6, -3},
    {-4, -5, 1, 14, 27, 36, 35, 24, 10, -1, -6, -3},
    {-4, -5, 0, 13, 26, 35, 35, 25, 11, 0, -5, -3},
    {-4, -6, 0, 12, 26, 36, 36, 26, 12, 0, -6, -4},
    {-3, -5, 0, 11, 25, 35, 35, 26, 13, 0, -5, -4},
    {-3, -6, -1, 10, 24, 35, 36, 27, 14, 1, -5, -4},
    {-3, -6, -2, 10, 24, 34, 36, 28, 15, 2, -5, -5},
    {-2, -6, -2, 9, 23, 34, 36, 29, 15, 2, -5, -5},
    {-2, -6, -2, 8, 22, 33, 36, 30, 16, 3, -5, -5},
    {-1, -6, -3, 7, 21, 33, 36, 30, 17, 4, -4, -6},
    {-1, -6, -3, 6, 20, 32, 37, 31, 18, 4, -4, -6}
  },
  // 20/7 < ratio <= 15/4
  {
    {-9, 0, 9, 20, 28, 32, 28, 20, 9, 0, -9, 0},
    {-9, 0, 8, 19, 28, 32, 29, 20, 10, 0, -4, -5},
    {-9, -1, 8, 18, 28, 32, 29, 21, 10, 1, -4, -5},
    {-9, -1, 7, 18, 27, 32, 30, 22, 11, 1, -4, -6},
    {-8, -2, 6, 17, 27, 32, 30, 22, 12, 2, -4, -6},
    {-8, -2, 6, 16, 26, 32, 31, 23, 12, 2, -4, -6},
    {-8, -2, 5, 16, 26, 31, 31, 23, 13, 3, -3, -7},
    {-8, -3, 5, 15, 25, 31, 31, 24, 14, 4, -3, -7},
    {-7, -3, 4, 14, 25, 31, 31, 25, 14, 4, -3, -7},
    {-7, -3, 4, 14, 24, 31, 31, 25, 15, 5, -3, -8},
    {-7, -3, 3, 13, 23, 31, 31, 26, 16, 5, -2, -8},
    {-6, -4, 2, 12, 23, 31, 32, 26, 16, 6, -2, -8},
    {-6, -4, 2, 12, 22, 30, 32, 27, 17, 6, -2, -8},
    {-6, -4, 1, 11, 22, 30, 32, 27, 18, 7, -1, -9},
    {-5, -4, 1, 10, 21, 29, 32, 28, 18, 8, -1, -9},
    {-5, -4, 0, 10, 20, 29, 32, 28, 19, 8, 0, -9}
  },
  // ratio > 15/4
  {
    {-8, 7, 13, 18, 22, 24, 22, 18, 13, 7, 2, -10},
    {-8, 7, 13, 18, 22, 23, 22, 19, 13, 7, 2, -10},
    {-8, 6, 12, 18, 22, 23, 22, 19, 14, 8, 2, -10},
    {-9, 6, 12, 17, 22, 23, 23, 19, 14, 8, 3, -10},
    {-9, 6, 12, 17, 21, 23, 23, 19, 14, 9, 3, -10},
    {-9, 5, 11, 17, 21, 23, 23, 20, 15, 9, 3, -10},
    {-9, 5, 11, 16, 21, 23, 23, 20, 15, 9, 4, -10},
    {-9, 5, 10, 16, 21, 23, 23, 20, 15, 10, 4, -10},
    {-10, 5, 10, 16, 20, 23, 23, 20, 16, 10, 5, -10},
    {-10, 4, 10, 15, 20, 23, 23, 21, 16, 10, 5, -9},
    {-10, 4, 9, 15, 20, 23, 23, 21, 16, 11, 5, -9},
    {-10, 3, 9, 15, 20, 23, 23, 21, 17, 11, 5, -9},
    {-10, 3, 9, 14, 19, 23, 23, 21, 17, 12, 6, -9},
    {-10, 3, 8, 14, 19, 23, 23, 22, 17, 12, 6, -9},
    {-10, 2, 8, 14, 19, 22, 23, 22, 18, 12, 6, -8},
    {-10, 2, 7, 13, 19, 22, 23, 22, 18, 13, 7, -8}
  }
};

static const int lumaUpFilter[16][8] = {
  {0, 0, 0, 64, 0, 0, 0, 0},
  {0, 1, -3, 63, 4, -2, 1, 0},
  {-1, 2, -5, 62, 8, -3, 1, 0},
  {-1, 3, -8, 60, 13, -4, 1, 0},
  {-1, 4, -10, 58, 17, -5, 1, 0},
  {-1, 4, -11, 52, 26, -8, 3, -1},
  {-1, 3, -9, 47, 31, -10, 4, -1},
  {-1, 4, -11, 45, 34, -10, 4, -1},
  {-1, 4, -11, 40, 40, -11, 4, -1},
  {-1, 4, -10, 34, 45, -11, 4, -1},
  {-1, 4, -10, 31, 47, -9, 3, -1},
  {-1, 3, -8, 26, 52, -11, 4, -1},
  {0, 1, -5, 17, 58, -10, 4, -1},
  {0, 1, -4, 13, 60, -8, 3, -1},
  {0, 1, -3, 8, 62, -5, 2, -1},
  {0, 1, -2, 4, 63, -3, 1, 0}
};

static const int chromaUpFilter[16][4] = {
  {0, 64, 0, 0},
  {-2, 62, 4, 0},
  {-2, 58, 10, -2},
  {-4, 56, 14, -2},
  {-4, 54, 16, -2},
  {-6, 52, 20, -2},
  {-6, 46, 28, -4},
  {-4, 42, 30, -4},
  {-4, 36, 36, -4},
  {-4, 30, 42, -4},
  {-4, 28, 46, -6},
  {-2, 20, 52, -6},
  {-2, 16, 54, -4},
  {-2, 14, 56, -4},
  {-2, 10, 58, -2},
  {0, 4, 62, -2}
};

//Used for clipping values
static int clip(int val, int min, int max)
{
  if (val <= min)
    return min;
  if (val >= max)
    return max;

  return val;
}

pic_buffer_t* kvz_newPictureBuffer(int width, int height, int has_tmp_row)
{
  pic_buffer_t* buffer = (pic_buffer_t*)malloc(sizeof(pic_buffer_t));
  if (buffer == NULL) {
    return NULL; //TODO: Add error message?
  }

  //Allocate enough memory to fit a width-by-height picture
  buffer->data = (pic_data_t*)malloc(sizeof(pic_data_t) * width * height);

  buffer->width = width;
  buffer->height = height;

  //Initialize tmp_row or set as NULL
  if (has_tmp_row) {
    //Use max dim for size
    int max_dim = MAX(width, height);
    buffer->tmp_row = (pic_data_t*)malloc(sizeof(pic_data_t) * max_dim);
  }
  else {
    buffer->tmp_row = NULL;
  }

  return buffer;
}

yuv_buffer_t* kvz_newYuvBuffer(int width, int height , chroma_format_t format, int has_tmp_row)
{
  yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));
  if (yuv == NULL) {
    return NULL; //TODO: Add error message?
  }
  yuv->format =format;
  yuv->y = kvz_newPictureBuffer(width, height, has_tmp_row);

  int w_factor = 0;
  int h_factor = 0;

  switch (format) {
    case CHROMA_400:
      {
        //No chroma
        width = height = 0;
        break;
      }
    case CHROMA_420:
      {
        w_factor = 1;
        h_factor = 1;
        break;
      }
    case CHROMA_422:
      {
        w_factor = 1;
        break;
      }
    case CHROMA_444:
      {
        break;
      }
    default:
      assert(0);//Unsupported format
  }

  width = width >> w_factor;
  height = height >> h_factor;

  yuv->u = kvz_newPictureBuffer( width, height, has_tmp_row);
  yuv->v = kvz_newPictureBuffer( width, height, has_tmp_row);

  return yuv;
}

// ======================= newPictureBuffer_ ==================================
//TODO: DO something about the lack of overloading?
/**
* \brief Create/Initialize a Picture buffer. Width/height should be the width/height of the data. The caller is responsible for deallocation.
*/
//static pic_buffer_t* newPictureBuffer_double(const double* const data, int width, int height, int has_tmp_row)
//{
//  pic_buffer_t* buffer = kvz_newPictureBuffer(width, height, has_tmp_row);
//
//  //If data is null skip initializing
//  if (data == NULL) return buffer;
//
//  //Initialize buffer
//  for (int i = width * height - 1; i >= 0; i--) {
//    buffer->data[i] = (int)data[i];
//  }
//
//  return buffer;
//}
//
///**
//* \brief Create/Initialize a Picture buffer. Width/height should be the width/height of the data. The caller is responsible for deallocation.
//*/
//static pic_buffer_t* newPictureBuffer_uint8(const uint8_t* const data, int width, int height, int has_tmp_row)
//{
//  pic_buffer_t* buffer = kvz_newPictureBuffer(width, height, has_tmp_row);
//
//  //If data is null skip initializing
//  if (data == NULL) return buffer;
//
//  //Initialize buffer
//  for (int i = width * height - 1; i >= 0; i--) {
//    buffer->data[i] = (int)data[i];
//  }
//
//  return buffer;
//}

/**
* \brief Create/Initialize a Picture buffer. Width/height should be the width/height of the final buffer. Stride should be the width of the input (padded image). The caller is responsible for deallocation
*/
static pic_buffer_t* newPictureBuffer_padded_uint8(const uint8_t* const data, int width, int height, int stride, int has_tmp_row)
{
  pic_buffer_t* buffer = kvz_newPictureBuffer(width, height, has_tmp_row);

  //If data is null skip initializing
  if (data == NULL) return buffer;

  //Initialize buffer
  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {
      buffer->data[col + row*width] =  (pic_data_t)data[col + row*stride];
    }
  }

  return buffer;
}

// ==============================================================================
/**
 * \brief Deallocate a picture buffer.
 */
static void deallocatePictureBuffer(pic_buffer_t* buffer)
{
  if (buffer != NULL) {
    free(buffer->data);
    free(buffer->tmp_row);
  }
  free(buffer);
}

/**
 * \brief Copies data from one buffer to the other.
 * \param src is the source buffer
 * \param dst is the destination buffer
 * \param fill signals if the inds in dst not overlapped by src should be filled
*    with values adjacent to the said index.
 */
static void copyPictureBuffer(const pic_buffer_t* const src, const pic_buffer_t* const dst, int fill)
{
  //TODO: add checks. Check if fill is necessary
  //max_dim_* is chosen so that no over indexing happenes (src/dst)
  //min_dim_* is chosen so that no over indexing happenes (src), but all inds in dst get a value
  int max_dim_x = fill ? dst->width : MIN(src->width, dst->width);
  int max_dim_y = fill ? dst->height : MIN(src->height, dst->height);
  int min_dim_x = fill ? src->width : max_dim_x;
  int min_dim_y = fill ? src->height : max_dim_y;

  int dst_row = 0;
  int src_row = 0;

  //Copy loop
  for (int i = 0; i < max_dim_y; i++) {
    if (i < min_dim_y) {
      for (int j = 0; j < max_dim_x; j++) {
        //If outside min x, copy adjacent value.
        dst->data[dst_row + j] = (j < min_dim_x) ? src->data[src_row + j] : dst->data[dst_row + j - 1];
      }
    }
    //Handle extra rows if needed
    else {
      for (int j = 0; j < max_dim_x; j++) {
        dst->data[dst_row + j] = dst->data[dst_row + j - dst->width];
      }
    }
    dst_row += dst->width; //switch to the next row
    src_row += src->width; //switch to the next row
  }
}

/**
* \brief Copies data from one buffer to the other.
* \param src is the source buffer
* \param dst is the destination buffer
* \param block_x is the x-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_y is the y-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_width is the width for the sub-block (needs to be valid for both buffers).
* \param block_height is the height for the sub-block (needs to be valid for both buffers).
*/
static void copyPictureBufferBlock(const pic_buffer_t* const src, const pic_buffer_t* const dst, const int block_x, const int block_y, const int block_width, const int block_height )
{
  for( int row = block_y; row < block_height; row++ ){
    for( int col = block_x; col < block_width; col++ ){
      dst->data[col + row*dst->width] = src->data[col + row*src->width];
    }
  }
}

/**
* \brief Copies data from one buffer to the other.
* \param src is the source buffer
* \param dst is the destination buffer
* \param fill signals if the inds in dst not overlapped by src should be filled
*    with values adjacent to the said index.
*/
static void copyYuvBuffer(const yuv_buffer_t* const src, const yuv_buffer_t* const dst, int fill)
{
  copyPictureBuffer(src->y, dst->y, fill);
  copyPictureBuffer(src->u, dst->u, fill);
  copyPictureBuffer(src->v, dst->v, fill);
}

/**
* \brief Copies data from a sub-block from one buffer to the other.
* \param src is the source buffer
* \param dst is the destination buffer
* \param block_x is the x-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_y is the y-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_width is the width for the sub-block (needs to be valid for both buffers).
* \param block_height is the height for the sub-block (needs to be valid for both buffers).
* \param w_factor is how much chroma sizes are scaled (width).
* \param h_factor is how much chroma sizes are scaled (heigth).
*/
static void copyYuvBufferBlock(const yuv_buffer_t* const src, const yuv_buffer_t* const dst, const int block_x, const int block_y, const int block_width, const int block_height, const int w_factor, const int h_factor)
{
  copyPictureBufferBlock(src->y, dst->y, block_x, block_y, block_width, block_height);
  copyPictureBufferBlock(src->u, dst->u, SHIFT(block_x, w_factor), SHIFT(block_y, h_factor), SHIFT(block_width, w_factor), SHIFT(block_height, h_factor));
  copyPictureBufferBlock(src->v, dst->v, SHIFT(block_x, w_factor), SHIFT(block_y, h_factor), SHIFT(block_width, w_factor), SHIFT(block_height, h_factor));
}

// ======================= newYuvBuffer_ ==================================
//static yuv_buffer_t* newYuvBuffer_double(const double* const y_data, const double* const u_data, const double* const v_data, int width, int height, chroma_format_t format, int has_tmp_row)
//{
//  yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));
//  yuv->format = format;
//
//  //Allocate y pic_buffer
//  yuv->y = newPictureBuffer_double(y_data, width, height, has_tmp_row);
//
//  //Allocate u and v buffers
//  int w_factor = 0;
//  int h_factor = 0;
//
//  switch (format) {
//    case CHROMA_400:
//      {
//        //No chroma
//        width = height = 0;
//        break;
//      }
//    case CHROMA_420:
//      {
//        w_factor = 1;
//        h_factor = 1;
//        break;
//      }
//    case CHROMA_422:
//      {
//        w_factor = 1;
//        break;
//      }
//    case CHROMA_444:
//      {
//        break;
//      }
//    default:
//      assert(0);//Unsupported format
//  }
//
//  width = width >> w_factor;
//  height = height >> h_factor;
//  yuv->u = newPictureBuffer_double(u_data, width, height, has_tmp_row);
//  yuv->v = newPictureBuffer_double(v_data, width, height, has_tmp_row);
//
//  return yuv;
//}
//
//static yuv_buffer_t* newYuvBuffer_uint8(const uint8_t* const y_data, const uint8_t* const u_data, const uint8_t* const v_data, int width, int height, chroma_format_t format, int has_tmp_row)
//{
//  yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));
//
//  //Allocate y pic_buffer
//  yuv->y = newPictureBuffer_uint8(y_data, width, height, has_tmp_row);
//  yuv->format = format;
//
//  //Allocate u and v buffers
//  int w_factor = 0;
//  int h_factor = 0;
//
//  switch (format) {
//    case CHROMA_400: {
//      //No chroma
//      width = height = 0;
//      break;
//    }
//    case CHROMA_420: {
//      w_factor = 1;
//      h_factor = 1;
//      break;
//    }
//    case CHROMA_422: {
//      w_factor = 1;
//      break;
//    }
//    case CHROMA_444: {
//      break;
//    }
//    default:
//      assert(0);//Unsupported format
//  }
//
//  width = width >> w_factor;
//  height = height >> h_factor;
//  yuv->u = newPictureBuffer_uint8(u_data, width, height, has_tmp_row);
//  yuv->v = newPictureBuffer_uint8(v_data, width, height, has_tmp_row);
//
//  return yuv;
//}

yuv_buffer_t* kvz_newYuvBuffer_padded_uint8(const uint8_t* const y_data, const uint8_t* const u_data, const uint8_t* const v_data, int width, int height, int stride, chroma_format_t format, int has_tmp_row)
{
  yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));
  if (yuv == NULL) {
    return NULL; //TODO: Add error message?
  }
  //Allocate y pic_buffer
  yuv->y = newPictureBuffer_padded_uint8(y_data, width, height, stride, has_tmp_row);

  //Allocate u and v buffers
  int w_factor = 0;
  int h_factor = 0;

  switch (format) {
    case CHROMA_400: {
      //No chroma
      width = height = 0;
      break;
    }
    case CHROMA_420: {
      w_factor = 1;
      h_factor = 1;
      break;
    }
    case CHROMA_422: {
      w_factor = 1;
      break;
    }
    case CHROMA_444: {
      break;
    }
    default:
      assert(0);//Unsupported format
  }

  width = width >> w_factor;
  height = height >> h_factor;
  stride = stride >> w_factor;
  yuv->u = newPictureBuffer_padded_uint8(u_data, width, height, stride, has_tmp_row);
  yuv->v = newPictureBuffer_padded_uint8(v_data, width, height, stride, has_tmp_row);

  return yuv;
}

// ==============================================================================

/**
* \brief Clone the given pic buffer
*/
static pic_buffer_t* clonePictureBuffer(const pic_buffer_t* const pic)
{
  pic_buffer_t* ret = malloc(sizeof(pic_buffer_t));
  if (ret == NULL) {
    return NULL; //TODO: Add error message?
  }
  int size = pic->width * pic->height;

  *ret = *pic;
  ret->data = malloc(sizeof(pic_data_t) * size);
  if (ret->data == NULL) {
    free(ret);
    return NULL; //TODO: Add error message?
  }
  memcpy(ret->data, pic->data, size * sizeof(pic_data_t));

  if (pic->tmp_row) {
    int tmp_size = MAX(pic->width, pic->height);
    ret->tmp_row = malloc(sizeof(pic_buffer_t) * tmp_size);
    if (ret->tmp_row == NULL) {
      deallocatePictureBuffer(ret);
      return NULL; //TODO: Add error message?
    }
    memcpy(ret->tmp_row, pic->tmp_row, tmp_size);
  }

  return ret;
}

yuv_buffer_t* kvz_cloneYuvBuffer(const yuv_buffer_t* const yuv)
{
  yuv_buffer_t* ret = malloc(sizeof(yuv_buffer_t));
  if (ret == NULL) {
    return NULL; //TODO: Add error message?
  }
  ret->y = clonePictureBuffer(yuv->y);
  ret->u = clonePictureBuffer(yuv->u);
  ret->v = clonePictureBuffer(yuv->v);

  return ret;
}

void kvz_deallocateYuvBuffer(yuv_buffer_t* yuv)
{
  if (yuv == NULL) return;

  deallocatePictureBuffer(yuv->y);
  deallocatePictureBuffer(yuv->u);
  deallocatePictureBuffer(yuv->v);

  free(yuv);
}


//Helper function for choosing the correct filter
//Returns the size of the filter and the filter param is set to the correct filter
static int getFilter(const int** const filter, int is_upsampling, int is_luma, int phase, int filter_ind)
{
  if (is_upsampling) {
    //Upsampling so use 8- or 4-tap filters
    if (is_luma) {
      *filter = lumaUpFilter[phase];
      return sizeof(lumaUpFilter[0]) / sizeof(lumaUpFilter[0][0]);
    }

    *filter = chromaUpFilter[phase];
    return sizeof(chromaUpFilter[0]) / sizeof(chromaUpFilter[0][0]);
  }

  //Downsampling so use 12-tap filter
  *filter = filter16[filter_ind][phase];
  return (sizeof(filter16[0][0]) / sizeof(filter16[0][0][0]));
}

//Resampling is done here per buffer
static void _resample(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
{
  //TODO: Add cropping etc.

  //Choose best filter to use when downsampling
  //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
  int ver_filter = 0;
  int hor_filter = 0;

  int src_height = param->src_height;
  int src_width = param->src_width;
  int trgt_height = param->trgt_height;//param->rnd_trgt_height;
  int trgt_width = param->trgt_width;//param->rnd_trgt_width;

  if (!is_upscaling) {
    int crop_width = src_width - param->right_offset; //- param->left_offset;
    int crop_height = src_height - param->bottom_offset; //- param->top_offset;

    if (4 * crop_height > 15 * trgt_height)
      ver_filter = 7;
    else if (7 * crop_height > 20 * trgt_height)
      ver_filter = 6;
    else if (2 * crop_height > 5 * trgt_height)
      ver_filter = 5;
    else if (1 * crop_height > 2 * trgt_height)
      ver_filter = 4;
    else if (3 * crop_height > 5 * trgt_height)
      ver_filter = 3;
    else if (4 * crop_height > 5 * trgt_height)
      ver_filter = 2;
    else if (19 * crop_height > 20 * trgt_height)
      ver_filter = 1;

    if (4 * crop_width > 15 * trgt_width)
      hor_filter = 7;
    else if (7 * crop_width > 20 * trgt_width)
      hor_filter = 6;
    else if (2 * crop_width > 5 * trgt_width)
      hor_filter = 5;
    else if (1 * crop_width > 2 * trgt_width)
      hor_filter = 4;
    else if (3 * crop_width > 5 * trgt_width)
      hor_filter = 3;
    else if (4 * crop_width > 5 * trgt_width)
      hor_filter = 2;
    else if (19 * crop_width > 20 * trgt_width)
      hor_filter = 1;
  }

  int shift_x = param->shift_x - 4;
  int shift_y = param->shift_y - 4;

  pic_data_t* tmp_row = buffer->tmp_row;

  // Horizontal downsampling
  for (int i = 0; i < src_height; i++) {
    pic_data_t* src_row = &buffer->data[i * buffer->width];

    for (int j = 0; j < trgt_width; j++) {
      //Calculate reference position in src pic
      int ref_pos_16 = (int)((unsigned int)(j * param->scale_x + param->add_x) >> shift_x)  - param->delta_x;
      int phase = ref_pos_16 & 15;
      int ref_pos = ref_pos_16 >> 4;

      //Choose filter
      const int* filter;
      int size = getFilter(&filter, is_upscaling, is_luma, phase, hor_filter);

      //Apply filter
      tmp_row[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = clip(ref_pos + k - (size >> 1) + 1, 0, src_width - 1);
        tmp_row[j] += filter[k] * src_row[m];
      }
    }
    //Copy tmp row to data
    memcpy(src_row, tmp_row, sizeof(pic_data_t) * trgt_width);
  }

  pic_data_t* tmp_col = tmp_row; //rename for clarity

  // Vertical downsampling
  for (int i = 0; i < trgt_width; i++) {
    pic_data_t* src_col = &buffer->data[i];
    for (int j = 0; j < trgt_height; j++) {
      //Calculate ref pos
      int ref_pos_16 = (int)((unsigned int)(j * param->scale_y + param->add_y) >> shift_y) - param->delta_y;
      int phase = ref_pos_16 & 15;
      int ref_pos = ref_pos_16 >> 4;

      //Choose filter
      const int* filter;
      int size = getFilter(&filter, is_upscaling, is_luma, phase, ver_filter);

      //Apply filter
      tmp_col[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = clip(ref_pos + k - (size >> 1) + 1, 0, src_height - 1);
        tmp_col[j] += filter[k] * src_col[m * buffer->width];
      }
      //TODO: Why? Filter coefs summ up to 128 applied 2x 128*128= 2^14
      //TODO: Why?  Filter coefs summ up to 64 applied 2x 64*64= 2^12
      //Scale values back down
      tmp_col[j] = is_upscaling ? (tmp_col[j] + 2048) >> 12 : (tmp_col[j] + 8192) >> 14;
    }

    //Clip and move to buffer data
    for (int n = 0; n < trgt_height; n++) {
      src_col[n * buffer->width] = clip(tmp_col[n], 0, 255);
    }
  }
}

//Resampling is done here per buffer
static void resample(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
{
  //TODO: Add cropping etc.

  //Choose best filter to use when downsampling
  //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
  int ver_filter = 0;
  int hor_filter = 0;

  int src_width = param->src_width + param->src_padding_x;
  int src_height = param->src_height + param->src_padding_y;
  int trgt_width = param->rnd_trgt_width;
  int trgt_height = param->rnd_trgt_height;

  if (!is_upscaling) {
    int crop_width = src_width - param->right_offset; //- param->left_offset;
    int crop_height = src_height - param->bottom_offset; //- param->top_offset;

    if (4 * crop_height > 15 * trgt_height)
      ver_filter = 7;
    else if (7 * crop_height > 20 * trgt_height)
      ver_filter = 6;
    else if (2 * crop_height > 5 * trgt_height)
      ver_filter = 5;
    else if (1 * crop_height > 2 * trgt_height)
      ver_filter = 4;
    else if (3 * crop_height > 5 * trgt_height)
      ver_filter = 3;
    else if (4 * crop_height > 5 * trgt_height)
      ver_filter = 2;
    else if (19 * crop_height > 20 * trgt_height)
      ver_filter = 1;

    if (4 * crop_width > 15 * trgt_width)
      hor_filter = 7;
    else if (7 * crop_width > 20 * trgt_width)
      hor_filter = 6;
    else if (2 * crop_width > 5 * trgt_width)
      hor_filter = 5;
    else if (1 * crop_width > 2 * trgt_width)
      hor_filter = 4;
    else if (3 * crop_width > 5 * trgt_width)
      hor_filter = 3;
    else if (4 * crop_width > 5 * trgt_width)
      hor_filter = 2;
    else if (19 * crop_width > 20 * trgt_width)
      hor_filter = 1;
  }

  int shift_x = param->shift_x - 4;
  int shift_y = param->shift_y - 4;

  pic_data_t* tmp_row = buffer->tmp_row;

  // Horizontal resampling
  for (int i = 0; i < src_height; i++) {
    pic_data_t* src_row = &buffer->data[i * buffer->width];

    for (int j = 0; j < trgt_width; j++) {
      //Calculate reference position in src pic
      int ref_pos_16 = (int)((unsigned int)(j * param->scale_x + param->add_x) >> shift_x) - param->delta_x;
      int phase = ref_pos_16 & 15;
      int ref_pos = ref_pos_16 >> 4;

      //Choose filter
      const int* filter;
      int size = getFilter(&filter, is_upscaling, is_luma, phase, hor_filter);

      //Apply filter
      tmp_row[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = clip(ref_pos + k - (size >> 1) + 1, 0, src_width - 1);
        tmp_row[j] += filter[k] * src_row[m];
      }
    }
    //Copy tmp row to data
    memcpy(src_row, tmp_row, sizeof(pic_data_t) * trgt_width);
  }

  pic_data_t* tmp_col = tmp_row; //rename for clarity

  // Vertical resampling
  for (int i = 0; i < trgt_width; i++) {
    pic_data_t* src_col = &buffer->data[i];
    for (int j = 0; j < trgt_height; j++) {
      //Calculate ref pos
      int ref_pos_16 = (int)((unsigned int)(j * param->scale_y + param->add_y) >> shift_y) - param->delta_y;
      int phase = ref_pos_16 & 15;
      int ref_pos = ref_pos_16 >> 4;

      //Choose filter
      const int* filter;
      int size = getFilter(&filter, is_upscaling, is_luma, phase, ver_filter);

      //Apply filter
      tmp_col[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = clip(ref_pos + k - (size >> 1) + 1, 0, src_height - 1);
        tmp_col[j] += filter[k] * src_col[m * buffer->width];
      }
      //TODO: Why? Filter coefs summ up to 128 applied 2x 128*128= 2^14
      //TODO: Why?  Filter coefs summ up to 64 applied 2x 64*64= 2^12
      //Scale values back down
      tmp_col[j] = is_upscaling ? (tmp_col[j] + 2048) >> 12 : (tmp_col[j] + 8192) >> 14;
    }

    //Clip and move to buffer data
    for (int n = 0; n < trgt_height; n++) {
      src_col[n * buffer->width] = clip(tmp_col[n], 0, 255);
    }
  }
}

//Do the resampling in one pass using 2D-convolution. TODO: Allow doing the resampling on the specified sub-blocks of the target (the input should be a pointer to the full image buffer and a target buffer also pointing to the full target image).
static void resampleBlock( const pic_buffer_t* const src_buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const pic_buffer_t *const trgt_buffer, const int block_x, const int block_y, const int block_width, const int block_height )
{
  //TODO: Add cropping etc.

  //Choose best filter to use when downsampling
  //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
  int ver_filter = 0;
  int hor_filter = 0;

  int src_width = param->src_width + param->src_padding_x;
  int src_height = param->src_height + param->src_padding_y;
  int trgt_width = param->rnd_trgt_width;
  int trgt_height = param->rnd_trgt_height;

  if (!is_upscaling) {
    int crop_width = src_width - param->right_offset; //- param->left_offset;
    int crop_height = src_height - param->bottom_offset; //- param->top_offset;

    if (4 * crop_height > 15 * trgt_height)
      ver_filter = 7;
    else if (7 * crop_height > 20 * trgt_height)
      ver_filter = 6;
    else if (2 * crop_height > 5 * trgt_height)
      ver_filter = 5;
    else if (1 * crop_height > 2 * trgt_height)
      ver_filter = 4;
    else if (3 * crop_height > 5 * trgt_height)
      ver_filter = 3;
    else if (4 * crop_height > 5 * trgt_height)
      ver_filter = 2;
    else if (19 * crop_height > 20 * trgt_height)
      ver_filter = 1;

    if (4 * crop_width > 15 * trgt_width)
      hor_filter = 7;
    else if (7 * crop_width > 20 * trgt_width)
      hor_filter = 6;
    else if (2 * crop_width > 5 * trgt_width)
      hor_filter = 5;
    else if (1 * crop_width > 2 * trgt_width)
      hor_filter = 4;
    else if (3 * crop_width > 5 * trgt_width)
      hor_filter = 3;
    else if (4 * crop_width > 5 * trgt_width)
      hor_filter = 2;
    else if (19 * crop_width > 20 * trgt_width)
      hor_filter = 1;
  }

  int shift_x = param->shift_x - 4;
  int shift_y = param->shift_y - 4;

  //Get the pointer to the target and source data.
  pic_data_t *src = src_buffer->data;
  pic_data_t *trgt = trgt_buffer->data; //&trgt_buffer->data[block_x + block_y*trgt_buffer->width];
  
  //Loop over every pixel in the target block and calculate the 2D-convolution to get the resampled value for the given pixel
  for( int y = block_y; y < (block_y + block_height); y++ ){
    for( int x = block_x; x < (block_x + block_width); x++ ){
      
      //Calculate reference position in src pic
      int ref_pos_x_16 = (int)((unsigned int)(x * param->scale_x + param->add_x) >> shift_x) - param->delta_x;
      int ref_pos_y_16 = (int)((unsigned int)(y * param->scale_y + param->add_y) >> shift_y) - param->delta_y;

      int phase_x = ref_pos_x_16 & 15;
      int phase_y = ref_pos_y_16 & 15;      

      int ref_pos_x = ref_pos_x_16 >> 4;
      int ref_pos_y = ref_pos_y_16 >> 4;

      //Choose filter
      const int *filter_x;
      const int *filter_y;
      const int size_x = getFilter(&filter_x, is_upscaling, is_luma, phase_x, hor_filter);
      const int size_y = getFilter(&filter_y, is_upscaling, is_luma, phase_y, ver_filter);

      pic_data_t new_val = 0; //Accumulate the new pixel value here

      //Convolution kernel, where (x,y)<=>(0,0)
      //Size of kernel depends on the filter size
      for( int j = 0; j < size_y; j++ ){
        //Calculate src sample position for kernel element (i,j)
        int m_y = clip(ref_pos_y + j - (size_y >> 1) + 1, 0, src_height - 1);

        for (int i = 0; i < size_x; i++) {
          //Calculate src sample position for kernel element (i,j)
          int m_x = clip( ref_pos_x + i - (size_x >> 1) + 1, 0, src_width - 1);

          //Calculate filter value in the 2D-filter for pos (i,j) and apply to sample (m_x,m_y)
          new_val += filter_x[i]*filter_y[j] * src[m_x + m_y*src_buffer->width];
        }
      }

      //Scale and clip values and save to trgt buffer.
      trgt[x + y*trgt_buffer->width] = clip(is_upscaling ? (new_val + 2048) >> 12 : (new_val + 8192) >> 14, 0, 255); //TODO: account for different bit dept / different filters etc
    }  
  }
}

//Calculate scaling parameters and update param. Factor determines if certain values are 
// divided eg. with chroma. 0 for no factor and -1 for halving stuff and 1 for doubling etc.
//Calculations from SHM
static void calculateParameters(scaling_parameter_t* const param, const int w_factor, const int h_factor, const int is_chroma)
{
  //First shift widths/height by an appropriate factor
  param->src_width = SHIFT(param->src_width, w_factor);
  param->src_height = SHIFT(param->src_height, h_factor);
  param->trgt_width = SHIFT(param->trgt_width, w_factor);
  param->trgt_height = SHIFT(param->trgt_height, h_factor);
  param->scaled_src_width = SHIFT(param->scaled_src_width, w_factor);
  param->scaled_src_height = SHIFT(param->scaled_src_height, h_factor);
  param->rnd_trgt_width = SHIFT(param->rnd_trgt_width, w_factor);
  param->rnd_trgt_height = SHIFT(param->rnd_trgt_height, h_factor);

  //Calculate sample positional parameters
  param->right_offset = param->src_width - param->scaled_src_width; //- left_offset
  param->bottom_offset = param->src_height - param->scaled_src_height; //- top_offset

  //TODO: Make dependant on width/height?
  param->shift_x = 16;
  param->shift_y = 16;

  param->scale_x = (((unsigned int)param->scaled_src_width << param->shift_x) + (param->rnd_trgt_width >> 1)) / param->rnd_trgt_width;
  param->scale_y = (((unsigned int)param->scaled_src_height << param->shift_y) + (param->rnd_trgt_height >> 1)) / param->rnd_trgt_height;

  //Phase calculations
  //param->phase_x = 0;
  //param->phase_y = 0;
  int phase_x = 0;
  int phase_y = 0;
  //Hardcode phases for chroma, values from SHM. TODO: Find out why these values?
  if( is_chroma != 0 && param->chroma!=CHROMA_444 ) {
    //param->phase_y = 1;
    phase_y = 1;
  }

  //TODO: Is delta_? strictly necessary?
  param->add_x = (((param->scaled_src_width * phase_x) << (param->shift_x - 2)) + (param->rnd_trgt_width >> 1)) / param->rnd_trgt_width + (1 << (param->shift_x - 5));
  param->add_y = (((param->scaled_src_height * phase_y) << (param->shift_y - 2)) + (param->rnd_trgt_height >> 1)) / param->rnd_trgt_height + (1 << (param->shift_y - 5));
  //param->add_x = -(((phase_x * param->scale_x + 8) >> 4 ) - (1 << (param->shift_x - 5)));
  //param->add_y = -(((phase_y * param->scale_y + 8) >> 4 ) - (1 << (param->shift_y - 5)));

  param->delta_x = 4 * phase_x; //- (left_offset << 4)
  param->delta_y = 4 * phase_y; //- (top_offset << 4)
}

scaling_parameter_t kvz_newScalingParameters(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma)
{
  scaling_parameter_t param = {
    .src_width = src_width,
    .src_height = src_height,
    .trgt_width = trgt_width,
    .trgt_height = trgt_height,
    .chroma = chroma,
    .src_padding_x = 0,
    .src_padding_y = 0,
    .trgt_padding_x = 0,
    .trgt_padding_y = 0
  };

  //Calculate Resampling parameters
  //Calculations from SHM
  int hor_div = param.trgt_width << 1;
  int ver_div = param.trgt_height << 1;

  param.rnd_trgt_width = ((param.trgt_width + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16
  param.rnd_trgt_height = ((param.trgt_height + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16

  //Round to multiple of 2
  //TODO: Why MAX? Try using src
  int scaled_src_width = param.src_width;//MAX(param.src_width, param.trgt_width); //Min/max of source/target values
  int scaled_src_height = param.src_height;//MAX(param.src_height, param.trgt_height); //Min/max of source/target values
  param.scaled_src_width = ((scaled_src_width * param.rnd_trgt_width + (hor_div >> 1)) / hor_div) << 1;
  param.scaled_src_height = ((scaled_src_height * param.rnd_trgt_height + (ver_div >> 1)) / ver_div) << 1;

  //Pre-Calculate other parameters
  calculateParameters(&param, 0, 0, 0);

  return param;
}

scaling_parameter_t kvz_newScalingParameters_(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma)
{
  scaling_parameter_t param = {
    .src_width = src_width,
    .src_height = src_height,
    .trgt_width = trgt_width,
    .trgt_height = trgt_height,
    .chroma = chroma
  };

  //Calculate Resampling parameters
  //Calculations from SHM
  int hor_div = param.trgt_width << 1;
  int ver_div = param.trgt_height << 1;

  param.rnd_trgt_width = param.trgt_width;//((param.trgt_width + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16
  param.rnd_trgt_height = param.trgt_height;//((param.trgt_height + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16

  //Round to multiple of 2
  //TODO: Why MAX? Try using src
  int scaled_src_width = param.src_width;//MAX(param.src_width, param.trgt_width); //Min/max of source/target values
  int scaled_src_height = param.src_height;//MAX(param.src_height, param.trgt_height); //Min/max of source/target values
  param.scaled_src_width = ((scaled_src_width * param.rnd_trgt_width + (hor_div >> 1)) / hor_div) << 1;
  param.scaled_src_height = ((scaled_src_height * param.rnd_trgt_height + (ver_div >> 1)) / ver_div) << 1;

  //Pre-Calculate other parameters
  calculateParameters(&param, 0, 0, 0);

  return param;
}


chroma_format_t kvz_getChromaFormat(int luma_width, int luma_height, int chroma_width, int chroma_height)
{
  if (chroma_width == 0 && chroma_height == 0) {
    return CHROMA_400;
  }
  if (chroma_width == luma_width && chroma_height == luma_height) {
    return CHROMA_444;
  }
  if (chroma_width == (luma_width >> 1) && chroma_height == (luma_height)) {
    return CHROMA_422;
  }
  //If not CHROMA_420, not a supported format
  assert(chroma_width == (luma_width >> 1) && chroma_height == (luma_height >> 1));

  return CHROMA_420;
}


yuv_buffer_t* kvz_yuvScaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param,
                         yuv_buffer_t* dst)
{
  /*========== Basic Initialization ==============*/
  //Initialize basic parameters
  scaling_parameter_t param = *base_param;

  //How much to scale the luma sizes to get the chroma sizes
  int w_factor = 0;
  int h_factor = 0;
  switch (param.chroma) {
    case CHROMA_400: {
      //No chroma
      assert(yuv->u->height == 0 && yuv->u->width == 0 && yuv->v->height == 0 && yuv->v->width == 0);
      break;
    }
    case CHROMA_420: {
      assert(yuv->u->height == (yuv->y->height >> 1) && yuv->u->width == (yuv->y->width >> 1)
        && yuv->v->height == (yuv->y->height >> 1) && yuv->v->width == (yuv->y->width >> 1));
      w_factor = -1;
      h_factor = -1;
      break;
    }
    case CHROMA_422: {
      assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width >> 1)
        && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width >> 1));
      w_factor = -1;
      break;
    }
    case CHROMA_444: {
      assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width)
        && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width));
      break;
    }
    default:
      assert(0); //Unsupported chroma type
  }

  //Check if base param and yuv buffer are the same size, if yes we can asume parameters are initialized
  if (yuv->y->width != param.src_width + param.src_padding_x || yuv->y->height != param.src_height + param.src_padding_y) {
    param.src_width = yuv->y->width - param.src_padding_x;
    param.src_height = yuv->y->height - param.src_padding_y;
    calculateParameters(&param, 0, 0, 0);
  }

  //Check if we need to allocate a yuv buffer for the new image or re-use dst.
  //Make sure the sizes match
  if (dst == NULL || dst->y->width != param.trgt_width || dst->y->height != param.trgt_height
    || dst->u->width != SHIFT(param.trgt_width, w_factor) || dst->u->height != SHIFT(param.trgt_height, h_factor)
    || dst->v->width != SHIFT(param.trgt_width, w_factor) || dst->v->height != SHIFT(param.trgt_height, h_factor)) {

    kvz_deallocateYuvBuffer(dst); //Free old buffer if it exists

    dst = kvz_newYuvBuffer(param.trgt_width, param.trgt_height, param.chroma, 0);
  }

  //Calculate if we are upscaling or downscaling
  int is_downscaled_width = base_param->src_width > base_param->trgt_width;
  int is_downscaled_height = base_param->src_height > base_param->trgt_height;
  int is_equal_width = base_param->src_width == base_param->trgt_width;
  int is_equal_height = base_param->src_height == base_param->trgt_height;

  int is_upscaling = 1;

  //both dimensions need to be either up- or downscaled
  if ((is_downscaled_width && !is_downscaled_height && !is_equal_height) ||
      (is_downscaled_height && !is_downscaled_width && !is_equal_width)) {
    return NULL;
  }
  if (is_equal_height && is_equal_width) {
    //If equal just return source
    copyYuvBuffer(yuv, dst, 0);
    return dst;
  }
  if (is_downscaled_width || is_downscaled_height) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }
  /*=================================*/

  //Allocate a pic_buffer to hold the component data while the downscaling is done
  //Size calculation from SHM. TODO: Figure out why. Use yuv as buffer instead?
  int max_width = MAX(param.src_width+param.src_padding_x, param.trgt_width);
  int max_height = MAX(param.src_height+param.src_padding_y, param.trgt_height);
  int min_width = MIN(param.src_width, param.trgt_width);
  int min_height = MIN(param.src_height, param.trgt_height);
  int min_width_rnd16 = ((min_width + 15) >> 4) << 4;
  int min_height_rnd32 = ((min_height + 31) >> 5) << 5;
  int buffer_width = ((max_width * min_width_rnd16 + (min_width << 4) - 1) / (min_width << 4)) << 4;
  int buffer_height = ((max_height * min_height_rnd32 + (min_height << 4) - 1) / (min_height << 4)) << 4;;
  pic_buffer_t* buffer = kvz_newPictureBuffer(buffer_width, buffer_height, 1);


  /*==========Start Resampling=============*/
  //Resample y
  copyPictureBuffer(yuv->y, buffer, 1);
  resample(buffer, &param, is_upscaling, 1);
  copyPictureBuffer(buffer, dst->y, 0);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //Resample u
    copyPictureBuffer(yuv->u, buffer, 1);
    resample(buffer, &param, is_upscaling, 0);
    copyPictureBuffer(buffer, dst->u, 0);

    //Resample v
    copyPictureBuffer(yuv->v, buffer, 1);
    resample(buffer, &param, is_upscaling, 0);
    copyPictureBuffer(buffer, dst->v, 0);
  }

  //Deallocate buffer
  deallocatePictureBuffer(buffer);

  return dst;
}

//Use yuv and dst as the buffer instead of allocating a new buffer. Also use unrounded sizes
//yuv is not quaranteet to contain the original data.
yuv_buffer_t* kvz_yuvScaling_(yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param,
                          yuv_buffer_t* dst)
{
  /*========== Basic Initialization ==============*/
  //Initialize basic parameters
  scaling_parameter_t param = *base_param;

  //How much to scale the luma sizes to get the chroma sizes
  int w_factor = 0;
  int h_factor = 0;
  switch (param.chroma) {
    case CHROMA_400: {
      //No chroma
      assert(yuv->u->height == 0 && yuv->u->width == 0 && yuv->v->height == 0 && yuv->v->width == 0);
      break;
    }
    case CHROMA_420: {
      assert(yuv->u->height == (yuv->y->height >> 1) && yuv->u->width == (yuv->y->width >> 1)
        && yuv->v->height == (yuv->y->height >> 1) && yuv->v->width == (yuv->y->width >> 1));
      w_factor = -1;
      h_factor = -1;
      break;
    }
    case CHROMA_422: {
      assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width >> 1)
        && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width >> 1));
      w_factor = -1;
      break;
    }
    case CHROMA_444: {
      assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width)
        && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width));
      break;
    }
    default:
      assert(0); //Unsupported chroma type
  }

  //Check if base param and yuv buffer are the same size, if yes we can asume parameters are initialized
  if (yuv->y->width != param.src_width || yuv->y->height != param.src_height) {
    param.src_width = yuv->y->width;
    param.src_height = yuv->y->height;
    calculateParameters(&param, w_factor, h_factor, 0);
  }

  //Check if we need to allocate a yuv buffer for the new image or re-use dst.
  //Make sure the sizes match
  if (dst == NULL || dst->y->width != param.trgt_width || dst->y->height != param.trgt_height
    || dst->u->width != SHIFT(param.trgt_width, w_factor) || dst->u->height != SHIFT(param.trgt_height, h_factor)
    || dst->v->width != SHIFT(param.trgt_width, w_factor) || dst->v->height != SHIFT(param.trgt_height, h_factor)) {

    kvz_deallocateYuvBuffer(dst); //Free old buffer if it exists

    dst = kvz_newYuvBuffer(param.trgt_width, param.trgt_height, param.chroma, 0);
  }

  //Calculate if we are upscaling or downscaling
  int is_downscaled_width = base_param->src_width > base_param->trgt_width;
  int is_downscaled_height = base_param->src_height > base_param->trgt_height;
  int is_equal_width = base_param->src_width == base_param->trgt_width;
  int is_equal_height = base_param->src_height == base_param->trgt_height;

  int is_upscaling = 1;

  //both dimensions need to be either up- or downscaled
  if ((is_downscaled_width && !is_downscaled_height && !is_equal_height) ||
      (is_downscaled_height && !is_downscaled_width && !is_equal_width)) {
    return NULL;
  }
  if (is_equal_height && is_equal_width) {
    //If equal just return source
    copyYuvBuffer(yuv, dst, 0);
    return dst;
  }
  if (is_downscaled_width || is_downscaled_height) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }
  /*=================================*/

  //Allocate a pic_buffer to hold the component data while the downscaling is done
  //Size calculation from SHM. TODO: Figure out why. Use yuv as buffer instead?
  /*int max_width = MAX(param.src_width, param.trgt_width);
  int max_height = MAX(param.src_height, param.trgt_height);
  int min_width = MIN(param.src_width, param.trgt_width);
  int min_height = MIN(param.src_height, param.trgt_height);
  int min_width_rnd16 = ((min_width + 15) >> 4) << 4;
  int min_height_rnd32 = ((min_height + 31) >> 5) << 5;
  int buffer_width = ((max_width * min_width_rnd16 + (min_width << 4) - 1) / (min_width << 4)) << 4;
  int buffer_height = ((max_height * min_height_rnd32 + (min_height << 4) - 1) / (min_height << 4)) << 4;;
  pic_buffer_t* buffer = kvz_newPictureBuffer(buffer_width, buffer_height, 1);*/
  //TODO: Clean up this part and implement properly
  //param.rnd_trgt_height = param.trgt_height;
  //param.rnd_trgt_width = param.trgt_width;
  yuv_buffer_t* buffer = is_upscaling ? dst : yuv;//malloc(sizeof(pic_buffer_t)); //Choose bigger buffer
  if (buffer->y->tmp_row == NULL) {
    buffer->y->tmp_row = malloc(sizeof(pic_data_t) * (MAX(buffer->y->width, buffer->y->height)));
  }
  if (buffer->u->tmp_row == NULL) {
    buffer->u->tmp_row = malloc(sizeof(pic_data_t) * (MAX(buffer->u->width, buffer->u->height)));
  }
  if (buffer->v->tmp_row == NULL) {
    buffer->v->tmp_row = malloc(sizeof(pic_data_t) * (MAX(buffer->v->width, buffer->v->height)));
  }

  /*==========Start Resampling=============*/
  //Resample y
  if (is_upscaling) copyPictureBuffer(yuv->y, buffer->y, 1);
  _resample(buffer->y, &param, is_upscaling, 1);
  if (!is_upscaling) copyPictureBuffer(buffer->y, dst->y, 0);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //Resample u
    if (is_upscaling) copyPictureBuffer(yuv->u, buffer->u, 1);
    _resample(buffer->u, &param, is_upscaling, 0);
    if (!is_upscaling) copyPictureBuffer(buffer->u, dst->u, 0);

    //Resample v
    if (is_upscaling) copyPictureBuffer(yuv->v, buffer->v, 1);
    _resample(buffer->v, &param, is_upscaling, 0);
    if (!is_upscaling) copyPictureBuffer(buffer->v, dst->v, 0);
  }

  //Deallocate buffer
  //deallocatePictureBuffer(buffer);

  return dst;
}

// yuv buffer should not be modified
int kvz_yuvBlockScaling( const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst, const int block_x, const int block_y, const int block_width, const int block_height )
{
  /*========== Basic Initialization ==============*/

  //Check that block parameters are valid
  if( block_x < 0 || block_y < 0 || block_x + block_width > base_param->trgt_width || block_y + block_height > base_param->trgt_height ){
    fprintf(stderr, "Specified block outside given target picture size.");
    return 0;
  }

  //Initialize basic parameters
  scaling_parameter_t param = *base_param;

  //How much to scale the luma sizes to get the chroma sizes
  int w_factor = 0;
  int h_factor = 0;
  switch (param.chroma) {
  case CHROMA_400: {
    //No chroma
    assert(yuv->u->height == 0 && yuv->u->width == 0 && yuv->v->height == 0 && yuv->v->width == 0);
    break;
  }
  case CHROMA_420: {
    assert(yuv->u->height == (yuv->y->height >> 1) && yuv->u->width == (yuv->y->width >> 1)
      && yuv->v->height == (yuv->y->height >> 1) && yuv->v->width == (yuv->y->width >> 1));
    w_factor = -1;
    h_factor = -1;
    break;
  }
  case CHROMA_422: {
    assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width >> 1)
      && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width >> 1));
    w_factor = -1;
    break;
  }
  case CHROMA_444: {
    assert(yuv->u->height == (yuv->y->height) && yuv->u->width == (yuv->y->width)
      && yuv->v->height == (yuv->y->height) && yuv->v->width == (yuv->y->width));
    break;
  }
  default:
    assert(0); //Unsupported chroma type
  }

  //Check if the buffers are large enough for the given parameters and destination is set.
  if (yuv == NULL || yuv->y->width < param.src_width + param.src_padding_x || yuv->y->height < param.src_height + param.src_padding_y || yuv->u->width < SHIFT(param.src_width + param.src_padding_x, w_factor) || yuv->u->height < SHIFT(param.src_height + param.src_padding_y, w_factor) || yuv->v->width < SHIFT(param.src_width + param.src_padding_x, w_factor) || yuv->v->height < SHIFT(param.src_height + param.src_padding_y, w_factor)) {
    fprintf(stderr, "Source buffer smaller than specified in the scaling parameters.\n");
    return 0;
  }
 
  if (dst == NULL || dst->y->width < param.trgt_width || dst->y->height < param.trgt_height
    || dst->u->width < SHIFT(param.trgt_width, w_factor) || dst->u->height < SHIFT(param.trgt_height, h_factor)
    || dst->v->width < SHIFT(param.trgt_width, w_factor) || dst->v->height < SHIFT(param.trgt_height, h_factor)) {

    fprintf(stderr, "Destination buffer smaller than specified in the scaling parameters,\n");
    return 0;
  }

  //Calculate if we are upscaling or downscaling
  int is_downscaled_width = param.src_width > param.trgt_width;
  int is_downscaled_height = param.src_height > param.trgt_height;
  int is_equal_width = param.src_width == param.trgt_width;
  int is_equal_height = param.src_height == param.trgt_height;

  int is_upscaling = 1;

  //both dimensions need to be either up- or downscaled
  if ((is_downscaled_width && !is_downscaled_height && !is_equal_height) ||
    (is_downscaled_height && !is_downscaled_width && !is_equal_width)) {
    fprintf(stderr, "Both demensions need to be either upscaled or downscaled");
    return 0;
  }
  if (is_equal_height && is_equal_width) {
    //If equal just copy block from src
    copyYuvBufferBlock(yuv, dst, block_x, block_y, block_width, block_height, w_factor, h_factor);
    return 1;
  }
  if (is_downscaled_width || is_downscaled_height) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }
  /*=================================*/

  /*==========Start Resampling=============*/
  //Resample y
  resampleBlock(yuv->y, &param, is_upscaling, 1, dst->y, block_x, block_y, block_width, block_height);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //Resample u
    resampleBlock(yuv->u, &param, is_upscaling, 0, dst->u, SHIFT(block_x, w_factor), SHIFT(block_y, h_factor), SHIFT(block_width, w_factor), SHIFT(block_height, h_factor) );

    //Resample v
    resampleBlock(yuv->v, &param, is_upscaling, 0, dst->v, SHIFT(block_x, w_factor), SHIFT(block_y, h_factor), SHIFT(block_width, w_factor), SHIFT(block_height, h_factor) );
  }

  return 1;
}

static void blockScalingSrcRange( int range[2], const int scale, const int add, const int shift, const int delta, const int block_low, const int block_high, int is_upsampling )
{
  //Get filter size
  int size = is_upsampling ? sizeof(lumaUpFilter[0]) / sizeof(lumaUpFilter[0][0]) : sizeof(filter16[0][0]) / sizeof(filter16[0][0][0]);

  //Calculate lower bound
  range[0] = ((int)((unsigned int)((block_low * scale + add) >> (shift - 4)) - delta) >> 4) - (size >> 1) + 1;

  //Calculate upper bound
  range[1] = ((int)((unsigned int)((block_high * scale + add) >> (shift - 4)) - delta) >> 4) - (size >> 1) + size;
}

void kvz_blockScalingSrcWidthRange(int range[2], const scaling_parameter_t * const base_param, const int block_x, const int block_width, int is_upsampling)
{
  //Calculate parameters
  calculateParameters(base_param, 0, 0, 0);

  blockScalingSrcRange(range, base_param->scale_x, base_param->add_x, base_param->shift_x, base_param->delta_x, block_x, block_x + block_width - 1, is_upsampling);
}

void kvz_blockScalingSrcHeightRange(int range[2], const scaling_parameter_t * const base_param, const int block_y, const int block_height, int is_upsampling)
{
  //Calculate parameters
  calculateParameters(base_param, 0, 0, 0);

  blockScalingSrcRange(range, base_param->scale_y, base_param->add_y, base_param->shift_y, base_param->delta_y, block_y, block_y + block_height - 1, is_upsampling);
}
