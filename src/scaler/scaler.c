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
#include <memory.h>

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

#define SHIFT(x,y) (((y) < 0) ? ((x)>>(-(y))) : ((x)<<(y)))

//Define filters for scaling operations
//Values from SHM
static const int filter16[8][16][12] = {
  // ratio <= 20/19
    {
      { 0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 2, -6, 127, 7, -2, 0, 0, 0, 0 },
      { 0, 0, 0, 3, -12, 125, 16, -5, 1, 0, 0, 0 },
      { 0, 0, 0, 4, -16, 120, 26, -7, 1, 0, 0, 0 },
      { 0, 0, 0, 5, -18, 114, 36, -10, 1, 0, 0, 0 },
      { 0, 0, 0, 5, -20, 107, 46, -12, 2, 0, 0, 0 },
      { 0, 0, 0, 5, -21, 99, 57, -15, 3, 0, 0, 0 },
      { 0, 0, 0, 5, -20, 89, 68, -18, 4, 0, 0, 0 },
      { 0, 0, 0, 4, -19, 79, 79, -19, 4, 0, 0, 0 },
      { 0, 0, 0, 4, -18, 68, 89, -20, 5, 0, 0, 0 },
      { 0, 0, 0, 3, -15, 57, 99, -21, 5, 0, 0, 0 },
      { 0, 0, 0, 2, -12, 46, 107, -20, 5, 0, 0, 0 },
      { 0, 0, 0, 1, -10, 36, 114, -18, 5, 0, 0, 0 },
      { 0, 0, 0, 1, -7, 26, 120, -16, 4, 0, 0, 0 },
      { 0, 0, 0, 1, -5, 16, 125, -12, 3, 0, 0, 0 },
      { 0, 0, 0, 0, -2, 7, 127, -6, 2, 0, 0, 0 }
    },
    // 20/19 < ratio <= 5/4
    {
      { 0, 2, 0, -14, 33, 86, 33, -14, 0, 2, 0, 0 },
      { 0, 1, 1, -14, 29, 85, 38, -13, -1, 2, 0, 0 },
      { 0, 1, 2, -14, 24, 84, 43, -12, -2, 2, 0, 0 },
      { 0, 1, 2, -13, 19, 83, 48, -11, -3, 2, 0, 0 },
      { 0, 0, 3, -13, 15, 81, 53, -10, -4, 3, 0, 0 },
      { 0, 0, 3, -12, 11, 79, 57, -8, -5, 3, 0, 0 },
      { 0, 0, 3, -11, 7, 76, 62, -5, -7, 3, 0, 0 },
      { 0, 0, 3, -10, 3, 73, 65, -2, -7, 3, 0, 0 },
      { 0, 0, 3, -9, 0, 70, 70, 0, -9, 3, 0, 0 },
      { 0, 0, 3, -7, -2, 65, 73, 3, -10, 3, 0, 0 },
      { 0, 0, 3, -7, -5, 62, 76, 7, -11, 3, 0, 0 },
      { 0, 0, 3, -5, -8, 57, 79, 11, -12, 3, 0, 0 },
      { 0, 0, 3, -4, -10, 53, 81, 15, -13, 3, 0, 0 },
      { 0, 0, 2, -3, -11, 48, 83, 19, -13, 2, 1, 0 },
      { 0, 0, 2, -2, -12, 43, 84, 24, -14, 2, 1, 0 },
      { 0, 0, 2, -1, -13, 38, 85, 29, -14, 1, 1, 0 }
    },
    // 5/4 < ratio <= 5/3
    {
      { 0, 5, -6, -10, 37, 76, 37, -10, -6, 5, 0, 0 },
      { 0, 5, -4, -11, 33, 76, 40, -9, -7, 5, 0, 0 },
      { -1, 5, -3, -12, 29, 75, 45, -7, -8, 5, 0, 0 },
      { -1, 4, -2, -13, 25, 75, 48, -5, -9, 5, 1, 0 },
      { -1, 4, -1, -13, 22, 73, 52, -3, -10, 4, 1, 0 },
      { -1, 4, 0, -13, 18, 72, 55, -1, -11, 4, 2, -1 },
      { -1, 4, 1, -13, 14, 70, 59, 2, -12, 3, 2, -1 },
      { -1, 3, 1, -13, 11, 68, 62, 5, -12, 3, 2, -1 },
      { -1, 3, 2, -13, 8, 65, 65, 8, -13, 2, 3, -1 },
      { -1, 2, 3, -12, 5, 62, 68, 11, -13, 1, 3, -1 },
      { -1, 2, 3, -12, 2, 59, 70, 14, -13, 1, 4, -1 },
      { -1, 2, 4, -11, -1, 55, 72, 18, -13, 0, 4, -1 },
      { 0, 1, 4, -10, -3, 52, 73, 22, -13, -1, 4, -1 },
      { 0, 1, 5, -9, -5, 48, 75, 25, -13, -2, 4, -1 },
      { 0, 0, 5, -8, -7, 45, 75, 29, -12, -3, 5, -1 },
      { 0, 0, 5, -7, -9, 40, 76, 33, -11, -4, 5, 0 },
    },
    // 5/3 < ratio <= 2
    {
      { 2, -3, -9, 6, 39, 58, 39, 6, -9, -3, 2, 0 },
      { 2, -3, -9, 4, 38, 58, 43, 7, -9, -4, 1, 0 },
      { 2, -2, -9, 2, 35, 58, 44, 9, -8, -4, 1, 0 },
      { 1, -2, -9, 1, 34, 58, 46, 11, -8, -5, 1, 0 },
      { 1, -1, -8, -1, 31, 57, 47, 13, -7, -5, 1, 0 },
      { 1, -1, -8, -2, 29, 56, 49, 15, -7, -6, 1, 1 },
      { 1, 0, -8, -3, 26, 55, 51, 17, -7, -6, 1, 1 },
      { 1, 0, -7, -4, 24, 54, 52, 19, -6, -7, 1, 1 },
      { 1, 0, -7, -5, 22, 53, 53, 22, -5, -7, 0, 1 },
      { 1, 1, -7, -6, 19, 52, 54, 24, -4, -7, 0, 1 },
      { 1, 1, -6, -7, 17, 51, 55, 26, -3, -8, 0, 1 },
      { 1, 1, -6, -7, 15, 49, 56, 29, -2, -8, -1, 1 },
      { 0, 1, -5, -7, 13, 47, 57, 31, -1, -8, -1, 1 },
      { 0, 1, -5, -8, 11, 46, 58, 34, 1, -9, -2, 1 },
      { 0, 1, -4, -8, 9, 44, 58, 35, 2, -9, -2, 2 },
      { 0, 1, -4, -9, 7, 43, 58, 38, 4, -9, -3, 2 },
    },
    // 2 < ratio <= 5/2
    {
      { -2, -7, 0, 17, 35, 43, 35, 17, 0, -7, -5, 2 },
      { -2, -7, -1, 16, 34, 43, 36, 18, 1, -7, -5, 2 },
      { -1, -7, -1, 14, 33, 43, 36, 19, 1, -6, -5, 2 },
      { -1, -7, -2, 13, 32, 42, 37, 20, 3, -6, -5, 2 },
      { 0, -7, -3, 12, 31, 42, 38, 21, 3, -6, -5, 2 },
      { 0, -7, -3, 11, 30, 42, 39, 23, 4, -6, -6, 1 },
      { 0, -7, -4, 10, 29, 42, 40, 24, 5, -6, -6, 1 },
      { 1, -7, -4, 9, 27, 41, 40, 25, 6, -5, -6, 1 },
      { 1, -6, -5, 7, 26, 41, 41, 26, 7, -5, -6, 1 },
      { 1, -6, -5, 6, 25, 40, 41, 27, 9, -4, -7, 1 },
      { 1, -6, -6, 5, 24, 40, 42, 29, 10, -4, -7, 0 },
      { 1, -6, -6, 4, 23, 39, 42, 30, 11, -3, -7, 0 },
      { 2, -5, -6, 3, 21, 38, 42, 31, 12, -3, -7, 0 },
      { 2, -5, -6, 3, 20, 37, 42, 32, 13, -2, -7, -1 },
      { 2, -5, -6, 1, 19, 36, 43, 33, 14, -1, -7, -1 },
      { 2, -5, -7, 1, 18, 36, 43, 34, 16, -1, -7, -2 }
    },
    // 5/2 < ratio <= 20/7
    {
      { -6, -3, 5, 19, 31, 36, 31, 19, 5, -3, -6, 0 },
      { -6, -4, 4, 18, 31, 37, 32, 20, 6, -3, -6, -1 },
      { -6, -4, 4, 17, 30, 36, 33, 21, 7, -3, -6, -1 },
      { -5, -5, 3, 16, 30, 36, 33, 22, 8, -2, -6, -2 },
      { -5, -5, 2, 15, 29, 36, 34, 23, 9, -2, -6, -2 },
      { -5, -5, 2, 15, 28, 36, 34, 24, 10, -2, -6, -3 },
      { -4, -5, 1, 14, 27, 36, 35, 24, 10, -1, -6, -3 },
      { -4, -5, 0, 13, 26, 35, 35, 25, 11, 0, -5, -3 },
      { -4, -6, 0, 12, 26, 36, 36, 26, 12, 0, -6, -4 },
      { -3, -5, 0, 11, 25, 35, 35, 26, 13, 0, -5, -4 },
      { -3, -6, -1, 10, 24, 35, 36, 27, 14, 1, -5, -4 },
      { -3, -6, -2, 10, 24, 34, 36, 28, 15, 2, -5, -5 },
      { -2, -6, -2, 9, 23, 34, 36, 29, 15, 2, -5, -5 },
      { -2, -6, -2, 8, 22, 33, 36, 30, 16, 3, -5, -5 },
      { -1, -6, -3, 7, 21, 33, 36, 30, 17, 4, -4, -6 },
      { -1, -6, -3, 6, 20, 32, 37, 31, 18, 4, -4, -6 }
    },
    // 20/7 < ratio <= 15/4
    {
      { -9, 0, 9, 20, 28, 32, 28, 20, 9, 0, -9, 0 },
      { -9, 0, 8, 19, 28, 32, 29, 20, 10, 0, -4, -5 },
      { -9, -1, 8, 18, 28, 32, 29, 21, 10, 1, -4, -5 },
      { -9, -1, 7, 18, 27, 32, 30, 22, 11, 1, -4, -6 },
      { -8, -2, 6, 17, 27, 32, 30, 22, 12, 2, -4, -6 },
      { -8, -2, 6, 16, 26, 32, 31, 23, 12, 2, -4, -6 },
      { -8, -2, 5, 16, 26, 31, 31, 23, 13, 3, -3, -7 },
      { -8, -3, 5, 15, 25, 31, 31, 24, 14, 4, -3, -7 },
      { -7, -3, 4, 14, 25, 31, 31, 25, 14, 4, -3, -7 },
      { -7, -3, 4, 14, 24, 31, 31, 25, 15, 5, -3, -8 },
      { -7, -3, 3, 13, 23, 31, 31, 26, 16, 5, -2, -8 },
      { -6, -4, 2, 12, 23, 31, 32, 26, 16, 6, -2, -8 },
      { -6, -4, 2, 12, 22, 30, 32, 27, 17, 6, -2, -8 },
      { -6, -4, 1, 11, 22, 30, 32, 27, 18, 7, -1, -9 },
      { -5, -4, 1, 10, 21, 29, 32, 28, 18, 8, -1, -9 },
      { -5, -4, 0, 10, 20, 29, 32, 28, 19, 8, 0, -9 }
    },
    // ratio > 15/4
    {
      { -8, 7, 13, 18, 22, 24, 22, 18, 13, 7, 2, -10 },
      { -8, 7, 13, 18, 22, 23, 22, 19, 13, 7, 2, -10 },
      { -8, 6, 12, 18, 22, 23, 22, 19, 14, 8, 2, -10 },
      { -9, 6, 12, 17, 22, 23, 23, 19, 14, 8, 3, -10 },
      { -9, 6, 12, 17, 21, 23, 23, 19, 14, 9, 3, -10 },
      { -9, 5, 11, 17, 21, 23, 23, 20, 15, 9, 3, -10 },
      { -9, 5, 11, 16, 21, 23, 23, 20, 15, 9, 4, -10 },
      { -9, 5, 10, 16, 21, 23, 23, 20, 15, 10, 4, -10 },
      { -10, 5, 10, 16, 20, 23, 23, 20, 16, 10, 5, -10 },
      { -10, 4, 10, 15, 20, 23, 23, 21, 16, 10, 5, -9 },
      { -10, 4, 9, 15, 20, 23, 23, 21, 16, 11, 5, -9 },
      { -10, 3, 9, 15, 20, 23, 23, 21, 17, 11, 5, -9 },
      { -10, 3, 9, 14, 19, 23, 23, 21, 17, 12, 6, -9 },
      { -10, 3, 8, 14, 19, 23, 23, 22, 17, 12, 6, -9 },
      { -10, 2, 8, 14, 19, 22, 23, 22, 18, 12, 6, -8 },
      { -10, 2, 7, 13, 19, 22, 23, 22, 18, 13, 7, -8 }
    }
};

//Used for clipping values
int clip(int val, int min, int max)
{
  if (val <= min)
    return min;
  if (val >= max)
    return max;

  return val;
}

 pic_buffer_t* newPictureBuffer(int width, int height, int has_tmp_row)
 {
   //TODO: Add error checking
   pic_buffer_t* buffer = (pic_buffer_t*)malloc(sizeof(pic_buffer_t));
   
   //Allocate enough memory to fit a width-by-height picture
   buffer->data = (pic_data_t*)malloc(sizeof(pic_data_t)*width*height);

   buffer->width = width;
   buffer->height = height;

   //Initialize tmp_row or set as NULL
   if (has_tmp_row) {
     buffer->tmp_row = (pic_data_t*)malloc(sizeof(pic_data_t)*width);
   }
   else {
     buffer->tmp_row = NULL;
   }

   return buffer;
 }

 // ======================= newPictureBuffer_ ==================================
 //TODO: DO something about the lack of overloading?
 pic_buffer_t* newPictureBuffer_double(double* data, int width, int height, int has_tmp_row)
 {
   pic_buffer_t* buffer = newPictureBuffer(width, height, has_tmp_row);

   //Initialize buffer
   for (int i = width*height-1; i >= 0; i--) {
     buffer->data[i] = (int)data[i];
   }

   return buffer;
 }

 pic_buffer_t* newPictureBuffer_uint8(uint8_t* data, int width, int height, int has_tmp_row)
 {
   pic_buffer_t* buffer = newPictureBuffer(width, height, has_tmp_row);

   //Initialize buffer
   for (int i = width*height-1; i >= 0; i--) {
     buffer->data[i] = (int)data[i];
   }

   return buffer;
 }
 // ==============================================================================

 void deallocatePictureBuffer(pic_buffer_t* buffer)
 {
   free(buffer->data);
   free(buffer->tmp_row);
   free(buffer);
 }

 void copyPictureBuffer(pic_buffer_t* src, pic_buffer_t* dst, int fill)
 {
   //TODO: add checks. Check if fill is necessary
   //max_dim_* is chosen so that no over indexing happenes (src/dst)
   //min_dim_* is chosen so that no over indexing happenes (src), but all inds in dst get a value
   int max_dim_x = fill ? dst->width : MIN(src->width,dst->width);
   int max_dim_y = fill ? dst->height : MIN(src->height,dst->height);
   int min_dim_x = fill ? src->width : max_dim_x;
   int min_dim_y = fill ? src->height : max_dim_y;

   int dst_row = 0;
   int src_row = 0;

   //Copy loop
   for (int i = 0; i < max_dim_y ; i++) {
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

 // ======================= newYuvBuffer_ ==================================
 yuv_buffer_t* newYuvBuffer_double(double* y_data, double* u_data, double* v_data, int width, int height, int is_420)
 {
   yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));

   //Allocate y pic_buffer
   yuv->y = newPictureBuffer_double(y_data, width, height, 0);
   
   //Allocate u and v buffers
   is_420 = is_420 ? 1 : 0;
   width >>= is_420;
   height >>= is_420;
   yuv->u = newPictureBuffer_double(u_data, width, height, 0);
   yuv->v = newPictureBuffer_double(v_data, width, height, 0);

   return yuv;
 }

 yuv_buffer_t* newYuvBuffer_uint8(uint8_t* y_data, uint8_t* u_data, uint8_t* v_data, int width, int height, int is_420)
 {
   yuv_buffer_t* yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));

   //Allocate y pic_buffer
   yuv->y = newPictureBuffer_uint8(y_data, width, height, 0);

   //Allocate u and v buffers
   is_420 = is_420 ? 1 : 0;
   width >>= is_420;
   height >>= is_420;
   yuv->u = newPictureBuffer_uint8(u_data, width, height, 0);
   yuv->v = newPictureBuffer_uint8(v_data, width, height, 0);

   return yuv;
 }
 // ==============================================================================

 void deallocateYuvBuffer(yuv_buffer_t* yuv)
 {
   deallocatePictureBuffer(yuv->y);
   deallocatePictureBuffer(yuv->u);
   deallocatePictureBuffer(yuv->v);

   free(yuv);
 }

 //Actual downscaling is done here
 void downsample(pic_buffer_t* buffer, scaling_parameter_t* param)
 {
   //TODO: Add cropping etc.
   
   //Choose best filter to use
   //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
   int ver_filter = 0;
   int hor_filter = 0;

   int src_height = param->src_height;
   int src_width = param->src_width;
   int trgt_height = param->rnd_trgt_height;
   int trgt_width = param->rnd_trgt_width;
   
   int crop_width = src_width - param->right_offset; //- param->left_offset;
   int crop_height = src_height - param->bottom_offset; //- param->top_offset;

   if (4 * crop_height > 15 * trgt_height)
     ver_filter = 7;
   else if (7 * crop_height > 20 * trgt_height)
     ver_filter = 6;
   else if ( 2 * crop_height > 5 * trgt_height)
     ver_filter = 5;
   else if ( 1 * crop_height > 2 * trgt_height)
     ver_filter = 4;
   else if ( 3 * crop_height > 5 * trgt_height)
     ver_filter = 3;
   else if ( 4 * crop_height > 5 * trgt_height)
     ver_filter = 2;
   else if ( 19 * crop_height > 20 * trgt_height)
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

   int shift_x = param->shift_x - 4;
   int shift_y = param->shift_y - 4;

   pic_data_t* tmp_row = buffer->tmp_row;

   // Horizontal downsampling
   for (int i = 0; i < src_height; i++) {
     pic_data_t* src_row = &buffer->data[i*buffer->width];

     for (int j = 0; j < trgt_width; j++) {
       //Calculate reference position in src pic
       int ref_pos_16 = (int)((unsigned int)(j*param->scale_x + param->add_x) >> shift_x);
       int phase = ref_pos_16 & 15;
       int ref_pos = ref_pos_16 >> 4;

       //Apply filter
       tmp_row[j] = 0;
       for (int k = 0; k < 12; k++) {
         int m = clip(ref_pos + k - 5, 0, src_width - 1);
         tmp_row[j] += filter16[hor_filter][phase][k] * src_row[m];
       }
     }
     //Copy tmp row to data
     memcpy(src_row, tmp_row, sizeof(pic_data_t)*trgt_width);
   }

   // Vertical downsampling
   for (int i = 0; i < trgt_width; i++) {
     pic_data_t* src_col = &buffer->data[i];
     for (int j = 0; j < trgt_height; j++) {
       //Calculate ref pos
       int ref_pos_16 = (int)((unsigned int)(j*param->scale_y + param->add_y) >> shift_y);
       int phase = ref_pos_16 & 15;
       int ref_pos = ref_pos_16 >> 4;

       //Apply filter
       tmp_row[j] = 0;
       for (int k = 0; k < 12; k++) {
         int m = clip(ref_pos + k - 5, 0, src_height - 1);
         tmp_row[j] += filter16[ver_filter][phase][k] * src_col[m*buffer->width];
       }
       //TODO: Why?
       //Scale values back down
       tmp_row[j] = (tmp_row[j] + 8192) >> 14;
     }

     //Clip and move to buffer data
     for (int n = 0; n < trgt_height; n++) {
       src_col[n*buffer->width] = clip(tmp_row[n], 0, 255);
     }
   }
 }


//Calculate scaling parameters and update param. Factor determines if certain values are 
// divided eg. with chroma. 0 for no factor and -1 for halving stuff and 1 for doubling etc.
//Calculations from SHM
void calculateParameters(scaling_parameter_t* param, const int w_factor, const int h_factor)
{
  //First shift widths/height by an appropriate factor
  param->src_width = SHIFT(param->src_width, w_factor);
  param->src_height = SHIFT(param->src_height, h_factor);
  param->trgt_width = SHIFT(param->trgt_width, w_factor);
  param->trgt_height = SHIFT(param->trgt_height, h_factor);
  param->rnd_src_width = SHIFT(param->rnd_src_width, w_factor);
  param->rnd_src_height = SHIFT(param->rnd_src_height, h_factor);
  param->rnd_trgt_width = SHIFT(param->rnd_trgt_width, w_factor);
  param->rnd_trgt_height = SHIFT(param->rnd_trgt_height, h_factor);

  //Calculate sample positional parameters
  param->right_offset = param->src_width - param->rnd_src_width; //- left_offset
  param->bottom_offset = param->src_height - param->rnd_src_height; //- top_offset

  //TODO: Make dependant on width/heiht?
  param->shift_x = 16;
  param->shift_y = 16;

  param->scale_x = (((unsigned int)param->rnd_src_width << param->shift_x) + (param->rnd_trgt_width >> 1)) / param->rnd_trgt_width;
  param->scale_y = (((unsigned int)param->rnd_src_height << param->shift_y) + (param->rnd_trgt_height >> 1)) / param->rnd_trgt_height;

  //TODO: Add dependace to phase?
  param->add_x = (param->rnd_trgt_width >> 1) / param->rnd_trgt_width + (1 << (param->shift_x - 5));
  param->add_y = (param->rnd_trgt_height >> 1) / param->rnd_trgt_height + (1 << (param->shift_y - 5));

}

scaling_parameter_t newScalingParameters(int src_width, int src_height, int trgt_width, int trgt_height)
{
  scaling_parameter_t param = {
    .src_width = src_width,
    .src_height = src_height,
    .trgt_width = trgt_width,
    .trgt_height = trgt_height
  };

  //Calculate Resampling parameters
  //Calculations from SHM
  int hor_div = param.trgt_width;
  int ver_div = param.trgt_height;

  param.rnd_trgt_height = ((param.trgt_height + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16
  param.rnd_trgt_width = ((param.trgt_width + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16

  //Round to multiple of 2
  param.rnd_src_height = ((param.src_height * param.rnd_trgt_height + (ver_div >> 1)) / ver_div) << 1;
  param.rnd_src_width = ((param.src_width * param.rnd_trgt_width + (hor_div >> 1)) / hor_div) << 1;

  //Pre-Calculate other parameters
  calculateParameters(&param, 0, 0);

  return param;
}

 
 yuv_buffer_t* yuvDownscaling( const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, int is_420)
 {
   //Initialize basic parameters
   scaling_parameter_t param = *base_param;

   //Calculate scaling parameters
   int w_factor = 0;
   int h_factor = 0;
   //calculateParamaeters(&param, w_factor, h_factor);

   //Check if base param and yuv buffer are the same size, if yes we can asume parameters are initialized
   if (yuv->y->width != param.src_width || yuv->y->height != param.src_height) {
     param.src_width = yuv->y->width;
     param.src_height = yuv->y->height;
     calculateParameters(&param, w_factor, h_factor);
   }

   //is_420 = is_420 ? 1 : 0;

   //Allocate a pic_buffer to hold the component data while the downscaling is done
   //Size calculation from SHM. TODO: Figure out why
   int max_width = MAX(param.src_width, param.trgt_width);
   int max_height = MAX(param.trgt_height, param.src_height);
   int min_width = MIN(param.src_width, param.trgt_width);
   int min_height = MIN(param.trgt_height, param.src_height);
   int min_width_rnd16 = ((min_width + 15) >> 4) << 4;
   int min_height_rnd32 = ((min_height + 31) >> 5) << 5;
   int buffer_width = ((max_width * min_width_rnd16 + (min_width << 4) - 1) / (min_width << 4)) << 4;
   int buffer_height = ((max_height * min_height_rnd32 + (min_height << 4) - 1) / (min_height << 4)) << 4;;
   pic_buffer_t* buffer = newPictureBuffer(buffer_width, buffer_height, 1);

   //Downscale y-component
   //Target buffers for components
   pic_buffer_t* y = newPictureBuffer(param.trgt_width, param.trgt_height, 0); //TODO: Need to use max height/width?

   copyPictureBuffer(yuv->y, buffer, 1);
   downsample(buffer, &param);
   copyPictureBuffer(buffer, y, 0);

   //TODO: Maybe make factor choosing less broad and have something to do with actual type.
   //Downscale u
   //Check if 420 or 444 etc.
   int recalc = 0;
   if (yuv->y->width != yuv->u->width) {
     w_factor = yuv->y->width < yuv->u->width ? 1 : -1;
     recalc = 1;
   }
   if (yuv->y->height != yuv->u->height) {
     h_factor = yuv->y->height < yuv->u->height ? 1 : -1;
     recalc = 1;
   }
   if (recalc) {
     //Re-calculate parameters
     calculateParameters(&param, w_factor, h_factor);
   }

   pic_buffer_t* u = newPictureBuffer(param.trgt_width, param.trgt_height, 0);
   copyPictureBuffer(yuv->u, buffer, 1);
   downsample(buffer, &param);
   copyPictureBuffer(buffer, u, 0);

   //Downscale v
   //Check if v same size as u
   recalc = 0;
   w_factor = 0;
   h_factor = 0;
   if (yuv->v->width != yuv->u->width) {
     w_factor = yuv->v->width > yuv->u->width ? 1 : -1;
     recalc = 1;
   }
   if (yuv->v->height != yuv->u->height) {
     h_factor = yuv->v->height > yuv->u->width ? 1 : -1;
     recalc = 1;
   }
   if (recalc) {
     //Re-calculate parameters
     calculateParameters(&param, w_factor, h_factor);
   }

   pic_buffer_t* v = newPictureBuffer(param.trgt_width, param.trgt_height, 0);
   copyPictureBuffer(yuv->v, buffer, 1);
   downsample(buffer, &param);
   copyPictureBuffer(buffer, v, 0);

   //Make final yuv_buffer
   yuv_buffer_t* new_yuv = (yuv_buffer_t*)malloc(sizeof(yuv_buffer_t));
   new_yuv->y = y;
   new_yuv->u = u;
   new_yuv->v = v;

   //Deallocate buffer. TODO: Reuse buffer?
   deallocatePictureBuffer(buffer);
   
   return new_yuv;
 }
