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
#include "scaler-util.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#define DEFAULT_RESAMPLE_BLOCK_STEP_FUNC resampleBlockStep
#define DEFAULT_RESAMPLE_FUNC resample
#define ALT_RESAMPLE_FUNC resample2resampleBlockStep_default
#define OPAQUE_RESAMPLE_BLOCK_STEP_FUNC opaqueResampleBlockStep_adapter


pic_buffer_t* kvz_newPictureBuffer(int width, int height, int has_tmp_row)
{
  pic_buffer_t* buffer = (pic_buffer_t*)malloc(sizeof(pic_buffer_t));
  if (buffer == NULL) {
    return NULL; //TODO: Add error message?
  }

  //Allocate enough memory to fit a width-by-height picture
  buffer->data = (pic_data_t*)malloc(sizeof(pic_data_t) * width * height + SCALER_BUFFER_PADDING);

  buffer->width = width;
  buffer->height = height;

  //Initialize tmp_row or set as NULL
  if (has_tmp_row) {
    //Use max dim for size
    int max_dim = SCALER_MAX(width, height);
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
  yuv->format = format;
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

opaque_pic_buffer_t * kvz_newOpaquePictureBuffer(void *const data, int width, int height, int stride, const unsigned alloc_depth)
{
  opaque_pic_buffer_t *buffer = (opaque_pic_buffer_t*)malloc(sizeof(opaque_pic_buffer_t));
  if (buffer == NULL)
  {
    return NULL;
  }
  
  if (data == NULL && alloc_depth != 0 ) {
    buffer->data = malloc(width * height * alloc_depth);
  } else {
    buffer->data = data;
  }

  buffer->width = width;
  buffer->height = height;
  buffer->stride = stride;
  buffer->depth = alloc_depth;

  return buffer;
}

void kvz_setOpaqueYuvBuffer(opaque_yuv_buffer_t * buffer, void * const y_data, void * const u_data, void * const v_data, const unsigned depth)
{
  kvz_setOpaquePicBuffer(buffer->y, y_data, depth);
  kvz_setOpaquePicBuffer(buffer->u, u_data, depth);
  kvz_setOpaquePicBuffer(buffer->v, v_data, depth);
}

void kvz_setOpaquePicBuffer(opaque_pic_buffer_t * buffer, void * const data, const unsigned depth)
{
  buffer->data = data;
  buffer->depth = depth;
}

opaque_yuv_buffer_t * kvz_newOpaqueYuvBuffer(void * const y_data, void * const u_data, void * const v_data, int width, int height, int stride, chroma_format_t format, const unsigned alloc_depth)
{
  opaque_yuv_buffer_t *buffer = (opaque_yuv_buffer_t*)malloc(sizeof(opaque_yuv_buffer_t));
  if (buffer == NULL)
  {
    return NULL;
  }

  buffer->format = format;

  buffer->y = kvz_newOpaquePictureBuffer(y_data, width, height, stride, alloc_depth);

  int w_factor = 0;
  int h_factor = 0;

  switch (format) {
    case CHROMA_400:
    {
      //No chroma
      width = height = stride = 0;
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
  stride = stride >> w_factor;

  buffer->u = kvz_newOpaquePictureBuffer(u_data, width, height, stride, alloc_depth);
  buffer->v = kvz_newOpaquePictureBuffer(v_data, width, height, stride, alloc_depth);

  return buffer;
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
void kvz_deallocatePictureBuffer(pic_buffer_t* buffer)
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
  int max_dim_x = fill ? dst->width : SCALER_MIN(src->width, dst->width);
  int max_dim_y = fill ? dst->height : SCALER_MIN(src->height, dst->height);
  int min_dim_x = fill ? src->width : max_dim_x;
  int min_dim_y = fill ? src->height : max_dim_y;

  int dst_row = 0;
  int src_row = 0;

  //Copy loop
  for (int i = 0; i < max_dim_y; i++) {
    if (i < min_dim_y) {
      for (int j = 0; j < max_dim_x; j++) {
        //If outside min o_ind, copy adjacent value.
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
* \param src/dst_x is the x-coordinate for the sub-block (needs to be valid for both buffers).
* \param src/dst_y is the y-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_width is the width for the sub-block (needs to be valid for both buffers).
* \param block_height is the height for the sub-block (needs to be valid for both buffers).
*/
static void copyPictureBufferBlock(const pic_buffer_t* const src, const pic_buffer_t* const dst, const int src_x, const int src_y, const int dst_x, const int dst_y, const int block_width, const int block_height )
{
  for( int sy = src_y, dy = dst_y; sy < block_height && dy < block_height; sy++, dy++ ){
    for(int sx = src_x, dx = dst_x; sx < block_height && dx < block_height; sx++, dx++){
      dst->data[dx + dy*dst->width] = src->data[sx + sy*src->width];
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
* \param src/dst_x is the x-coordinate for the sub-block (needs to be valid for both buffers).
* \param src/dst_y is the y-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_width is the width for the sub-block (needs to be valid for both buffers).
* \param block_height is the height for the sub-block (needs to be valid for both buffers).
* \param w_factor is how much chroma sizes are scaled (width).
* \param h_factor is how much chroma sizes are scaled (heigth).
*/
static void copyYuvBufferBlock(const yuv_buffer_t* const src, const yuv_buffer_t* const dst, const int src_x, const int src_y, const int dst_x, const int dst_y, const int block_width, const int block_height, const int w_factor, const int h_factor)
{
  copyPictureBufferBlock(src->y, dst->y, src_x, src_y, dst_x, dst_y, block_width, block_height);
  copyPictureBufferBlock(src->u, dst->u, SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(block_width, w_factor), SCALER_SHIFT(block_height, h_factor));
  copyPictureBufferBlock(src->v, dst->v, SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(block_width, w_factor), SCALER_SHIFT(block_height, h_factor));
}


//Copy memory blocks between different types
static void copyMemBlock( void * const dst, const void * const src, const int dst_sizeof, const int src_sizeof, const int dst_stride, const int src_stride, const int dst_x, const int dst_y, const int src_x, const int src_y, const int block_width, const int block_height)
{
  //Cast to char pointers to allow indexing
  char * dst_char = &((char *)dst)[(dst_x + dst_y * dst_stride) * dst_sizeof];
  const char * src_char = &((const char *)src)[(src_x + src_y * src_stride) * src_sizeof];

  assert(sizeof(char)==1); //May not work if char is not one byte

  //Loop over rows
  for(int y = 0; y < block_height; y++){
    //Init dst to zero
    memset(dst_char, 0, block_width * dst_sizeof);
    //Copy row
    for(int x = 0; x < block_width; x++){
      memcpy(&dst_char[x*dst_sizeof], &src_char[x*src_sizeof], SCALER_MIN(dst_sizeof,src_sizeof));
    }
    
    //Move to next row
    dst_char += dst_stride * dst_sizeof;
    src_char += src_stride * src_sizeof;
  }
}

/**
* \brief Copies data from a sub-block from one opaque buffer to the other.
* \param src is the source buffer
* \param dst is the destination buffer
* \param src/dst_x is the x-coordinate for the sub-block (needs to be valid for both buffers).
* \param src/dst_y is the y-coordinate for the sub-block (needs to be valid for both buffers).
* \param block_width is the width for the sub-block (needs to be valid for both buffers).
* \param block_height is the height for the sub-block (needs to be valid for both buffers).
* \param w_factor is how much chroma sizes are scaled (width).
* \param h_factor is how much chroma sizes are scaled (heigth).
*/
static void copyOpaqueYuvBufferBlock(const opaque_yuv_buffer_t* const src, const opaque_yuv_buffer_t* const dst, const int src_x, const int src_y, const int dst_x, const int dst_y, const int block_width, const int block_height, const int w_factor, const int h_factor)
{
  copyMemBlock(dst->y->data, src->y->data, dst->y->depth, src->y->depth, dst->y->stride, src->y->stride, dst_x, dst_y, src_x, src_y, block_width, block_height);
  copyMemBlock(dst->u->data, src->u->data, dst->u->depth, src->u->depth, dst->u->stride, src->u->stride, SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_SHIFT(block_width, w_factor), SCALER_SHIFT(block_height, h_factor));
  copyMemBlock(dst->v->data, src->v->data, dst->v->depth, src->v->depth, dst->v->stride, src->v->stride, SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_SHIFT(block_width, w_factor), SCALER_SHIFT(block_height, h_factor));
}


void kvz_copy_uint8_block_to_YuvBuffer(const yuv_buffer_t* dst, const uint8_t* const y, const uint8_t* const u, const uint8_t* const v, const int luma_stride, const int dst_x, const int dst_y, const int src_x, const int src_y, const int block_width, const int block_height, const int w_factor, const int h_factor){
  //TODO: Sanity checks

  copyMemBlock(dst->y->data, y, sizeof(pic_data_t), sizeof(uint8_t), dst->y->width, luma_stride, dst_x, dst_y, src_x, src_y, block_width, block_height);
  copyMemBlock(dst->u->data, u, sizeof(pic_data_t), sizeof(uint8_t), dst->u->width, SCALER_SHIFT(luma_stride, w_factor), SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));
  copyMemBlock(dst->v->data, v, sizeof(pic_data_t), sizeof(uint8_t), dst->v->width, SCALER_SHIFT(luma_stride, w_factor), SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));
}


void kvz_copy_YuvBuffer_block_to_uint8(uint8_t* const y, uint8_t* const u, uint8_t* const v, const int luma_stride, const yuv_buffer_t * const src, const int dst_x, const int dst_y, const int src_x, const int src_y, const int block_width, const int block_height, const int w_factor, const int h_factor){
  //TODO: Sanity checks
  copyMemBlock(y, src->y->data, sizeof(uint8_t), sizeof(pic_data_t), luma_stride, src->y->width, dst_x, dst_y, src_x, src_y, block_width, block_height);
  copyMemBlock(u, src->u->data, sizeof(uint8_t), sizeof(pic_data_t), SCALER_SHIFT(luma_stride, w_factor), src->u->width, SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));
  copyMemBlock(v, src->v->data, sizeof(uint8_t), sizeof(pic_data_t), SCALER_SHIFT(luma_stride, w_factor), src->v->width, SCALER_SHIFT(dst_x, w_factor), SCALER_SHIFT(dst_y, h_factor), SCALER_SHIFT(src_x, w_factor), SCALER_SHIFT(src_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));
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
  pic_buffer_t* ret = kvz_newPictureBuffer(pic->width, pic->height, pic->tmp_row != NULL);
  if (ret == NULL) {
    return NULL;
  }

  memcpy(ret->data, pic->data, pic->width * pic->height * sizeof(pic_data_t));

  if (pic->tmp_row) {
    int tmp_size = SCALER_MAX(pic->width, pic->height);
    memcpy(ret->tmp_row, pic->tmp_row, tmp_size);
  }

  return ret;
}

yuv_buffer_t* kvz_cloneYuvBuffer(const yuv_buffer_t* const yuv)
{
  yuv_buffer_t* ret = malloc(sizeof(yuv_buffer_t));
  if (ret == NULL || yuv == NULL) {
    free(ret);
    return NULL; //TODO: Add error message?
  }
  ret->y = clonePictureBuffer(yuv->y);
  ret->u = clonePictureBuffer(yuv->u);
  ret->v = clonePictureBuffer(yuv->v);

  return ret;
}


opaque_yuv_buffer_t * kvz_copyOpaqueYuvBuffer(const opaque_yuv_buffer_t * const yuv)
{
  if (yuv == NULL) {
    return NULL; //TODO: Add error message?
  }
  
  return kvz_newOpaqueYuvBuffer(yuv->y->data, yuv->u->data, yuv->v->data, yuv->y->width, yuv->y->height, yuv->y->stride, yuv->format, yuv->y->depth);
}

void kvz_deallocateYuvBuffer(yuv_buffer_t* yuv)
{
  if (yuv == NULL) return;

  kvz_deallocatePictureBuffer(yuv->y);
  kvz_deallocatePictureBuffer(yuv->u);
  kvz_deallocatePictureBuffer(yuv->v);

  free(yuv);
}

void kvz_deallocateOpaquePictureBuffer(opaque_pic_buffer_t * buffer, const int free_buffer)
{
  //Don't deallocate data here if free buffer is not set
  if (free_buffer)
  {
    free(buffer->data);
  }
  free(buffer);
}

void kvz_deallocateOpaqueYuvBuffer(opaque_yuv_buffer_t * yuv, const int free_buffer)
{
  if (yuv == NULL) return;

  kvz_deallocateOpaquePictureBuffer(yuv->y, free_buffer);
  kvz_deallocateOpaquePictureBuffer(yuv->u, free_buffer);
  kvz_deallocateOpaquePictureBuffer(yuv->v, free_buffer);

  free(yuv);
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

  const int *filter_hor;
  const int filter_size_hor = prepareFilter(&filter_hor, is_upscaling, is_luma, hor_filter);
  const int *filter_ver;
  const int filter_size_ver = prepareFilter(&filter_ver, is_upscaling, is_luma, ver_filter);

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
      //int size = getFilter(&filter, is_upscaling, is_luma, phase, hor_filter);
      int size = filter_size_hor;
      filter = &getFilterCoeff(filter_hor, size, phase, 0);

      //Apply filter
      tmp_row[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = SCALER_CLIP(ref_pos + k - (size >> 1) + 1, 0, src_width - 1);
        tmp_row[j] += filter[k] * src_row[m];
      }
    }
    //Copy tmp row to data
    memcpy(src_row, tmp_row, sizeof(pic_data_t) * trgt_width);
  }

  pic_data_t* tmp_col = tmp_row; //rename for clarity

  // Vertical resampling. TODO: Why this order? should loop over rows not cols?
  for (int i = 0; i < trgt_width; i++) {
    pic_data_t* src_col = &buffer->data[i];
    for (int j = 0; j < trgt_height; j++) {
      //Calculate ref pos
      int ref_pos_16 = (int)((unsigned int)(j * param->scale_y + param->add_y) >> shift_y) - param->delta_y;
      int phase = ref_pos_16 & 15;
      int ref_pos = ref_pos_16 >> 4;

      //Choose filter
      const int* filter;
      //int size = getFilter(&filter, is_upscaling, is_luma, phase, ver_filter);
      int size = filter_size_ver;
      filter = &getFilterCoeff(filter_ver, size, phase, 0);

      //Apply filter
      tmp_col[j] = 0;
      for (int k = 0; k < size; k++) {
        int m = SCALER_CLIP(ref_pos + k - (size >> 1) + 1, 0, src_height - 1);
        tmp_col[j] += filter[k] * src_col[m * buffer->width];
      }
      //TODO: Why? Filter coefs summ up to 128 applied 2x 128*128= 2^14
      //TODO: Why?  Filter coefs summ up to 64 applied 2x 64*64= 2^12
      //Scale values back down
      tmp_col[j] = is_upscaling ? (tmp_col[j] + 2048) >> 12 : (tmp_col[j] + 8192) >> 14;
    }

    //Clip and move to buffer data
    for (int n = 0; n < trgt_height; n++) {
      src_col[n * buffer->width] = SCALER_CLIP(tmp_col[n], 0, 255);
    }
  }
}

//Do only the vertical/horizontal resampling step in the given block
static void resampleBlockStep(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height,  const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical){
  //TODO: Add cropping etc.

  //Choose best filter to use when downsampling
  //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
  int filter_phase = 0;

  const int src_size = is_vertical ? param->src_height + param->src_padding_y : param->src_width + param->src_padding_x;
  const int trgt_size = is_vertical ? param->rnd_trgt_height : param->rnd_trgt_width;

  if (!is_upscaling) {
    int crop_size = src_size - ( is_vertical ? param->bottom_offset : param->right_offset); //- param->left_offset/top_offset;
    if (4 * crop_size > 15 * trgt_size)
      filter_phase = 7;
    else if (7 * crop_size > 20 * trgt_size)
      filter_phase = 6;
    else if (2 * crop_size > 5 * trgt_size)
      filter_phase = 5;
    else if (1 * crop_size > 2 * trgt_size)
      filter_phase = 4;
    else if (3 * crop_size > 5 * trgt_size)
      filter_phase = 3;
    else if (4 * crop_size > 5 * trgt_size)
      filter_phase = 2;
    else if (19 * crop_size > 20 * trgt_size)
      filter_phase = 1;
  }

  const int shift = (is_vertical ? param->shift_y : param->shift_x) - 4;
  const int scale = is_vertical ? param->scale_y : param->scale_x;
  const int add = is_vertical ? param->add_y : param->add_x;
  const int delta = is_vertical ? param->delta_y : param->delta_x;

  //Set loop parameters based on the resampling dir
  const int *filter;
  const int filter_size = prepareFilter(&filter, is_upscaling, is_luma, filter_phase);
  const int outer_init = is_vertical ? 0 : block_x;
  const int outer_bound = is_vertical ? filter_size : block_x + block_width;
  const int inner_init = is_vertical ? block_x : 0;
  const int inner_bound = is_vertical ? block_x + block_width : filter_size;
  const int s_stride = is_vertical ? src_buffer->width : 1; //Multiplier to s_ind

  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < (block_y + block_height); y++) {

    pic_data_t* src = is_vertical ? src_buffer->data : &src_buffer->data[y * src_buffer->width];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width];

    //Outer loop:
    //  if is_vertical -> loop over k (filter inds)
    //  if !is_vertical -> loop over x (block width)
    for (int o_ind = outer_init; o_ind < outer_bound; o_ind++) {
      
      const int t_ind = is_vertical ? y : o_ind; //trgt_buffer row/col index for cur resampling dir

      //Calculate reference position in src pic
      const int ref_pos_16 = (int)((unsigned int)(t_ind * scale + add) >> shift) - delta;
      const int phase = ref_pos_16 & 15;
      const int ref_pos = ref_pos_16 >> 4;
      
      //Inner loop:
      //  if is_vertical -> loop over x (block width)
      //  if !is_vertical -> loop over k (filter inds)-
      for (int i_ind = inner_init; i_ind < inner_bound; i_ind++) {

        const int f_ind = is_vertical ? o_ind : i_ind; //Filter index
        const int t_col = is_vertical ? i_ind : o_ind; //trgt_buffer column
        
        //Choose filter
        //const int *filter;
        //const int f_size = getFilter(&filter, is_upscaling, is_luma, phase, filter_phase);

        //Set trgt buffer val to zero on first loop over filter
        if( f_ind == 0 ){
          trgt_row[t_col + trgt_offset] = 0;
        }

        const int s_ind = SCALER_CLIP(ref_pos + f_ind - (filter_size >> 1) + 1, 0, src_size - 1); //src_buffer row/col index for cur resampling dir

        //Move src pointer to correct position (correct column in vertical resampling)
        pic_data_t *src_col = src + (is_vertical ? i_ind : 0);
        trgt_row[t_col + trgt_offset] += getFilterCoeff(filter,filter_size,phase,f_ind) * src_col[s_ind * s_stride + src_offset];

        //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
        if (is_vertical && o_ind == outer_bound - 1) {
          trgt_row[t_col + trgt_offset] = SCALER_CLIP(is_upscaling ? (trgt_row[t_col + trgt_offset] + 2048) >> 12 : (trgt_row[t_col + trgt_offset] + 8192) >> 14, 0, 255);
        }
      }
    }
  }

}

#define OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(src_type, trgt_type, tmp_trgt_type) do { \
for (int y = block_y; y < (block_y + block_height); y++) {\
  src_type *src = is_vertical ? (src_type *)(src_buffer->data) : &((src_type *)(src_buffer->data))[y * src_buffer->stride];\
  trgt_type *trgt_row = &((trgt_type *)(trgt_buffer->data))[y * trgt_buffer->stride];\
  tmp_trgt_type *tmp_trgt_row = &tmp_trgt_data[(y - block_y) * tmp_trgt_data_stride];\
  for (int o_ind = outer_init; o_ind < outer_bound; o_ind++) {\
    const int t_ind = is_vertical ? y : o_ind;\
    const int ref_pos_16 = (int)((unsigned int)(t_ind * scale + add) >> shift) - delta;\
    const int phase = ref_pos_16 & 15;\
    const int ref_pos = ref_pos_16 >> 4;\
    for (int i_ind = inner_init; i_ind < inner_bound; i_ind++) {\
      const int f_ind = is_vertical ? o_ind : i_ind;\
      const int t_col = is_vertical ? i_ind : o_ind;\
      if (f_ind == 0) {\
        tmp_trgt_row[t_col - block_x + tmp_trgt_offset] = 0;\
      }\
      const int s_ind = SCALER_CLIP(ref_pos + f_ind - (filter_size >> 1) + 1, 0, src_size - 1);\
      src_type *src_col = src + (is_vertical ? i_ind : 0);\
      tmp_trgt_row[t_col - block_x + tmp_trgt_offset] += (tmp_trgt_type)(getFilterCoeff(filter, filter_size, phase, f_ind) * (pic_data_t)(src_col[s_ind * s_stride + src_offset]));\
      if (is_vertical && o_ind == outer_bound - 1) {\
        trgt_row[t_col + trgt_offset] = (trgt_type)(SCALER_CLIP(is_upscaling ? (tmp_trgt_row[t_col - block_x + tmp_trgt_offset] + 2048) >> 12 : (tmp_trgt_row[t_col - block_x + tmp_trgt_offset] + 8192) >> 14, 0, 255));\
      }\
    }\
  }\
}\
}while(0)

/**
*  \brief Do resampling on opaque data buffers. Supported bit-depths: 
*     horizontal step: {8,16,32}-bit -> {16,32}-bit
*     vertical step: {16,32}-bit -> {8,16,32}-bit
*/
static void opaqueResampleBlockStep_adapter(const opaque_pic_buffer_t* const src_buffer, const opaque_pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical) 
{  
  //Based on the src and trgt depths select the relevant function
  if (trgt_buffer->depth == sizeof(pic_data_t) && src_buffer->depth == sizeof(pic_data_t)) {
    //Use the basic algorithm
    resampleBlockStep((const pic_buffer_t *)src_buffer, (const pic_buffer_t *)trgt_buffer, src_offset, trgt_offset, block_x, block_y, block_width, block_height, param, is_upscaling, is_luma, is_vertical);
    return;
  }

  //TODO: Add cropping etc.

  //Choose best filter to use when downsampling
  //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
  int filter_phase = 0;

  const int src_size = is_vertical ? param->src_height + param->src_padding_y : param->src_width + param->src_padding_x;
  const int trgt_size = is_vertical ? param->rnd_trgt_height : param->rnd_trgt_width;

  if (!is_upscaling) {
    int crop_size = src_size - (is_vertical ? param->bottom_offset : param->right_offset); //- param->left_offset/top_offset;
    if (4 * crop_size > 15 * trgt_size)
      filter_phase = 7;
    else if (7 * crop_size > 20 * trgt_size)
      filter_phase = 6;
    else if (2 * crop_size > 5 * trgt_size)
      filter_phase = 5;
    else if (1 * crop_size > 2 * trgt_size)
      filter_phase = 4;
    else if (3 * crop_size > 5 * trgt_size)
      filter_phase = 3;
    else if (4 * crop_size > 5 * trgt_size)
      filter_phase = 2;
    else if (19 * crop_size > 20 * trgt_size)
      filter_phase = 1;
  }

  const int shift = (is_vertical ? param->shift_y : param->shift_x) - 4;
  const int scale = is_vertical ? param->scale_y : param->scale_x;
  const int add = is_vertical ? param->add_y : param->add_x;
  const int delta = is_vertical ? param->delta_y : param->delta_x;

  //Set loop parameters based on the resampling dir
  const int *filter;
  const int filter_size = prepareFilter(&filter, is_upscaling, is_luma, filter_phase);
  const int outer_init = is_vertical ? 0 : block_x;
  const int outer_bound = is_vertical ? filter_size : block_x + block_width;
  const int inner_init = is_vertical ? block_x : 0;
  const int inner_bound = is_vertical ? block_x + block_width : filter_size;
  const int s_stride = is_vertical ? src_buffer->stride : 1; //Multiplier to s_ind

  

  if (is_vertical && src_buffer->depth >= sizeof(short)) {

    //Set tmp target and trgt buffer the same by default
    pic_data_t *tmp_trgt_data = (pic_data_t *)malloc(sizeof(pic_data_t) * block_width * block_height);
    int tmp_trgt_data_stride = 0;
    int tmp_trgt_offset = block_width;

    if (src_buffer->depth == sizeof(pic_data_t)) {
      //Handle 32-bit input buffer case
      if (trgt_buffer->depth == sizeof(short)) {
        //Handle 16-bit output buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(pic_data_t, unsigned short, pic_data_t);
      }
      else if (trgt_buffer->depth == sizeof(char)) {
        //Handle 8-bit output buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(pic_data_t, unsigned char, pic_data_t);
      }
      else {
        //No valid handling for the given depth
        assert(0);
      }
    }
    else if (src_buffer->depth == sizeof(short)) {
      //Handle 16-bit input buffer case
      if (trgt_buffer->depth == sizeof(pic_data_t)) {
        //Handle 32-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned short, pic_data_t, pic_data_t);
      }
      else if (trgt_buffer->depth == sizeof(short)) {
        //Handle 16-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned short, unsigned short, pic_data_t);
      } 
      else if (trgt_buffer->depth == sizeof(char)) {
        //Handle 8-bit output buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned short, unsigned char, pic_data_t);
      }
      else {
        //No valid handling for the given depth
        assert(0);
      }
    } 
    else {
      //No valid handling for the given depth
      assert(0);
    }

    free(tmp_trgt_data);
  }
  else if (!is_vertical)
  {
    //Set tmp target and trgt buffer the same by default
    int tmp_trgt_data_stride = trgt_buffer->stride;
    int tmp_trgt_offset = trgt_offset;

    if (trgt_buffer->depth == sizeof(pic_data_t)) {
      //Handle 32-bit out buffer case
      pic_data_t *tmp_trgt_data = (pic_data_t *)trgt_buffer->data + block_y * trgt_buffer->stride + block_x;

      if (src_buffer->depth == sizeof(short)) {
        //Handle 16-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned short, pic_data_t, pic_data_t);
      } 
      else if (src_buffer->depth == sizeof(char)) {
        //Handle 8-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned char, pic_data_t, pic_data_t);
      } 
      else {
        //8-bit trgt buffer not supported
        assert(0);
      }
    }
    else if (trgt_buffer->depth == sizeof(short)) {
      //Handle 16-bit out buffer case
      unsigned short *tmp_trgt_data = (unsigned short *)trgt_buffer->data + block_y * trgt_buffer->stride + block_x;

      if (src_buffer->depth == sizeof(pic_data_t)) {
        //Handle 32-bit in buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(pic_data_t, unsigned short, unsigned short);
      } 
      else if (src_buffer->depth == sizeof(short)) {
        //Handle 16-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned short, unsigned short, unsigned short);
      } 
      else if (src_buffer->depth == sizeof(char)) {
        //Handle 8-bit input buffer case
        OPAQUE_RESAMPLE_BLOCK_STEP_TYPE_MACRO(unsigned char, unsigned short, unsigned short);
      } 
      else {
        //8-bit trgt buffer not supported
        assert(0);
      }
    } else {
      //No valid handling for the given depths
      assert(0);
    }
  }
  else {
    //Not 8-bit vertical src buffer not supported
    assert(0);
  }
}

//Do the resampling in one pass using 2D-convolution.
static void resampleBlock( const pic_buffer_t* const src_buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const pic_buffer_t *const trgt_buffer, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height )
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
  int trgt_stride = trgt_buffer->width;

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

      //Convolution kernel, where (o_ind,y)<=>(0,0)
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
      //trgt offset is used to reposition the block in trgt in case buffer size is not target image size
      trgt[x + y*trgt_stride + trgt_offset] = clip(is_upscaling ? (new_val + 2048) >> 12 : (new_val + 8192) >> 14, 0, 255); //TODO: account for different bit dept / different filters etc
    }  
  }
}

//Calculate scaling parameters and update param. Factor determines if certain values are 
// divided eg. with chroma. 0 for no factor and -1 for halving stuff and 1 for doubling etc.
//Calculations from SHM
static void calculateParameters(scaling_parameter_t* const param, const int w_factor, const int h_factor, const int is_chroma)
{
  if(param->is_calculated){
    return;
  }

  //First shift widths/height by an appropriate factor
  param->src_width = SCALER_SHIFT(param->src_width, w_factor);
  param->src_height = SCALER_SHIFT(param->src_height, h_factor);
  param->trgt_width = SCALER_SHIFT(param->trgt_width, w_factor);
  param->trgt_height = SCALER_SHIFT(param->trgt_height, h_factor);
  param->scaled_src_width = SCALER_SHIFT(param->scaled_src_width, w_factor);
  param->scaled_src_height = SCALER_SHIFT(param->scaled_src_height, h_factor);
  param->rnd_trgt_width = SCALER_SHIFT(param->rnd_trgt_width, w_factor);
  param->rnd_trgt_height = SCALER_SHIFT(param->rnd_trgt_height, h_factor);
  param->src_padding_x = SCALER_SHIFT(param->src_padding_x, w_factor);
  param->src_padding_y = SCALER_SHIFT(param->src_padding_y, h_factor);
  param->trgt_padding_x = SCALER_SHIFT(param->trgt_padding_x, w_factor);
  param->trgt_padding_y = SCALER_SHIFT(param->trgt_padding_y, h_factor);

  //Calculate sample positional parameters
  param->right_offset = param->src_width - param->scaled_src_width; //- left_offset
  param->bottom_offset = param->src_height - param->scaled_src_height; //- top_offset

  //TODO: Make dependant on width/height?
  param->shift_x = SCALER_SHIFT_CONST;
  param->shift_y = SCALER_SHIFT_CONST;

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

  param->is_calculated = 1;
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
  //TODO: Why SCALER_MAX? Try using src
  int scaled_src_width = param.src_width;//SCALER_MAX(param.src_width, param.trgt_width); //Min/max of source/target values
  int scaled_src_height = param.src_height;//SCALER_MAX(param.src_height, param.trgt_size); //Min/max of source/target values
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
  param.rnd_trgt_height = param.trgt_height;//((param.trgt_size + (1 << 4) - 1) >> 4) << 4; //Round to multiple of 16

  //Round to multiple of 2
  //TODO: Why SCALER_MAX? Try using src
  int scaled_src_width = param.src_width;//SCALER_MAX(param.src_width, param.trgt_width); //Min/max of source/target values
  int scaled_src_height = param.src_height;//SCALER_MAX(param.src_height, param.trgt_size); //Min/max of source/target values
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

yuv_buffer_t* kvz_yuvScaling_adapter(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst, resample_func *const resample_func)
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
    param.is_calculated = 0;
    calculateParameters(&param, 0, 0, 0);
  }

  //Check if we need to allocate a yuv buffer for the new image or re-use dst.
  //Make sure the sizes match
  if (dst == NULL || dst->y->width != param.trgt_width || dst->y->height != param.trgt_height
    || dst->u->width != SCALER_SHIFT(param.trgt_width, w_factor) || dst->u->height != SCALER_SHIFT(param.trgt_height, h_factor)
    || dst->v->width != SCALER_SHIFT(param.trgt_width, w_factor) || dst->v->height != SCALER_SHIFT(param.trgt_height, h_factor)) {

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
  int max_width = SCALER_MAX(param.src_width + param.src_padding_x, param.trgt_width);
  int max_height = SCALER_MAX(param.src_height + param.src_padding_y, param.trgt_height);
  int min_width = SCALER_MIN(param.src_width, param.trgt_width);
  int min_height = SCALER_MIN(param.src_height, param.trgt_height);
  int min_width_rnd16 = ((min_width + 15) >> 4) << 4;
  int min_height_rnd32 = ((min_height + 31) >> 5) << 5;
  int buffer_width = ((max_width * min_width_rnd16 + (min_width << 4) - 1) / (min_width << 4)) << 4;
  int buffer_height = ((max_height * min_height_rnd32 + (min_height << 4) - 1) / (min_height << 4)) << 4;;
  pic_buffer_t* buffer = kvz_newPictureBuffer(buffer_width, buffer_height, 1);


  /*==========Start Resampling=============*/
  //Resample y
  copyPictureBuffer(yuv->y, buffer, 1);
  /*resample*/resample_func(buffer, &param, is_upscaling, 1);
  copyPictureBuffer(buffer, dst->y, 0);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      param.is_calculated = 0;
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //Resample u
    copyPictureBuffer(yuv->u, buffer, 1);
    /*resample*/resample_func(buffer, &param, is_upscaling, 0);
    copyPictureBuffer(buffer, dst->u, 0);

    //Resample v
    copyPictureBuffer(yuv->v, buffer, 1);
    /*resample*/resample_func(buffer, &param, is_upscaling, 0);
    copyPictureBuffer(buffer, dst->v, 0);
  }

  //Deallocate buffer
  kvz_deallocatePictureBuffer(buffer);

  return dst;
}

yuv_buffer_t* kvz_yuvScaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param,
                         yuv_buffer_t* dst)
{
  return kvz_yuvScaling_adapter(yuv, base_param, dst, kvz_default_resample_func);
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
    param.is_calculated = 0;
    calculateParameters(&param, w_factor, h_factor, 0);
  }

  //Check if we need to allocate a yuv buffer for the new image or re-use dst.
  //Make sure the sizes match
  if (dst == NULL || dst->y->width != param.trgt_width || dst->y->height != param.trgt_height
    || dst->u->width != SCALER_SHIFT(param.trgt_width, w_factor) || dst->u->height != SCALER_SHIFT(param.trgt_height, h_factor)
    || dst->v->width != SCALER_SHIFT(param.trgt_width, w_factor) || dst->v->height != SCALER_SHIFT(param.trgt_height, h_factor)) {

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
  /*int max_width = SCALER_MAX(param.src_width, param.trgt_width);
  int max_height = SCALER_MAX(param.src_height, param.trgt_size);
  int min_width = SCALER_MIN(param.src_width, param.trgt_width);
  int min_height = SCALER_MIN(param.src_height, param.trgt_size);
  int min_width_rnd16 = ((min_width + 15) >> 4) << 4;
  int min_height_rnd32 = ((min_height + 31) >> 5) << 5;
  int buffer_width = ((max_width * min_width_rnd16 + (min_width << 4) - 1) / (min_width << 4)) << 4;
  int buffer_height = ((max_height * min_height_rnd32 + (min_height << 4) - 1) / (min_height << 4)) << 4;;
  pic_buffer_t* buffer = kvz_newPictureBuffer(buffer_width, buffer_height, 1);*/
  //TODO: Clean up this part and implement properly
  //param.rnd_trgt_height = param.trgt_size;
  //param.rnd_trgt_width = param.trgt_width;
  yuv_buffer_t* buffer = is_upscaling ? dst : yuv;//malloc(sizeof(pic_buffer_t)); //Choose bigger buffer
  if (buffer->y->tmp_row == NULL) {
    buffer->y->tmp_row = malloc(sizeof(pic_data_t) * (SCALER_MAX(buffer->y->width, buffer->y->height)));
  }
  if (buffer->u->tmp_row == NULL) {
    buffer->u->tmp_row = malloc(sizeof(pic_data_t) * (SCALER_MAX(buffer->u->width, buffer->u->height)));
  }
  if (buffer->v->tmp_row == NULL) {
    buffer->v->tmp_row = malloc(sizeof(pic_data_t) * (SCALER_MAX(buffer->v->width, buffer->v->height)));
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
      param.is_calculated = 0;
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
  //kvz_deallocatePictureBuffer(buffer);

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
  if (yuv == NULL || yuv->y->width < param.src_width + param.src_padding_x || yuv->y->height < param.src_height + param.src_padding_y || yuv->u->width < SCALER_SHIFT(param.src_width + param.src_padding_x, w_factor) || yuv->u->height < SCALER_SHIFT(param.src_height + param.src_padding_y, w_factor) || yuv->v->width < SCALER_SHIFT(param.src_width + param.src_padding_x, w_factor) || yuv->v->height < SCALER_SHIFT(param.src_height + param.src_padding_y, w_factor)) {
    fprintf(stderr, "Source buffer smaller than specified in the scaling parameters.\n");
    return 0;
  }
 
  //Calculate a dst offset depending on wheather dst is the whole image or just the block
  // if dst is smaller than the specified trgt size, the scaled block is written starting from (0,0)
  // if dst is the size of the specified trgt, the scaled block is written starting from (block_x,block_y)
  int dst_offset_luma = 0;
  int dst_offset_chroma = 0;
  if (dst == NULL || dst->y->width < param.trgt_width || dst->y->height < param.trgt_height
    || dst->u->width < SCALER_SHIFT(param.trgt_width, w_factor) || dst->u->height < SCALER_SHIFT(param.trgt_height, h_factor)
    || dst->v->width < SCALER_SHIFT(param.trgt_width, w_factor) || dst->v->height < SCALER_SHIFT(param.trgt_height, h_factor)) {

    //Check that dst is large enough to hold the block
    if (dst == NULL || dst->y->width < block_width || dst->y->height < block_height
      || dst->u->width < SCALER_SHIFT(block_width, w_factor) || dst->u->height < SCALER_SHIFT(block_height, h_factor)
      || dst->v->width < SCALER_SHIFT(block_width, w_factor) || dst->v->height < SCALER_SHIFT(block_height, h_factor)) {
      fprintf(stderr, "Destination buffer not large enough to hold block\n");
      return 0;
    }
    //Set dst offset so that the block is written to the correct pos
    dst_offset_luma = block_x + block_y * dst->y->width;
    dst_offset_chroma = SCALER_SHIFT(block_x, w_factor) + SCALER_SHIFT(block_y, h_factor) * dst->u->width;
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
    fprintf(stderr, "Both dimensions need to be either upscaled or downscaled");
    return 0;
  }
  if (is_equal_height && is_equal_width) {
    //If equal just copy block from src
    copyYuvBufferBlock(yuv, dst, block_x, block_y, dst_offset_luma != 0 ? 0 : block_x, dst_offset_luma != 0 ? 0 : block_y, block_width, block_height, w_factor, h_factor);
    return 1;
  }
  if (is_downscaled_width || is_downscaled_height) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }
  /*=================================*/

  /*==========Start Resampling=============*/

  //Resample y
  resampleBlock(yuv->y, &param, is_upscaling, 1, dst->y, -dst_offset_luma, block_x, block_y, block_width, block_height);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      param.is_calculated = 0;
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //In order to scale blocks not divisible by 2 correctly, need to do some tricks
    //Resample u
    resampleBlock(yuv->u, &param, is_upscaling, 0, dst->u, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));

    //Resample v
    resampleBlock(yuv->v, &param, is_upscaling, 0, dst->v, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor));
  }

  return 1;
}

static void blockScalingSrcRange( int range[2], const int scale, const int add, const int shift, const int delta, const int block_low, const int block_high, const int src_size )
{
  //Check if equal size
  if(scale == SCALER_UNITY_SCALE_CONST){
    range[0] = block_low;
    range[1] = block_high;
    return;
  }

  //Get filter size
  int size = scale < SCALER_UNITY_SCALE_CONST ? sizeof(lumaUpFilter[0]) / sizeof(lumaUpFilter[0][0]) : sizeof(downFilter[0][0]) / sizeof(downFilter[0][0][0]);

  //Calculate lower bound
  range[0] = ((int)((unsigned int)((block_low * scale + add) >> (shift - 4)) - delta) >> 4) - (size >> 1) + 1;

  //Calculate upper bound
  range[1] = ((int)((unsigned int)((block_high * scale + add) >> (shift - 4)) - delta) >> 4) - (size >> 1) + size;

  //clip the ranges so that they are within the pic
  range[0] = SCALER_CLIP(range[0], 0, src_size - 1);
  range[1] = SCALER_CLIP(range[1], 0, src_size - 1);

}

//Do validity checks and calculate offset parameters and scaling dir if needed
static int blockStepScalingChecks( const yuv_buffer_t *const dst, const yuv_buffer_t *const src, const scaling_parameter_t *const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical, int *dst_offset_luma, int *dst_offset_chroma, int *src_offset_luma, int *src_offset_chroma, int *scaling_dir, int *w_factor_out, int *h_factor_out)
{
  //Check that block parameters are valid
  int width_bound = base_param->trgt_width;
  int height_bound = is_vertical ? base_param->trgt_height : base_param->src_height + base_param->src_padding_y;
  if (block_x < 0 || block_y < 0 || block_x + block_width > width_bound || block_y + block_height > height_bound) {
    fprintf(stderr, "Specified block outside given target picture size.\n");
    return 0;
  }

  //Initialize basic parameters
  const scaling_parameter_t param = *base_param;

  //How much to scale the luma sizes to get the chroma sizes
  int w_factor = 0;
  int h_factor = 0;
  switch (param.chroma) {
  case CHROMA_400: {
    //No chroma
    assert(src->u->height == 0 && src->u->width == 0 && src->v->height == 0 && src->v->width == 0);
    break;
  }
  case CHROMA_420: {
    assert(src->u->height == (src->y->height >> 1) && src->u->width == (src->y->width >> 1)
      && src->v->height == (src->y->height >> 1) && src->v->width == (src->y->width >> 1));
    w_factor = -1;
    h_factor = -1;
    break;
  }
  case CHROMA_422: {
    assert(src->u->height == (src->y->height) && src->u->width == (src->y->width >> 1)
      && src->v->height == (src->y->height) && src->v->width == (src->y->width >> 1));
    w_factor = -1;
    break;
  }
  case CHROMA_444: {
    assert(src->u->height == (src->y->height) && src->u->width == (src->y->width)
      && src->v->height == (src->y->height) && src->v->width == (src->y->width));
    break;
  }
  default:
    assert(0); //Unsupported chroma type
  }

  if (w_factor_out) *w_factor_out = w_factor;
  if (h_factor_out) *h_factor_out = h_factor;

  //Calculate a src offset depending on wheather src is the whole image or just the block
  // if src is smaller than the specified src size, the src buffer is intepreted to contain the area specified by kvz_blockScaling*Range.
  // if src is the size of the specified src, the src buffer is accessed in the area specified by kvz_blockScaling*Range.
  width_bound = is_vertical ? param.trgt_width : param.src_width + param.src_padding_x;
  height_bound = param.src_height + param.src_padding_y;
  if (src == NULL || src->y->width < width_bound || src->y->height < height_bound || src->u->width < SCALER_SHIFT(width_bound, w_factor) || src->u->height < SCALER_SHIFT(height_bound, w_factor) || src->v->width < SCALER_SHIFT(width_bound, w_factor) || src->v->height < SCALER_SHIFT(height_bound, w_factor)) {

    //Get src range needed for scaling
    int range[4];
    kvz_blockScalingSrcWidthRange(range, base_param, block_x, block_width);
    kvz_blockScalingSrcHeightRange(range + 2, base_param, block_y, block_height);
    width_bound = is_vertical ? block_width : range[1] - range[0];
    height_bound = is_vertical ? range[3] - range[2] : block_height;

    //Check that src is large enough to hold the block
    if (src == NULL || src->y->width < width_bound || src->y->height < height_bound
      || src->u->width < SCALER_SHIFT(width_bound, w_factor) || src->u->height < SCALER_SHIFT(height_bound, h_factor)
      || src->v->width < SCALER_SHIFT(width_bound, w_factor) || src->v->height < SCALER_SHIFT(height_bound, h_factor)) {
      fprintf(stderr, "Source buffer smaller than specified in the scaling parameters.\n");
      return 0;
    }
    //Set src offset so that the block is written to the correct pos
    if(src_offset_luma) *src_offset_luma = is_vertical ? block_x + range[2] * src->y->width : range[0] + block_y * src->y->width;
    if(src_offset_chroma) *src_offset_chroma = is_vertical ? SCALER_SHIFT(block_x, w_factor) + SCALER_SHIFT(range[2], h_factor) * src->u->width : SCALER_SHIFT(range[0], w_factor) + SCALER_SHIFT(block_y, h_factor) * src->u->width;
  }

  //Calculate a dst offset depending on wheather dst is the whole image or just the block
  // if dst is smaller than the specified trgt size, the scaled block is written starting from (0,0)
  // if dst is the size of the specified trgt, the scaled block is written starting from (block_x,block_y)
  width_bound = param.trgt_width;
  height_bound = is_vertical ? param.trgt_height : param.src_height + param.src_padding_y;
  if (dst == NULL || dst->y->width < width_bound || dst->y->height < height_bound
    || dst->u->width < SCALER_SHIFT(width_bound, w_factor) || dst->u->height < SCALER_SHIFT(height_bound, h_factor)
    || dst->v->width < SCALER_SHIFT(width_bound, w_factor) || dst->v->height < SCALER_SHIFT(height_bound, h_factor)) {

    //Check that dst is large enough to hold the block
    if (dst == NULL || dst->y->width < block_width || dst->y->height < block_height
      || dst->u->width < SCALER_SHIFT(block_width, w_factor) || dst->u->height < SCALER_SHIFT(block_height, h_factor)
      || dst->v->width < SCALER_SHIFT(block_width, w_factor) || dst->v->height < SCALER_SHIFT(block_height, h_factor)) {
      fprintf(stderr, "Destination buffer not large enough to hold block\n");
      return 0;
    }
    //Set dst offset so that the block is written to the correct pos
    if (dst_offset_luma) *dst_offset_luma = block_x + block_y * dst->y->width;
    if (dst_offset_chroma) *dst_offset_chroma = SCALER_SHIFT(block_x, w_factor) + SCALER_SHIFT(block_y, h_factor) * dst->u->width;
  }

  //Calculate if we are upscaling or downscaling
  int is_downscaled_width = param.src_width > param.trgt_width;
  int is_downscaled_height = param.src_height > param.trgt_height;
  int is_equal_width = param.src_width == param.trgt_width;
  int is_equal_height = param.src_height == param.trgt_height;

  //both dimensions need to be either up- or downscaled
  if ((is_downscaled_width && !is_downscaled_height && !is_equal_height) ||
    (is_downscaled_height && !is_downscaled_width && !is_equal_width)) {
    fprintf(stderr, "Both dimensions need to be either upscaled or downscaled\n");
    return 0;
  }

  if (scaling_dir == NULL) {
    //Don't set scaling dir
    return 1;
  }

  if (is_equal_height && is_equal_width) {
    *scaling_dir = 0; //No scaling
  }
  else if (is_downscaled_width || is_downscaled_height) {
    //Atleast one dimension is downscaled
    *scaling_dir = -1;
  } else {
    *scaling_dir = 1; //upscaling
  }

  return 1;
}

void kvz_blockScalingSrcWidthRange(int range[2], const scaling_parameter_t * const base_param, const int block_x, const int block_width)
{
  //Calculate parameters
  calculateParameters((scaling_parameter_t*)base_param, 0, 0, 0);

  blockScalingSrcRange(range, base_param->scale_x, base_param->add_x, base_param->shift_x, base_param->delta_x, block_x, block_x + block_width - 1, base_param->src_width + base_param->src_padding_x);
}

int kvz_opaqueYuvBlockStepScaling_adapter(opaque_yuv_buffer_t* const dst, const opaque_yuv_buffer_t* const src, const scaling_parameter_t* const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical, opaque_resample_block_step_func * const resample_func)
{
  //Do validity checks
  int scaling_dir = 0;
  int src_offset_luma = 0;
  int src_offset_chroma = 0;
  int dst_offset_luma = 0;
  int dst_offset_chroma = 0;

  int w_factor = 0;
  int h_factor = 0;

  if (!blockStepScalingChecks((yuv_buffer_t *)dst, (yuv_buffer_t *)src, base_param, block_x, block_y, block_width, block_height, is_vertical, &dst_offset_luma, &dst_offset_chroma, &src_offset_luma, &src_offset_chroma, &scaling_dir, &w_factor, &h_factor))
  {
    return 0;
  }

  if (scaling_dir == 0) {
    //If equal just copy block from src
    copyOpaqueYuvBufferBlock(src, dst, src_offset_luma != 0 ? 0 : block_x, src_offset_luma != 0 ? 0 : block_y, dst_offset_luma != 0 ? 0 : block_x, dst_offset_luma != 0 ? 0 : block_y, block_width, block_height, w_factor, h_factor);
    return 1;
  }

  int is_upscaling = 1;
  if (scaling_dir == -1) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }
  
  scaling_parameter_t param = *base_param;

  //Do resampling

  //Resample y
  /*resampleBlockStep*/resample_func(src->y, dst->y, -src_offset_luma, -dst_offset_luma, block_x, block_y, block_width, block_height, &param, is_upscaling, 1, is_vertical);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      param.is_calculated = 0;
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //In order to scale blocks not divisible by 2 correctly, need to do some tricks
    //Resample u
    /*resampleBlockStep*/resample_func(src->u, dst->u, -src_offset_chroma, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor), &param, is_upscaling, 0, is_vertical);

    //Resample v
    /*resampleBlockStep*/resample_func(src->v, dst->v, -src_offset_chroma, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor), &param, is_upscaling, 0, is_vertical);
  }

  return 1;
}

int kvz_opaqueYuvBlockStepScaling(opaque_yuv_buffer_t* const dst, const opaque_yuv_buffer_t* const src, const scaling_parameter_t* const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical)
{
  return kvz_opaqueYuvBlockStepScaling_adapter(dst, src, base_param, block_x, block_y, block_width, block_height, is_vertical, kvz_opaque_block_step_resample_func);
}

void kvz_blockScalingSrcHeightRange(int range[2], const scaling_parameter_t * const base_param, const int block_y, const int block_height)
{
  //Calculate parameters
  calculateParameters((scaling_parameter_t*)base_param, 0, 0, 0);

  blockScalingSrcRange(range, base_param->scale_y, base_param->add_y, base_param->shift_y, base_param->delta_y, block_y, block_y + block_height - 1, base_param->src_height + base_param->src_padding_y);
}

// Do block scaling in one direction. yuv buffer should not be modified.
int kvz_yuvBlockStepScaling_adapter(yuv_buffer_t* const dst, const yuv_buffer_t* const src, const scaling_parameter_t* const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical, resample_block_step_func * const resample_func)
{
  /*========== Basic Initialization ==============*/

    //Do validity checks
  int scaling_dir = 0;
  int src_offset_luma = 0;
  int src_offset_chroma = 0;
  int dst_offset_luma = 0;
  int dst_offset_chroma = 0;

  int w_factor = 0;
  int h_factor = 0;

  if (!blockStepScalingChecks(dst, src, base_param, block_x, block_y, block_width, block_height, is_vertical, &dst_offset_luma, &dst_offset_chroma, &src_offset_luma, &src_offset_chroma, &scaling_dir, &w_factor, &h_factor))
  {
    return 0;
  }

  if (scaling_dir == 0) {
    //If equal just copy block from src
    copyYuvBufferBlock(src, dst, src_offset_luma != 0 ? 0 : block_x, src_offset_luma != 0 ? 0 : block_y, dst_offset_luma != 0 ? 0 : block_x, dst_offset_luma != 0 ? 0 : block_y, block_width, block_height, w_factor, h_factor);
    return 1;
  }

  int is_upscaling = 1;
  if (scaling_dir == -1) {
    //Atleast one dimension is downscaled
    is_upscaling = 0;
  }

  scaling_parameter_t param = *base_param;
  /*=================================*/

  /*==========Start Resampling=============*/

  //Resample y
  /*resampleBlockStep*/resample_func(src->y, dst->y, -src_offset_luma, -dst_offset_luma, block_x, block_y, block_width, block_height, &param, is_upscaling, 1, is_vertical);

  //Skip chroma if CHROMA_400
  if (param.chroma != CHROMA_400) {
    //If chroma size differs from luma size, we need to recalculate the parameters
    if (h_factor != 0 || w_factor != 0) {
      param.is_calculated = 0;
      calculateParameters(&param, w_factor, h_factor, 1);
    }

    //In order to scale blocks not divisible by 2 correctly, need to do some tricks
    //Resample u
    /*resampleBlockStep*/resample_func(src->u, dst->u, -src_offset_chroma, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor), &param, is_upscaling, 0, is_vertical);

    //Resample v
    /*resampleBlockStep*/resample_func(src->v, dst->v, -src_offset_chroma, -dst_offset_chroma, SCALER_SHIFT(block_x, w_factor), SCALER_SHIFT(block_y, h_factor), SCALER_ROUND_SHIFT(block_width, w_factor), SCALER_ROUND_SHIFT(block_height, h_factor), &param, is_upscaling, 0, is_vertical);
  }

  return 1;
}

int kvz_yuvBlockStepScaling(yuv_buffer_t * const dst, const yuv_buffer_t * const src, const scaling_parameter_t * const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical)
{
  return kvz_yuvBlockStepScaling_adapter(dst, src, base_param, block_x, block_y, block_width, block_height, is_vertical, kvz_default_block_step_resample_func);
}

static void resample2resampleBlockStep_default(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
{
  pic_buffer_t* tmp = kvz_newPictureBuffer(param->trgt_width + param->trgt_padding_x, param->src_height + param->src_padding_y, 0);

  //Vertical resampling
  kvz_default_block_step_resample_func(buffer, tmp, 0, 0, 0, 0, param->trgt_width + param->trgt_padding_x, param->src_height + param->src_padding_y, param, is_upscaling, is_luma, 0);

  //Horizontal resampling
  kvz_default_block_step_resample_func(tmp, buffer, 0, 0, 0, 0, param->trgt_width + param->trgt_padding_x, param->trgt_height + param->trgt_padding_y, param, is_upscaling, is_luma, 1);

  kvz_deallocatePictureBuffer(tmp);
}

//Set the default resample function
resample_block_step_func *const kvz_default_block_step_resample_func = &DEFAULT_RESAMPLE_BLOCK_STEP_FUNC;
resample_func *const kvz_default_resample_func = &DEFAULT_RESAMPLE_FUNC;
resample_func *const kvz_alt_resample_func  = &ALT_RESAMPLE_FUNC;
opaque_resample_block_step_func *const kvz_opaque_block_step_resample_func = &OPAQUE_RESAMPLE_BLOCK_STEP_FUNC;