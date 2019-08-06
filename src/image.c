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

#include "image.h"

#include <limits.h>
#include <stdlib.h>

#include "strategies/strategies-ipol.h"
#include "strategies/strategies-picture.h"
#include "threads.h"

#include "strategies/strategies-resample.h"

/**
* \brief Allocate a new image with 420.
* This function signature is part of the libkvz API.
* \return image pointer or NULL on failure
*/
kvz_picture * kvz_image_alloc_420(const int32_t width, const int32_t height)
{
  return kvz_image_alloc(KVZ_CSP_420, width, height);
}

/**
 * \brief Allocate a new image.
 * \return image pointer or NULL on failure
 */
kvz_picture * kvz_image_alloc(enum kvz_chroma_format chroma_format, const int32_t width, const int32_t height)
{
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  kvz_picture *im = MALLOC(kvz_picture, 1);
  if (!im) return NULL;

  unsigned int luma_size = width * height;
  unsigned chroma_sizes[] = { 0, luma_size / 4, luma_size / 2, luma_size };
  unsigned chroma_size = chroma_sizes[chroma_format];

  im->chroma_format = chroma_format;

  //Allocate memory
  im->fulldata = MALLOC(kvz_pixel, (luma_size + 2 * chroma_size + SCALER_BUFFER_PADDING));
  if (!im->fulldata) {
    free(im);
    return NULL;
  }

  im->base_image = im;
  im->refcount = 1; //We give a reference to caller
  im->width = width;
  im->height = height;
  im->stride = width;
  im->chroma_format = chroma_format;

  im->y = im->data[COLOR_Y] = &im->fulldata[0];

  if (chroma_format == KVZ_CSP_400) {
    im->u = im->data[COLOR_U] = NULL;
    im->v = im->data[COLOR_V] = NULL;
  } else {
    im->u = im->data[COLOR_U] = &im->fulldata[luma_size];
    im->v = im->data[COLOR_V] = &im->fulldata[luma_size + chroma_size];
  }

  im->pts = 0;
  im->dts = 0;

  im->interlacing = KVZ_INTERLACING_NONE;

  return im;
}

/**
 * \brief Free an image.
 *
 * Decrement reference count of the image and deallocate associated memory
 * if no references exist any more.
 *
 * \param im image to free
 */
void kvz_image_free(kvz_picture *const im)
{
  if (im == NULL) return;

  int32_t new_refcount = KVZ_ATOMIC_DEC(&(im->refcount));
  if (new_refcount > 0) {
    // There are still references so we don't free the data yet.
    return;
  }

  if (im->base_image != im) {
    // Free our reference to the base image.
    kvz_image_free(im->base_image);
  } else {
    free(im->fulldata);
  }

  // Make sure freed data won't be used.
  im->base_image = NULL;
  im->fulldata = NULL;
  im->y = im->u = im->v = NULL;
  im->data[COLOR_Y] = im->data[COLOR_U] = im->data[COLOR_V] = NULL;
  free(im);
}

/**
 * \brief Get a new pointer to an image.
 *
 * Increment reference count and return the image.
 */
kvz_picture *kvz_image_copy_ref(kvz_picture *im)
{
  int32_t new_refcount = KVZ_ATOMIC_INC(&im->refcount);
  // The caller should have had another reference and we added one
  // reference so refcount should be at least 2.
  assert(new_refcount >= 2);
  return im;
}

kvz_picture *kvz_image_make_subimage(kvz_picture *const orig_image,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height)
{
  // Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  assert((x_offset % 2) == 0);
  assert((y_offset % 2) == 0);

  assert(x_offset + width <= orig_image->width);
  assert(y_offset + height <= orig_image->height);

  kvz_picture *im = MALLOC(kvz_picture, 1);
  if (!im) return NULL;

  im->base_image = kvz_image_copy_ref(orig_image->base_image);
  im->refcount = 1; // We give a reference to caller
  im->width = width;
  im->height = height;
  im->stride = orig_image->stride;
  im->chroma_format = orig_image->chroma_format;

  im->y = im->data[COLOR_Y] = &orig_image->y[x_offset + y_offset * orig_image->stride];
  if (orig_image->chroma_format != KVZ_CSP_400) {
    im->u = im->data[COLOR_U] = &orig_image->u[x_offset / 2 + y_offset / 2 * orig_image->stride / 2];
    im->v = im->data[COLOR_V] = &orig_image->v[x_offset / 2 + y_offset / 2 * orig_image->stride / 2];
  }

  im->pts = 0;
  im->dts = 0;

  return im;
}

yuv_t * kvz_yuv_t_alloc(int luma_size, int chroma_size)
{
  yuv_t *yuv = (yuv_t *)malloc(sizeof(*yuv));
  yuv->size = luma_size;

  // Get buffers with separate mallocs in order to take advantage of
  // automatic buffer overrun checks.
  yuv->y = (kvz_pixel *)malloc(luma_size * sizeof(*yuv->y));
  if (chroma_size == 0) {
    yuv->u = NULL;
    yuv->v = NULL;
  } else {
    yuv->u = (kvz_pixel *)malloc(chroma_size * sizeof(*yuv->u));
    yuv->v = (kvz_pixel *)malloc(chroma_size * sizeof(*yuv->v));
  }
  
  return yuv;
}

void kvz_yuv_t_free(yuv_t *yuv)
{
  if (yuv) {
    FREE_POINTER(yuv->y);
    FREE_POINTER(yuv->u);
    FREE_POINTER(yuv->v);
  }
  FREE_POINTER(yuv);
}

hi_prec_buf_t * kvz_hi_prec_buf_t_alloc(int luma_size)
{
  // Get buffers with separate mallocs in order to take advantage of
  // automatic buffer overrun checks.
  hi_prec_buf_t *yuv = (hi_prec_buf_t *)malloc(sizeof(*yuv));
  yuv->y = (int16_t *)malloc(luma_size * sizeof(*yuv->y));
  yuv->u = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->u));
  yuv->v = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->v));
  yuv->size = luma_size;

  return yuv;
}

void kvz_hi_prec_buf_t_free(hi_prec_buf_t * yuv)
{
  free(yuv->y);
  free(yuv->u);
  free(yuv->v);
  free(yuv);
}


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
static unsigned cor_sad(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride)
{
  kvz_pixel ref = *ref_data;
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref);
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
static unsigned ver_sad(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref_data[x]);
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
static unsigned hor_sad(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride, unsigned ref_stride)
{
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref_data[y * ref_stride]);
    }
  }

  return sad;
}


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
static unsigned image_interpolated_sad(const kvz_picture *pic, const kvz_picture *ref,
                                 int pic_x, int pic_y, int ref_x, int ref_y,
                                 int block_width, int block_height)
{
  kvz_pixel *pic_data, *ref_data;

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
  pic_data = &pic->y[pic_y * pic->stride + pic_x];
  ref_data = &ref->y[ref_y * ref->stride + ref_x];

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
                      &ref_data[top * ref->stride + left],
                      left, top, pic->stride);
    result += ver_sad(&pic_data[left],
                      &ref_data[top * ref->stride + left],
                      block_width - left, top, pic->stride);
    result += hor_sad(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride + left],
                      left, block_height - top, pic->stride, ref->stride);
    result += kvz_reg_sad(&pic_data[top * pic->stride + left],
                      &ref_data[top * ref->stride + left],
                      block_width - left, block_height - top, pic->stride, ref->stride);
  } else if (top && right) {
    result += ver_sad(pic_data,
                      &ref_data[top * ref->stride],
                      block_width - right, top, pic->stride);
    result += cor_sad(&pic_data[block_width - right],
                      &ref_data[top * ref->stride + (block_width - right - 1)],
                      right, top, pic->stride);
    result += kvz_reg_sad(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride],
                      block_width - right, block_height - top, pic->stride, ref->stride);
    result += hor_sad(&pic_data[top * pic->stride + (block_width - right)],
                      &ref_data[top * ref->stride + (block_width - right - 1)],
                      right, block_height - top, pic->stride, ref->stride);
  } else if (bottom && left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height - bottom, pic->stride, ref->stride);
    result += kvz_reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height - bottom, pic->stride, ref->stride);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      left, bottom, pic->stride);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->stride + left],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      block_width - left, bottom, pic->stride);
  } else if (bottom && right) {
    result += kvz_reg_sad(pic_data,
                      ref_data,
                      block_width - right, block_height - bottom, pic->stride, ref->stride);
    result += hor_sad(&pic_data[block_width - right],
                      &ref_data[block_width - right - 1],
                      right, block_height - bottom, pic->stride, ref->stride);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride],
                      block_width - right, bottom, pic->stride);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->stride + block_width - right],
                      &ref_data[(block_height - bottom - 1) * ref->stride + block_width - right - 1],
                      right, bottom, pic->stride);
  } else if (top) {
    result += ver_sad(pic_data,
                      &ref_data[top * ref->stride],
                      block_width, top, pic->stride);
    result += kvz_reg_sad(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride],
                      block_width, block_height - top, pic->stride, ref->stride);
  } else if (bottom) {
    result += kvz_reg_sad(pic_data,
                      ref_data,
                      block_width, block_height - bottom, pic->stride, ref->stride);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride],
                      block_width, bottom, pic->stride);
  } else if (left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height, pic->stride, ref->stride);
    result += kvz_reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height, pic->stride, ref->stride);
  } else if (right) {
    result += kvz_reg_sad(pic_data,
                      ref_data,
                      block_width - right, block_height, pic->stride, ref->stride);
    result += hor_sad(&pic_data[block_width - right],
                      &ref_data[block_width - right - 1],
                      right, block_height, pic->stride, ref->stride);
  } else {
    result += kvz_reg_sad(pic_data, ref_data, block_width, block_height, pic->stride, ref->stride);
  }

  return result;
}


/**
* \brief Calculate interpolated SAD between two blocks.
*
* \param pic        Image for the block we are trying to find.
* \param ref        Image where we are trying to find the block.
*
* \returns          Sum of absolute differences
*/
unsigned kvz_image_calc_sad(const kvz_picture *pic,
                            const kvz_picture *ref,
                            int pic_x,
                            int pic_y,
                            int ref_x,
                            int ref_y,
                            int block_width,
                            int block_height)
{
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);

  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const kvz_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];
    const kvz_pixel *ref_data = &ref->y[ref_y * ref->stride + ref_x];
    return kvz_reg_sad(pic_data, ref_data, block_width, block_height, pic->stride, ref->stride)>>(KVZ_BIT_DEPTH-8);
  } else {
    // Call a routine that knows how to interpolate pixels outside the frame.
    return image_interpolated_sad(pic, ref, pic_x, pic_y, ref_x, ref_y, block_width, block_height) >> (KVZ_BIT_DEPTH - 8);
  }
}


/**
* \brief Calculate interpolated SATD between two blocks.
*
* \param pic        Image for the block we are trying to find.
* \param ref        Image where we are trying to find the block.
*/
unsigned kvz_image_calc_satd(const kvz_picture *pic,
                             const kvz_picture *ref,
                             int pic_x,
                             int pic_y,
                             int ref_x,
                             int ref_y,
                             int block_width,
                             int block_height)
{
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);

  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const kvz_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];
    const kvz_pixel *ref_data = &ref->y[ref_y * ref->stride + ref_x];
    return kvz_satd_any_size(block_width,
                             block_height,
                             pic_data,
                             pic->stride,
                             ref_data,
                             ref->stride) >> (KVZ_BIT_DEPTH - 8);
  } else {
    // Extrapolate pixels from outside the frame.
    kvz_extended_block block;
    kvz_get_extended_block(pic_x,
                           pic_y,
                           ref_x - pic_x,
                           ref_y - pic_y,
                           0,
                           0,
                           ref->y,
                           ref->width,
                           ref->height,
                           0,
                           block_width,
                           block_height,
                           &block);

    const kvz_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];

    unsigned satd = kvz_satd_any_size(block_width,
                                      block_height,
                                      pic_data,
                                      pic->stride,
                                      block.buffer,
                                      block.stride) >> (KVZ_BIT_DEPTH - 8);

    if (block.malloc_used) {
      FREE_POINTER(block.buffer);
    }

    return satd;
  }
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
#define BLIT_PIXELS_CASE(n) case n:\
  for (y = 0; y < n; ++y) {\
    memcpy(&dst[y*dst_stride], &orig[y*orig_stride], n * sizeof(kvz_pixel));\
  }\
  break;

void kvz_pixels_blit(const kvz_pixel * const orig, kvz_pixel * const dst,
                         const unsigned width, const unsigned height,
                         const unsigned orig_stride, const unsigned dst_stride)
{
  unsigned y;
  //There is absolutely no reason to have a width greater than the source or the destination stride.
  assert(width <= orig_stride);
  assert(width <= dst_stride);

#ifdef CHECKPOINTS
  char *buffer = malloc((3 * width + 1) * sizeof(char));
  for (y = 0; y < height; ++y) {
    int p;
    for (p = 0; p < width; ++p) {
      sprintf((buffer + 3*p), "%02X ", orig[y*orig_stride]);
    }
    buffer[3*width] = 0;
    CHECKPOINT("kvz_pixels_blit_avx2: %04d: %s", y, buffer);
  }
  FREE_POINTER(buffer);
#endif //CHECKPOINTS

  if (width == orig_stride && width == dst_stride) {
    memcpy(dst, orig, width * height * sizeof(kvz_pixel));
    return;
  }

  int nxn_width = (width == height) ? width : 0;
  switch (nxn_width) {
    BLIT_PIXELS_CASE(4)
    BLIT_PIXELS_CASE(8)
    BLIT_PIXELS_CASE(16)
    BLIT_PIXELS_CASE(32)
    BLIT_PIXELS_CASE(64)
  default:

    if (orig == dst) {
      //If we have the same array, then we should have the same stride
      assert(orig_stride == dst_stride);
      return;
    }
    assert(orig != dst || orig_stride == dst_stride);

    for (y = 0; y < height; ++y) {
      memcpy(&dst[y*dst_stride], &orig[y*orig_stride], width * sizeof(kvz_pixel));
    }
    break;
  }
}

// ***********************************************
  // Modified for SHVC
//Deallocates the given parameters
void kvz_image_scaler_worker( void *opaque_param)
{
  kvz_image_scaling_parameter_t *in_param = opaque_param;
  kvz_picture *pic_in = in_param->pic_in;
  kvz_picture *pic_out = in_param->pic_out;
  const scaling_parameter_t *const param = in_param->param;


  yuv_buffer_t* src_pic = kvz_newYuvBuffer_padded_uint8(pic_in->y, pic_in->u, pic_in->v,
                                                        param->src_width + param->src_padding_x,
                                                        param->src_height + param->src_padding_y,
                                                        pic_in->stride, param->chroma, 0);
  //yuv_buffer_t* src_pic = newYuvBuffer_uint8(pic_in->y, pic_in->u, pic_in->v, pic_in->width, pic_in->height, param->chroma, 0);
  
  yuv_buffer_t* trgt_pic = kvz_yuvScaling_adapter(src_pic, param, NULL, kvz_resample);//kvz_yuvScaling(src_pic, param, NULL );

  if (trgt_pic != NULL) {
    //Get out_img padding
    uint8_t padding_x = param->trgt_padding_x;
    uint8_t padding_y = param->trgt_padding_y;

    //Copy other information
    pic_out->dts = pic_in->dts;
    pic_out->pts = pic_in->pts;
    pic_out->interlacing = pic_in->interlacing;

    int chroma_shift = param->chroma == CHROMA_444 ? 0 : 1;
    pic_data_t* comp_list[] = { trgt_pic->y->data, trgt_pic->u->data, trgt_pic->v->data };
    int stride_list[] = { trgt_pic->y->width, trgt_pic->u->width, trgt_pic->v->width };
    int height_list[] = { trgt_pic->y->height, trgt_pic->u->height, trgt_pic->v->height };
    int padd_x[] = { padding_x, padding_x >> chroma_shift, padding_x >> chroma_shift };
    int padd_y[] = { padding_y, padding_y >> chroma_shift, padding_y >> chroma_shift };
    assert(sizeof(kvz_pixel) == sizeof(char)); //Image copy (memset) only works if the pixels are the same size as char 

    //Loop over components
    for (int comp = 0, i = 0; comp < sizeof(comp_list) / sizeof(*comp_list); comp++) {
      int comp_size = height_list[comp] * stride_list[comp];
      int pic_out_stride = pic_out->stride >> (comp < 1 ? 0 : chroma_shift);
      for (int src_ind = 0; src_ind < comp_size; i++, src_ind++) {
        //TODO: go over src image correctly
        //TODO: Make a better loop
        //Copy value normally
        pic_out->fulldata[i] = comp_list[comp][src_ind];

        if (padding_x != 0 && (src_ind % stride_list[comp] == stride_list[comp] - 1)) { //Padd end of row by copying last pixel
          memset(pic_out->fulldata + i + 1, pic_out->fulldata[i], padd_x[comp]);
          i += padd_x[comp];
        }
      }
      if (padd_y[comp] != 0) { //Padd image with lines copied from the prev row
        for (int j = 0; j < padd_y[comp]; j++) {
          memcpy(pic_out->fulldata + i, pic_out->fulldata + i - pic_out_stride, pic_out_stride);
          i += pic_out_stride;
        }
      }
    }
  }

  //Do deallocation
  kvz_deallocateYuvBuffer(src_pic);
  kvz_deallocateYuvBuffer(trgt_pic);
  kvz_image_free(pic_in);
  kvz_image_free(pic_out);
  free(in_param);
}

void kvz_block_scaler_worker(void * opaque_param)
{
  kvz_image_scaling_parameter_t *in_param = opaque_param;
  kvz_picture * const pic_in = in_param->pic_in;
  kvz_picture *pic_out = in_param->pic_out;
  const scaling_parameter_t *const param = in_param->param;


  yuv_buffer_t* src_pic = kvz_newYuvBuffer_padded_uint8(pic_in->y, pic_in->u, pic_in->v,
    param->src_width + param->src_padding_x,
    param->src_height + param->src_padding_y,
    pic_in->stride, param->chroma, 0);
  //yuv_buffer_t* src_pic = newYuvBuffer_uint8(pic_in->y, pic_in->u, pic_in->v, pic_in->width, pic_in->height, param->chroma, 0);

  yuv_buffer_t* trgt_pic = kvz_newYuvBuffer(in_param->block_width, in_param->block_height, param->chroma, 0);
    
  if( !kvz_yuvBlockScaling(src_pic, param, trgt_pic, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height) )
  {
    //TODO: Do error stuff?
    kvz_deallocateYuvBuffer(src_pic);
    kvz_deallocateYuvBuffer(trgt_pic);
    kvz_image_free(pic_in);
    kvz_image_free(pic_out);
    free(in_param);
    return;
  }

  if (trgt_pic != NULL) {
    
    //Copy other information
    pic_out->dts = pic_in->dts;
    pic_out->pts = pic_in->pts;
    pic_out->interlacing = pic_in->interlacing;

    int chroma_shift = param->chroma == CHROMA_444 ? 0 : 1;
    pic_buffer_t* comp_list[] = { trgt_pic->y, trgt_pic->u, trgt_pic->v };
    
    //Loop over components
    for (int c = 0; c < sizeof(comp_list)/sizeof(*comp_list); c++) {
      int pic_out_stride = pic_out->stride >> (c < 1 ? 0 : chroma_shift);
      int pic_out_offset_x = in_param->block_x >> (c < 1 ? 0 : chroma_shift);
      int pic_out_offset_y = in_param->block_y >> (c < 1 ? 0 : chroma_shift);
      
      //Copy block back to pic_out
      for( int y = 0; y < comp_list[c]->height; y++)
      {
        for( int x = 0; x < comp_list[c]->width; x++)
        {
          pic_out->data[c][(pic_out_offset_x + x) + (pic_out_offset_y + y) * pic_out_stride] = comp_list[c]->data[x + y * comp_list[c]->width];
        }
      }
    }
  }

  //Do deallocation
  kvz_deallocateYuvBuffer(src_pic);
  kvz_deallocateYuvBuffer(trgt_pic);
  kvz_image_free(pic_in);
  kvz_image_free(pic_out);
  free(in_param);

}

//Debug stuff for printing thread info
#if 0 && defined(linux)
#include <sys/types.h>
#include <sys/syscall.h>
#define PRINT_TID_JOB_INFO(x,y,w,h,dir) fprintf(stderr, "TID: %ld, pos: (%d,%d), size: (%d,%d), dir: %d\n", syscall(SYS_gettid), x, y, w, h, dir)
#define PRINT_JOB_EXTRA_INFO(msg,dx,dy,sx,sy,w,h) fprintf(stderr, "  %s: dst_pos: (%d,%d), src_pos: (%d,%d), size: (%d,%d)\n", msg, dx, dy, sx, sy, w, h)
#else
#define PRINT_TID_JOB_INFO(x,y,w,h,dir)
#define PRINT_JOB_EXTRA_INFO(msg,dx,dy,sx,sy,w,h)
#endif

/** \brief Handle hor/ver scaling steps
*  If tiles not used:
*    If pic_in is given, copy relevant block to src_buffer and run horizontal scaling step
*    If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out
*    If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
*  If tiles used:
*    If pic_in is given, copy relevant tile to src_buffer and run horizontal scaling step; src_buffer is filled starting from (0,0)
*    If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out; trgt_buffer is indexed starting from (0,0)
*    If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
*/
/*void kvz_block_step_scaler_worker(void * opaque_param)
{
  kvz_image_scaling_parameter_t *in_param = opaque_param;
  kvz_picture * const pic_in = in_param->pic_in;
  kvz_picture * const pic_out = in_param->pic_out;
  const scaling_parameter_t *const param = in_param->param;

  //TODO: account for chroma format properly
  int w_factor = -1;
  int h_factor = -1;

  PRINT_TID_JOB_INFO(in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, pic_in != NULL ? 1 : 0);

  //Hor Scaling
  if (pic_in != NULL) {
    
    int range[4] = {0, 0, 0, 0};
    if (in_param->use_tiles){
      kvz_blockScalingSrcWidthRange(range, param, in_param->block_x, in_param->block_width);    
    } else {
      //Get range that needs to be copied from pic_in
      //range[0:1] is prev blocks range and range[2:3] is the new block range
      //Only copy pixels not in range[0:1] (already copied by previous workers).
      //Getting the correct block_x of the previous block does not matter as only range[1] is used
      if (in_param->block_x - in_param->block_width < 0) {
        kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
        range[1] = range[2] - 1;
      }
      else {
        kvz_blockScalingSrcWidthRange(range, param, in_param->block_x - in_param->block_width, in_param->block_width);
        kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
      }
    }

    int cp_block_x = in_param->use_tiles ? range[0] : (range[1] + 1);
    int cp_block_y = in_param->block_y;
    int cp_block_width = (in_param->use_tiles ? range[1] : range[3]) - cp_block_x + 1;
    int cp_block_height = in_param->block_height;
    
    int hor_block_y = in_param->block_y;
    int hor_block_height = in_param->block_height;

    if (pic_out != NULL) {
      if (in_param->use_tiles){
        kvz_blockScalingSrcHeightRange(range, param, in_param->block_y, in_param->block_height);
        cp_block_y  = hor_block_y = range[0];
        cp_block_height = hor_block_height = range[1] - hor_block_y + 1;
      } else {
        //Do the same procedure as with horizontal range
        if (in_param->block_y - in_param->block_height < 0) {
          kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
          range[1] = range[2] - 1;
        }
        else {
          kvz_blockScalingSrcHeightRange(range, param, in_param->block_y - in_param->block_height, in_param->block_height);
          kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
        }
        cp_block_y = range[1] + 1;
        cp_block_height = range[3] - cp_block_y + 1;
        hor_block_y = range[2];
        hor_block_height = range[3] - range[2] + 1;
      }
    }

    //When using tiles, copy from in_pic to the src buffer (src buffer should hold only one tile and start from indexing (0,0))
    int cp_dst_x = in_param->use_tiles ? 0 : cp_block_x;
    int cp_dst_y = in_param->use_tiles ? 0 : cp_block_y;

    //Copy from in_pic to the src buffer
    kvz_copy_uint8_block_to_YuvBuffer(in_param->src_buffer, pic_in->y, pic_in->u, pic_in->v, pic_in->stride, cp_dst_x, cp_dst_y, cp_block_x, cp_block_y, cp_block_width, cp_block_height, w_factor, h_factor);

    PRINT_JOB_EXTRA_INFO("Copy to src buffer", cp_dst_x, cp_dst_y, cp_block_x, cp_block_y, cp_block_width, cp_block_height);
    PRINT_JOB_EXTRA_INFO("Hor scaling", in_param->block_x, hor_block_y, 0, 0, in_param->block_width, hor_block_height);

    //If both ver and hor done at the same time interpred in_param->block_y/height as the final output block and so we need to do hor scaling in the approriate range to accomodate the final block
    if (!kvz_yuvBlockStepScaling_adapter(in_param->ver_tmp_buffer, in_param->src_buffer, param, in_param->block_x, hor_block_y, in_param->block_width, hor_block_height, 0, kvz_resample_block_step)) {
      //TODO: Do error stuff?
      kvz_image_free(pic_in);
      kvz_image_free(pic_out);
      free(in_param);
      return;
    }
  }

  //Ver scaling
  if (pic_out != NULL) {
    //Do ver scaling step
    if (!kvz_yuvBlockStepScaling_adapter(in_param->trgt_buffer, in_param->ver_tmp_buffer, param, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, 1, kvz_resample_block_step)) {
      //TODO: Do error stuff?
      kvz_image_free(pic_in);
      kvz_image_free(pic_out);
      free(in_param);
      return;
    }

    const int dst_x = in_param->block_x;
    const int dst_y = in_param->block_y;
    const int src_x = in_param->use_tiles ? 0 : in_param->block_x;
    const int src_y = in_param->use_tiles ? 0 : in_param->block_y;

    //Copy results to pic_out
    kvz_copy_YuvBuffer_block_to_uint8(pic_out->y, pic_out->u, pic_out->v, pic_out->stride, in_param->trgt_buffer, dst_x, dst_y, src_x, src_y, in_param->block_width, in_param->block_height, w_factor, h_factor);

  }

  //Do deallocation
  kvz_image_free(pic_in);
  kvz_image_free(pic_out);
  free(in_param);
}*/

/** \brief Handle hor/ver scaling steps
*  If tiles not used:
*    If pic_in is given, copy relevant block to src_buffer and run horizontal scaling step
*    If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out
*    If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
*  If tiles used:
*    If pic_in is given, copy relevant tile to src_buffer and run horizontal scaling step; src_buffer is filled starting from (0,0)
*    If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out; trgt_buffer is indexed starting from (0,0)
*    If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
*/
void kvz_opaque_block_step_scaler_worker(void * opaque_param)
{
  kvz_image_scaling_parameter_t *in_param = opaque_param;
  kvz_picture * const pic_in = in_param->pic_in;
  kvz_picture * const pic_out = in_param->pic_out;
  const scaling_parameter_t *const param = in_param->param;

  PRINT_TID_JOB_INFO(in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, pic_in != NULL ? 1 : 0);

  //Hor Scaling
  if (pic_in != NULL) {

    int range[4] = { 0, 0, 0, 0 };
    if (in_param->use_tiles) {
      kvz_blockScalingSrcWidthRange(range, param, in_param->block_x, in_param->block_width);
    } else {
      //Get range that needs to be copied from pic_in
      //range[0:1] is prev blocks range and range[2:3] is the new block range
      //Only copy pixels not in range[0:1] (already copied by previous workers).
      //Getting the correct block_x of the previous block does not matter as only range[1] is used
      if (in_param->block_x - in_param->block_width < 0) {
        kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
        range[1] = range[2] - 1;
      } else {
        kvz_blockScalingSrcWidthRange(range, param, in_param->block_x - in_param->block_width, in_param->block_width);
        kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
      }
    }

    //int cp_block_x = in_param->use_tiles ? range[0] : (range[1] + 1);
    //int cp_block_y = in_param->block_y;
    //int cp_block_width = (in_param->use_tiles ? range[1] : range[3]) - cp_block_x + 1;
    //int cp_block_height = in_param->block_height;

    int hor_block_y = in_param->block_y;
    int hor_block_height = in_param->block_height;

    if (pic_out != NULL) {
      if (in_param->use_tiles) {
        kvz_blockScalingSrcHeightRange(range, param, in_param->block_y, in_param->block_height);
        //cp_block_y = hor_block_y = range[0];
        //cp_block_height = hor_block_height = range[1] - hor_block_y + 1;
      } else {
        //Do the same procedure as with horizontal range
        if (in_param->block_y - in_param->block_height < 0) {
          kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
          range[1] = range[2] - 1;
        } else {
          kvz_blockScalingSrcHeightRange(range, param, in_param->block_y - in_param->block_height, in_param->block_height);
          kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
        }
        //cp_block_y = range[1] + 1;
        //cp_block_height = range[3] - cp_block_y + 1;
        hor_block_y = range[2];
        hor_block_height = range[3] - range[2] + 1;
      }
    }

    //When using tiles, copy from in_pic to the src buffer (src buffer should hold only one tile and start from indexing (0,0))
    //int cp_dst_x = in_param->use_tiles ? 0 : cp_block_x;
    //int cp_dst_y = in_param->use_tiles ? 0 : cp_block_y;

    //Copy from in_pic to the src buffer
    //kvz_copy_uint8_block_to_YuvBuffer(in_param->src_buffer, pic_in->y, pic_in->u, pic_in->v, pic_in->stride, cp_dst_x, cp_dst_y, cp_block_x, cp_block_y, cp_block_width, cp_block_height, w_factor, h_factor);

    //PRINT_JOB_EXTRA_INFO("Copy to src buffer", cp_dst_x, cp_dst_y, cp_block_x, cp_block_y, cp_block_width, cp_block_height);
    PRINT_JOB_EXTRA_INFO("Hor scaling", in_param->block_x, hor_block_y, 0, 0, in_param->block_width, hor_block_height);

    //If both ver and hor done at the same time interpred in_param->block_y/height as the final output block and so we need to do hor scaling in the approriate range to accomodate the final block
    if (!kvz_opaqueYuvBlockStepScaling_adapter(in_param->ver_tmp_buffer, in_param->src_buffer, param, in_param->block_x, hor_block_y, in_param->block_width, hor_block_height, 0, kvz_opaque_resample_block_step)) {
      //TODO: Do error stuff?
      kvz_image_free(pic_in);
      kvz_image_free(pic_out);
      free(in_param);
      return;
    }
  }

  //Ver scaling
  if (pic_out != NULL) {
    //Do ver scaling step
    if (!kvz_opaqueYuvBlockStepScaling_adapter(in_param->trgt_buffer, in_param->ver_tmp_buffer, param, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, 1, kvz_opaque_resample_block_step)) {
      //TODO: Do error stuff?
      kvz_image_free(pic_in);
      kvz_image_free(pic_out);
      free(in_param);
      return;
    }

    /*const int dst_x = in_param->block_x;
    const int dst_y = in_param->block_y;
    const int src_x = in_param->use_tiles ? 0 : in_param->block_x;
    const int src_y = in_param->use_tiles ? 0 : in_param->block_y;*/

    //Copy results to pic_out
    //kvz_copy_YuvBuffer_block_to_uint8(pic_out->y, pic_out->u, pic_out->v, pic_out->stride, in_param->trgt_buffer, dst_x, dst_y, src_x, src_y, in_param->block_width, in_param->block_height, w_factor, h_factor);

  }

  //Do deallocation
  kvz_image_free(pic_in);
  kvz_image_free(pic_out);
  free(in_param);
}

//Handle hor/ver scaling steps:
//  If pic_in is given, copy relevant block to src_buffer and run horizontal scaling step
//  If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out
//  If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
//void kvz_block_step_scaler_worker(void * opaque_param)
//{
//  kvz_image_scaling_parameter_t *in_param = opaque_param;
//  kvz_picture * const pic_in = in_param->pic_in;
//  kvz_picture * const pic_out = in_param->pic_out;
//  const scaling_parameter_t *const param = in_param->param;
//
//  //TODO: account for chroma format properly
//  int w_factor = -1;
//  int h_factor = -1;
//
//  //Hor Scaling
//  if( pic_in != NULL ){
//    //Get range that needs to be copied from pic_in
//    //range[0:1] is prev blocks range and range[2:3] is the new block range
//    //Only copy pixels not in range[0:1] (already copied by previous workers).
//    //Getting the correct block_x of the previous block does not matter as only range[1] is used
//    int range[4];
//    if (in_param->block_x - in_param->block_width < 0){
//      kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
//      range[1] = range[2] - 1;
//    } else {
//      kvz_blockScalingSrcWidthRange(range, param, in_param->block_x - in_param->block_width, in_param->block_width);
//      kvz_blockScalingSrcWidthRange(range + 2, param, in_param->block_x, in_param->block_width);
//    }
//
//    int cp_block_x = range[1] + 1;
//    int cp_block_y = in_param->block_y;
//    int cp_block_width = range[3] - cp_block_x + 1;
//    int cp_block_height = in_param->block_height;
//    int hor_block_y = in_param->block_y;
//    int hor_block_height = in_param->block_height;
//
//    if(pic_out != NULL ){
//      if (in_param->block_y - in_param->block_height < 0) {
//        kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
//        range[1] = range[2] - 1;
//      }
//      else {
//        kvz_blockScalingSrcHeightRange(range, param, in_param->block_y - in_param->block_height, in_param->block_height);
//        kvz_blockScalingSrcHeightRange(range + 2, param, in_param->block_y, in_param->block_height);
//      }
//      cp_block_y = range[1] + 1;
//      cp_block_height = range[3] - cp_block_y + 1;
//      hor_block_y = range[2];
//      hor_block_height = range[3] - range[2] + 1;
//    }
//
//    //Copy from in_pic to the src buffer
//    kvz_copy_uint8_block_to_YuvBuffer(in_param->src_buffer, pic_in->y, pic_in->u, pic_in->v, pic_in->stride, cp_block_x, cp_block_y, cp_block_x, cp_block_y, cp_block_width, cp_block_height, w_factor, h_factor);
//
//    //If both ver and hor done at the same time interpred in_param->block_y/height as the final output block and so we need to do hor scaling in the approriate range to accomodate the final block
//    if (!kvz_yuvBlockStepScaling(in_param->ver_tmp_buffer, in_param->src_buffer, param, in_param->block_x, hor_block_y, in_param->block_width, hor_block_height, 0)) {  
//      //TODO: Do error stuff?
//      kvz_image_free(pic_in);
//      kvz_image_free(pic_out);
//      free(in_param);
//      return;
//    }
//  }
//
//  //Ver scaling
//  if (pic_out != NULL) {
//    //Do ver scaling step
//    if (!kvz_yuvBlockStepScaling(in_param->trgt_buffer, in_param->ver_tmp_buffer, param, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, 1)) {
//      //TODO: Do error stuff?
//      kvz_image_free(pic_in);
//      kvz_image_free(pic_out);
//      free(in_param);
//      return;
//    }
//
//    //Copy results to pic_out
//    kvz_copy_YuvBuffer_block_to_uint8(pic_out->y, pic_out->u, pic_out->v, pic_out->stride, in_param->trgt_buffer, in_param->block_x, in_param->block_y, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, w_factor, h_factor);
//    
//  }
//
//  //Do deallocation
//  kvz_image_free(pic_in);
//  kvz_image_free(pic_out);
//  free(in_param);
//
//}

//Handle hor/ver scaling steps for tiles:
//  If pic_in is given, copy relevant tile to src_buffer and run horizontal scaling step; src_buffer is filled starting from (0,0)
//  If pic_out is given, run vertical scaling step and copy relevant block from trgt_buffer to pic_out; trgt_buffer is indexed starting from (0,0)
//  If both are given, do both directions and the given block is taken to mean the pic_out block that should be calculated
//void kvz_tile_step_scaler_worker(void * opaque_param)
//{
//  kvz_image_scaling_parameter_t *in_param = opaque_param;
//  kvz_picture * const pic_in = in_param->pic_in;
//  kvz_picture * const pic_out = in_param->pic_out;
//  const scaling_parameter_t *const param = in_param->param;
//
//  //TODO: account for chroma format properly
//  int w_factor = -1;
//  int h_factor = -1;
//
//  //Hor Scaling
//  if (pic_in != NULL) {
//    //Get range needed to be copied from pic_in
//    int range[2];
//    kvz_blockScalingSrcWidthRange(range, param, in_param->block_x, in_param->block_width);
//    
//    int cp_block_x = range[0];
//    int hor_block_y = in_param->block_y;
//    int cp_block_width = range[1] - cp_block_x + 1;
//    int hor_block_height = in_param->block_height;
//
//    if (pic_out != NULL) {
//      kvz_blockScalingSrcHeightRange(range, param, in_param->block_y, in_param->block_height);
//      hor_block_y = range[0];
//      hor_block_height = range[1] - hor_block_y + 1;
//    }
//
//    //Copy from in_pic to the src buffer (src buffer should hold only one tile and start from indexing (0,0))
//    kvz_copy_uint8_block_to_YuvBuffer(in_param->src_buffer, pic_in->y, pic_in->u, pic_in->v, pic_in->stride, 0, 0, cp_block_x, hor_block_y, cp_block_width, hor_block_height, w_factor, h_factor);
//
//    //If both ver and hor done at the same time interpred in_param->block_y/height as the final output block and so we need to do hor scaling in the approriate range to accomodate the final block
//    if (!kvz_yuvBlockStepScaling(in_param->ver_tmp_buffer, in_param->src_buffer, param, in_param->block_x, hor_block_y, in_param->block_width, hor_block_height, 0)) {
//      //TODO: Do error stuff?
//      kvz_image_free(pic_in);
//      kvz_image_free(pic_out);
//      free(in_param);
//      return;
//    }
//  }
//
//  //Ver scaling
//  if (pic_out != NULL) {
//    //Do ver scaling step
//    if (!kvz_yuvBlockStepScaling(in_param->trgt_buffer, in_param->ver_tmp_buffer, param, in_param->block_x, in_param->block_y, in_param->block_width, in_param->block_height, 1)) {
//      //TODO: Do error stuff?
//      kvz_image_free(pic_in);
//      kvz_image_free(pic_out);
//      free(in_param);
//      return;
//    }
//
//    //Copy results to pic_out (trgt_buffer should only hold one tile and start from (0,0)
//    kvz_copy_YuvBuffer_block_to_uint8(pic_out->y, pic_out->u, pic_out->v, pic_out->stride, in_param->trgt_buffer, in_param->block_x, in_param->block_y, 0, 0, in_param->block_width, in_param->block_height, w_factor, h_factor);
//
//  }
//
//  //Do deallocation
//  kvz_image_free(pic_in);
//  kvz_image_free(pic_out);
//  free(in_param);
//
//}

//TODO: Reuse buffers? Or not, who cares. Use a scaler struct to hold all relevant info for different layers?
//TODO: remove memory db stuff
//Create a new kvz picture based on pic_in with size given by width and height
kvz_picture* kvz_image_scaling(kvz_picture* const pic_in, const scaling_parameter_t *const param, uint8_t skip_same)
{

  if(pic_in == NULL) {
    return NULL;
  }

  //If no scaling needs to be done, just return pic_in
  if (skip_same && param->src_height == param->trgt_height && param->src_width == param->trgt_width) {
    return kvz_image_copy_ref(pic_in);
  }

  kvz_image_scaling_parameter_t *scaling_param = calloc(1, sizeof(kvz_image_scaling_parameter_t));

  kvz_picture* pic_out = kvz_image_alloc(pic_in->chroma_format,
                                          param->trgt_width + param->trgt_padding_x,
                                          param->trgt_height + param->trgt_padding_y);

  
  //Allocate scaling parameters to give to the worker. Worker should handle freeing.
  scaling_param->pic_in = kvz_image_copy_ref(pic_in);
  scaling_param->pic_out = kvz_image_copy_ref(pic_out);
  scaling_param->param = param;

  //Do scaling
  kvz_image_scaler_worker(scaling_param);

  //Pic out should now contain the scaled image
  return pic_out;
}

void kvz_propagate_image_scaling_parameters(kvz_image_scaling_parameter_t * const dst, const kvz_image_scaling_parameter_t * const src)
{
  //Make sure that if dst and src pics are the same, the pic is not deallocated.
  kvz_picture* tmp_in = dst->pic_in;
  kvz_picture* tmp_out = dst->pic_out;

  dst->pic_in = src->pic_in != NULL ? kvz_image_copy_ref(src->pic_in) : NULL;
  dst->pic_out = src->pic_out != NULL ? kvz_image_copy_ref(src->pic_out) : NULL;

  kvz_image_free(tmp_in);
  kvz_image_free(tmp_out);

  //propagate the opaque buffers
  if (dst->src_buffer != NULL) {
    kvz_setOpaqueYuvBuffer(dst->src_buffer, src->src_buffer->y->data, src->src_buffer->u->data, src->src_buffer->v->data, src->src_buffer->y->depth);
  } else {
    dst->src_buffer = src->src_buffer;
  }
  if (dst->trgt_buffer != NULL) {
    kvz_setOpaqueYuvBuffer(dst->trgt_buffer, src->trgt_buffer->y->data, src->trgt_buffer->u->data, src->trgt_buffer->v->data, src->trgt_buffer->y->depth);
  } else {
    dst->trgt_buffer = src->trgt_buffer;
  }

  if (dst->ver_tmp_buffer == NULL) dst->ver_tmp_buffer = src->ver_tmp_buffer;

  dst->block_x = src->block_x;
  dst->block_y = src->block_y;
  dst->block_width = src->block_width;
  dst->block_height = src->block_height;
  dst->param = src->param;

  dst->use_tiles = src->use_tiles;
}

void kvz_copy_image_scaling_parameters(kvz_image_scaling_parameter_t * const dst, const kvz_image_scaling_parameter_t * const src)
{
  //Make sure that if dst and src pics are the same, the pic is not deallocated.
  kvz_picture* tmp_in = dst->pic_in;
  kvz_picture* tmp_out = dst->pic_out;

  dst->pic_in = src->pic_in != NULL ? kvz_image_copy_ref(src->pic_in) : NULL;
  dst->pic_out = src->pic_out != NULL ? kvz_image_copy_ref(src->pic_out) : NULL;

  kvz_image_free(tmp_in);
  kvz_image_free(tmp_out);

  //Make sure buffers are not deallocated if they are the same
  if (dst->src_buffer != src->src_buffer) kvz_deallocateOpaqueYuvBuffer(dst->src_buffer, 0);
  if (dst->ver_tmp_buffer != src->ver_tmp_buffer) kvz_deallocateOpaqueYuvBuffer(dst->ver_tmp_buffer, 1);
  if (dst->trgt_buffer != src->trgt_buffer) kvz_deallocateOpaqueYuvBuffer(dst->trgt_buffer, 0);
  dst->src_buffer = src->src_buffer;
  dst->ver_tmp_buffer = src->ver_tmp_buffer;
  dst->trgt_buffer = src->trgt_buffer;

  dst->block_x = src->block_x;
  dst->block_y = src->block_y;
  dst->block_width = src->block_width;
  dst->block_height = src->block_height;
  dst->param = src->param;

  dst->use_tiles = src->use_tiles;
}
// ***********************************************
