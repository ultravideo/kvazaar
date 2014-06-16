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

#include "threads.h"
#include "image.h"
#include "strategyselector.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "checkpoint.h"
#include "sao.h"

/**
 * \brief Allocate new image
 * \return image pointer
 */
image *image_alloc(const int32_t width, const int32_t height, const int32_t poc)
{
  image *im = MALLOC(image, 1);
  
  unsigned int luma_size = width * height;
  unsigned int chroma_size = luma_size / 4;
  
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  if (!im) return NULL;
  
  im->width = width;
  im->height = height;
  im->stride = width;
  
  im->base_image = im;
  
  im->refcount = 1; //We give a reference to caller
  
  im->poc = poc;
  
  //Allocate memory
  im->fulldata = MALLOC(pixel, (luma_size + 2*chroma_size));
  im->y = im->data[COLOR_Y] = &im->fulldata[0];
  im->u = im->data[COLOR_U] = &im->fulldata[luma_size];
  im->v = im->data[COLOR_V] = &im->fulldata[luma_size + chroma_size];

  return im;
}

/**
 * \brief Free memory allocated to picture (if we have no reference left)
 * \param pic picture pointer
 * \return 1 on success, 0 on failure
 */
int image_free(image * const im)
{
  //Either we are the base image, or we should have no references
  assert(im->base_image == im || im->refcount == 0);
  
  int32_t new_refcount = ATOMIC_DEC(&(im->base_image->refcount));
  //If we're freeing a subimage, then we must free the pointer
  //Base image may be stored in image_list, and should not be freed
  //FIXME I don't find this very clean...
  if (new_refcount > 0 && im->base_image != im) free(im);
  if (new_refcount > 0) return 1;
  FREE_POINTER(im->base_image->fulldata);
  
  //Just to make the program crash when using those values after the free
  im->y = im->u = im->v = im->data[COLOR_Y] = im->data[COLOR_U] = im->data[COLOR_V] = NULL;
  
  free(im);

  return 1;
}


image *image_make_subimage(image * const orig_image, const unsigned int x_offset, const unsigned int y_offset, const unsigned int width, const unsigned int height)
{
  image *im = MALLOC(image, 1);
  if (!im) return NULL;
  
  im->base_image = orig_image->base_image;
  ATOMIC_INC(&(im->base_image->refcount));
  
  assert(x_offset + width <= orig_image->width);
  assert(y_offset + height <= orig_image->height);
  
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);
  
  assert((x_offset % 2) == 0);
  assert((y_offset % 2) == 0);
  
  im->stride = orig_image->stride;
  im->refcount = 0; //No references on subimages
  
  im->width = width;
  im->height = height;
  
  im->y = im->data[COLOR_Y] = &orig_image->y[x_offset + y_offset * orig_image->stride];
  im->u = im->data[COLOR_U] = &orig_image->u[x_offset/2 + y_offset/2 * orig_image->stride/2];
  im->v = im->data[COLOR_V] = &orig_image->v[x_offset/2 + y_offset/2 * orig_image->stride/2];

  return im;
}

yuv_t * yuv_t_alloc(int luma_size)
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

void yuv_t_free(yuv_t * yuv)
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
static unsigned cor_sad(const pixel *pic_data, const pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride)
{
  pixel ref = *ref_data;
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
static unsigned ver_sad(const pixel *pic_data, const pixel *ref_data,
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
static unsigned hor_sad(const pixel *pic_data, const pixel *ref_data,
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
static unsigned image_interpolated_sad(const image *pic, const image *ref,
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
    result += reg_sad(&pic_data[top * pic->stride + left],
                      &ref_data[top * ref->stride + left],
                      block_width - left, block_height - top, pic->stride, ref->stride);
  } else if (top && right) {
    result += ver_sad(pic_data,
                      &ref_data[top * ref->stride],
                      block_width - right, top, pic->stride);
    result += cor_sad(&pic_data[block_width - right],
                      &ref_data[top * ref->stride + (block_width - right - 1)],
                      right, top, pic->stride);
    result += reg_sad(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride],
                      block_width - right, block_height - top, pic->stride, ref->stride);
    result += hor_sad(&pic_data[top * pic->stride + (block_width - right)],
                      &ref_data[top * ref->stride + (block_width - right - 1)],
                      right, block_height - top, pic->stride, ref->stride);
  } else if (bottom && left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height - bottom, pic->stride, ref->stride);
    result += reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height - bottom, pic->stride, ref->stride);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      left, bottom, pic->stride);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->stride + left],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      block_width - left, bottom, pic->stride);
  } else if (bottom && right) {
    result += reg_sad(pic_data,
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
    result += reg_sad(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride],
                      block_width, block_height - top, pic->stride, ref->stride);
  } else if (bottom) {
    result += reg_sad(pic_data,
                      ref_data,
                      block_width, block_height - bottom, pic->stride, ref->stride);
    result += ver_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride],
                      block_width, bottom, pic->stride);
  } else if (left) {
    result += hor_sad(pic_data,
                      &ref_data[left],
                      left, block_height, pic->stride, ref->stride);
    result += reg_sad(&pic_data[left],
                      &ref_data[left],
                      block_width - left, block_height, pic->stride, ref->stride);
  } else if (right) {
    result += reg_sad(pic_data,
                      ref_data,
                      block_width - right, block_height, pic->stride, ref->stride);
    result += hor_sad(&pic_data[block_width - right],
                      &ref_data[block_width - right - 1],
                      right, block_height, pic->stride, ref->stride);
  } else {
    result += reg_sad(pic_data, ref_data, block_width, block_height, pic->stride, ref->stride);
  }

  return result;
}


unsigned image_calc_sad(const image *pic, const image *ref, int pic_x, int pic_y, int ref_x, int ref_y,
                        int block_width, int block_height) {
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);
  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];
    const pixel *ref_data = &ref->y[ref_y * ref->stride + ref_x];
    return reg_sad(pic_data, ref_data, block_width, block_height, pic->stride, ref->stride);
  } else {
    // Call a routine that knows how to interpolate pixels outside the frame.
    return image_interpolated_sad(pic, ref, pic_x, pic_y, ref_x, ref_y, block_width, block_height);
  }
}


unsigned pixels_calc_ssd(const pixel *const ref, const pixel *const rec,
                 const int ref_stride, const int rec_stride,
                 const int width)
{
  int ssd = 0;
  int y, x;

  for (y = 0; y < width; ++y) {
    for (x = 0; x < width; ++x) {
      int diff = ref[x + y * ref_stride] - rec[x + y * rec_stride];
      ssd += diff * diff;
    }
  }

  return ssd;
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
void pixels_blit(const pixel * const orig, pixel * const dst,
                         const unsigned width, const unsigned height,
                         const unsigned orig_stride, const unsigned dst_stride)
{
  unsigned y;
  //There is absolutely no reason to have a width greater than the source or the destination stride.
  assert(width <= orig_stride);
  assert(width <= dst_stride);

#ifdef CHECKPOINTS
  for (y = 0; y < height; ++y) {
    char buffer[3*width];
    int p;
    for (p = 0; p < width; ++p) {
      sprintf((buffer + 3*p), "%02X ", orig[y*orig_stride]);
    }
    buffer[3*width] = 0;
    CHECKPOINT("pixels_blit: %04d: %s", y, buffer);
  }
#endif //CHECKPOINTS

  if (orig == dst) {
    //If we have the same array, then we should have the same stride
    assert(orig_stride == dst_stride);
    return;
  }
  assert(orig != dst || orig_stride == dst_stride);

  for (y = 0; y < height; ++y) {
    memcpy(&dst[y*dst_stride], &orig[y*orig_stride], width * sizeof(pixel));
  }
}

