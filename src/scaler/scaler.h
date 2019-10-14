#ifndef SCALER_H_
#define SCALER_H_
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

#include <stdint.h>


/*=====================Scaling parameter definition=====================*/
//Format for specifying the ratio between chroma and luma
//Match the ones used in kvazaar.h
typedef enum
{
  CHROMA_400 = 0,
  CHROMA_420 = 1,
  CHROMA_422 = 2,
  CHROMA_444 = 3
} chroma_format_t;

//TODO: Move to .c?
//TODO: Add offsets/cropping
typedef struct
{
  //Original parameters
  int src_width;
  int src_height;

  int trgt_width;
  int trgt_height;

  chroma_format_t chroma;

  //Resampling parameters
  int rnd_trgt_width;
  int rnd_trgt_height;

  int scaled_src_width;
  int scaled_src_height;

  //Sample positional parameters
  int right_offset;
  int bottom_offset;

  int shift_x;
  int shift_y;

  int scale_x;
  int scale_y;

  int add_x;
  int add_y;

  int delta_x;
  int delta_y;

  //Phase parameters mostly for chroma
  //int phase_x;
  //int phase_y;

  //Padding for when the src image size differes from the parameters
  int src_padding_x;
  int src_padding_y;
  int trgt_padding_x;
  int trgt_padding_y;

  int is_calculated; //Flag that tells that parameters have been calculated. Needs to be set to false manually if values are changed

} scaling_parameter_t;

/*==========================================================================*/

#define SCALER_BUFFER_PADDING 32 //Define padding added to picture buffers. Keeps avx2 optimizations from accessing bad memory.

/*===========================Scaling parameter utility functions=================================*/
/**
* \brief Returns the appropriate chroma format for the given parameters
*/
chroma_format_t kvz_getChromaFormat(int luma_width, int luma_height, int chroma_width, int chroma_height);

/**
* \brief Function for getting initial scaling parameters given src and trgt size parameters.
*/
//TODO: Get rid of is_upsampling (it just toggles rounded width stuff that propably should not be used)
scaling_parameter_t kvz_newScalingParameters(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma, int is_upsampling);
/**
* \brief Experimental. Function for getting initial scaling parameters given src and trgt size parameters.
*/
scaling_parameter_t kvz_newScalingParameters_(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma);
/*=============================================================================================*/


/*===============Buffer datastructure definitions=============*/
typedef int pic_data_t; //Use some other type?

/**
 * \brief Picture buffer type for operating on image data.
 */
typedef struct
{
  int width;
  int height;
  pic_data_t* data; //Contain main data
  
  pic_data_t* tmp_row; //A temporary buffer row that may be used to hold data when operating on buffer
  
} pic_buffer_t;

/**
* \brief Picture buffer type for yuv frames.
*/
typedef struct
{
  pic_buffer_t* y;
  pic_buffer_t* u;
  pic_buffer_t* v;

  chroma_format_t format;
} yuv_buffer_t;

/**
* \brief Opaque version for passing arbitrary buffers
* Should mirror pic_buffer so that opaque_pic_buffer_t can be used like pic_buffer_t
*/
typedef struct
{
  int width;
  int height;

  void *data;

  int stride;
  unsigned depth; //Determine the bit-depth of data
} opaque_pic_buffer_t;

/**
* \brief Opaque version for passing arbitrary buffers
* Should mirror yuv_buffer_t so that opaque_yuv_buffer_t can be cast to yuv_buffer_t
*/
typedef struct
{
  opaque_pic_buffer_t* y;
  opaque_pic_buffer_t* u;
  opaque_pic_buffer_t* v;

  chroma_format_t format;
} opaque_yuv_buffer_t;

/*==========================================================*/
/*==================Buffer utility functions===============*/
/**
 * \brief Create a Picture buffer. The caller is responsible for deallocation
 */
pic_buffer_t* kvz_newPictureBuffer(int width, int height, int has_tmp_row);
yuv_buffer_t* kvz_newYuvBuffer(int width, int height , chroma_format_t format, int has_tmp_row);

/**
* \brief return a opaque buffer type using the given input
* \param ?_data: set data buffer to given opaque buffer. If data is NULL, allocate memory based on alloc_depth.
* \param alloc_depth: Set element size and allocate needed memory. If depth is 0 does not allocate memory for data
*/
opaque_yuv_buffer_t* kvz_newOpaqueYuvBuffer(void *const y_data, void *const u_data, void *const v_data, int width, int height, int stride, chroma_format_t format, const unsigned alloc_depth);
opaque_pic_buffer_t* kvz_newOpaquePictureBuffer(void *const data, int width, int height, int stride, const unsigned alloc_depth);

/**
* \brief set the opaque data buffers and depth value for an existing buffer
*/
void kvz_setOpaqueYuvBuffer(opaque_yuv_buffer_t *buffer, void *const y_data, void *const u_data, void *const v_data, const unsigned depth);
void kvz_setOpaquePicBuffer(opaque_pic_buffer_t *buffer, void *const data, const unsigned depth);

/**
* \brief Create/Initialize a yuv buffer. Width/height should be the width/height of the data. The caller is responsible for deallocation
*/
//yuv_buffer_t* newYuvBuffer_double(const double* const y_data, const double* const u_data, const double* const v_data, int width, int height, chroma_format_t format, int has_tmp_row);
//yuv_buffer_t* newYuvBuffer_uint8(const uint8_t* const y_data, const uint8_t* const u_data, const uint8_t* const v_data, int width, int height, chroma_format_t format, int has_tmp_row);

/**
* \brief Create/Initialize a yuv buffer. Width/height should be the width/height of the final buffer. Stride should be the width of the input (padded image). The caller is responsible for deallocation
*/
yuv_buffer_t* kvz_newYuvBuffer_padded_uint8(const uint8_t* const y_data, const uint8_t* const u_data, const uint8_t* const v_data, int width, int height, int stride, chroma_format_t format, int has_tmp_row);


/**
* \brief Clone the given yuv buffer (deep copy)
*/
yuv_buffer_t* kvz_cloneYuvBuffer(const yuv_buffer_t* const yuv);

/**
* \brief Copy the given yuv buffer (shallow copy)
*/
opaque_yuv_buffer_t* kvz_copyOpaqueYuvBuffer(const opaque_yuv_buffer_t* const yuv);


/**
* \brief Deallocate picture buffer
*/
void kvz_deallocatePictureBuffer(pic_buffer_t* buffer);

/**
* \brief Deallocate yuv buffer
*/
void kvz_deallocateYuvBuffer(yuv_buffer_t* yuv);

/**
* \brief Deallocate opaque picture buffer
*/
void kvz_deallocateOpaquePictureBuffer(opaque_pic_buffer_t* buffer, const int free_buffer);

/**
* \brief Deallocate opaque yuv buffer
*/
void kvz_deallocateOpaqueYuvBuffer(opaque_yuv_buffer_t* yuv, const int free_buffer);


/**
* \brief Copy a block from a uint8 src to a yuv buffer.
*/
void kvz_copy_uint8_block_to_YuvBuffer(const yuv_buffer_t* dst, const uint8_t* const y, const uint8_t* const u, const uint8_t* const v, const int luma_stride, const int dst_x, const int dst_y, const int src_x, const int src_y, const int block_width, const int block_height, const int w_factor, const int h_factor);
/**
* \brief Copy a block from a yuv buffer to a uint8 dst.
*/
void kvz_copy_YuvBuffer_block_to_uint8(uint8_t* const y, uint8_t* const u, uint8_t* const v, const int luma_stride, const yuv_buffer_t * const src, const int dst_x, const int dst_y, const int src_x, const int src_y, const int block_width, const int block_height, const int w_factor, const int h_factor);
/*=======================================================*/

/*============================Resample function pointer typedef================================*/
typedef void (resample_block_step_func)(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical);

typedef void (resample_func)(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma);

typedef void (opaque_resample_block_step_func)(const opaque_pic_buffer_t* const src_buffer, const opaque_pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical);

extern resample_block_step_func *const kvz_default_block_step_resample_func;
extern resample_func *const kvz_default_resample_func;
extern resample_func *const kvz_alt_resample_func;
extern opaque_resample_block_step_func *const kvz_opaque_block_step_resample_func;
/*=============================================================================================*/

/*================Main scaling functions========================*/
//TODO: Return/recycle the same buffer for the scaled yuv. Use yuv it self and not a separate buffer?
/**
* \brief Function for scaling an image given in a yuv buffer (can handle down- and upscaling).
*        Returns result in yuv buffer. If dst is null or incorrect size, allocate new buffer and return it (dst is deallocated). If dst is a usable buffer, returns the given dst
*/
yuv_buffer_t* kvz_yuvScaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst);

/**
* \brief Function for scaling an image given in a yuv buffer (can handle down- and upscaling).
* \param resample_func Pass the function to be used for resampling
*        Returns result in yuv buffer. If dst is null or incorrect size, allocate new buffer and return it (dst is deallocated). If dst is a usable buffer, returns the given dst
*/
yuv_buffer_t* kvz_yuvScaling_adapter(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst, resample_func *const resample_func);

/**
* \brief Experimental. Function for scaling an image given in a yuv buffer (can handle down- and upscaling).
*        Returns result in yuv buffer. If dst is null or incorrect size, allocate new buffer and return it (dst is deallocated). If dst is a usable buffer, returns the given dst
* \pre yuv and dst must have tmp rows that are either NULL or valid and guaranteed to be atleast SCALER_MAX(width,height) of the respective pic buffer.
* \post the larger of yuv and dst will have valid tmp rows in it's pic buffers. 
*/
yuv_buffer_t* kvz_yuvScaling_(yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst);

/**
* \brief Function for scaling an image given in a yuv buffer (can handle down- and upscaling).
* \pre dst should be a buffer of either size block_width-by-block_height or the size of the trgt image.
*        Result given in dst buffer.
*/
int kvz_yuvBlockScaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst, const int block_x, const int block_y, const int block_width, const int block_height);

/**
* \brief Function for scaling an image, given in a yuv buffer, in either the vertical or horizontal direction (can handle down- and upscaling).
* \pre dst should be a buffer of either size block_width-by-block_height or the size of the trgt image. And src should be large enough to accomodate the block schaling src range
*        Result given in dst buffer.
*/
int kvz_yuvBlockStepScaling( yuv_buffer_t* const dst, const yuv_buffer_t* const src, const scaling_parameter_t* const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical);

/**
* \brief Function for scaling an image, given in a yuv buffer, in either the vertical or horizontal direction (can handle down- and upscaling).
* \pre dst should be a buffer of either size block_width-by-block_height or the size of the trgt image. And src should be large enough to accomodate the block schaling src range
*        Result given in dst buffer.
* \param resample_func uses the given function pointer as the resample function
*/
int kvz_yuvBlockStepScaling_adapter(yuv_buffer_t* const dst, const yuv_buffer_t* const src, const scaling_parameter_t* const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical, resample_block_step_func *const resample_func);

/**
* \brief Function for scaling an image, given in a opaque yuv buffer, in either the vertical or horizontal direction (can handle down- and upscaling).
* \pre dst should be a buffer of either size block_width-by-block_height or the size of the trgt image. And src should be large enough to accomodate the block schaling src range
*        Result given in dst buffer.
* \param resample_func uses the given function pointer as the resample function
*/
int kvz_opaqueYuvBlockStepScaling_adapter(opaque_yuv_buffer_t * const dst, const opaque_yuv_buffer_t * const src, const scaling_parameter_t * const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical, opaque_resample_block_step_func * const resample_func);

/**
* \brief Function for scaling an image, given in a opaque yuv buffer, in either the vertical or horizontal direction (can handle down- and upscaling).
* \pre dst should be a buffer of either size block_width-by-block_height or the size of the trgt image. And src should be large enough to accomodate the block schaling src range
*        Result given in dst buffer.
*/
int kvz_opaqueYuvBlockStepScaling(opaque_yuv_buffer_t * const dst, const opaque_yuv_buffer_t * const src, const scaling_parameter_t * const base_param, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical);
/*=============================================================*/

/*================Block scaling helper functions========================*/
/**
* \brief Function for calculating the src pixel range used for block scaling when using the given parameters.
*        Result put in the range array.
*/
void kvz_blockScalingSrcWidthRange(int range[2], const scaling_parameter_t* const base_param, const int block_x, const int block_width);

/**
* \brief Function for calculating the src pixel range used for block scaling when using the given parameters.
*        Result put in the range array.
*/
void kvz_blockScalingSrcHeightRange(int range[2], const scaling_parameter_t* const base_param, const int block_y, const int block_height);
/*======================================================================*/

#endif
