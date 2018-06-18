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

} scaling_parameter_t;

/*==========================================================================*/
/*===========================Scaling parameter utility functions=================================*/
/**
* \brief Returns the appropriate chroma format for the given parameters
*/
chroma_format_t kvz_getChromaFormat(int luma_width, int luma_height, int chroma_width, int chroma_height);

/**
* \brief Function for getting initial scaling parameters given src and trgt size parameters.
*/
scaling_parameter_t kvz_newScalingParameters(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma);
/**
* \brief Experimental. Function for getting initial scaling parameters given src and trgt size parameters.
*/
scaling_parameter_t kvz_newScalingParameters_(int src_width, int src_height, int trgt_width, int trgt_height, chroma_format_t chroma);
/*=============================================================================================*/


/*===============Buffer datastructure definitions=============*/
typedef unsigned int pic_data_t; //Use some other type?

/**
 * \brief Picture buffer type for operating on image data.
 */
typedef struct
{
  pic_data_t* data; //Contain main data
  pic_data_t* tmp_row; //A temporary buffer row that may be used to hold data when operating on buffer

  int width;
  int height;
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

/*==========================================================*/
/*==================Buffer utility functions===============*/
/**
 * \brief Create a Picture buffer. The caller is responsible for deallocation
 */
pic_buffer_t* kvz_newPictureBuffer(int width, int height, int has_tmp_row);
yuv_buffer_t* kvz_newYuvBuffer(int width, int height , chroma_format_t format, int has_tmp_row);


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
* \brief Clone the given yuv buffer
*/
yuv_buffer_t* kvz_cloneYuvBuffer(const yuv_buffer_t* const yuv);

/**
* \brief Create/Initialize a Picture buffer. The caller is responsible for deallocation
*/
void kvz_deallocateYuvBuffer(yuv_buffer_t* yuv);

/*=======================================================*/


/*================Main scaling functions========================*/
//TODO: Return/recycle the same buffer for the scaled yuv. Use yuv it self and not a separate buffer?
/**
* \brief Function for scaling an image given in a yuv buffer (can handle down- and upscaling).
*        Returns result in yuv buffer. If dst is null or incorrect size, allocate new buffer and return it (dst is deallocated). If dst is a usable buffer, returns the given dst
*/
yuv_buffer_t* kvz_yuvScaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, yuv_buffer_t* dst);
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
int kvz_yuvBlockStepScaling(const yuv_buffer_t* const src, const scaling_parameter_t* const base_param, yuv_buffer_t* dst, const int block_x, const int block_y, const int block_width, const int block_height, const int is_vertical);

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
