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

#include "scaler-avx2.h"
#include "scaler-util.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>

#define DEFAULT_RESAMPLE_BLOCK_STEP_FUNC_AVX2 resampleBlockStep_avx2
#define DEFAULT_RESAMPLE_FUNC_AVX2 resample_avx2

#define B11011000 0xD8 //0b11011000

// Clip sum of add_val to each epi32 of lane
static __m256i clip_add_avx2(const int add_val, __m256i lane, const int min, const int max)
{
 __m256i min_epi32 = _mm256_set1_epi32(min);
 __m256i max_epi32 = _mm256_set1_epi32(max);

 __m256i add_values_epi32 = _mm256_set1_epi32(add_val);
 add_values_epi32 = _mm256_add_epi32(add_values_epi32, lane);

 // Compare integers so value is in range 0 - (src-width-1)
 add_values_epi32 = _mm256_max_epi32(min_epi32, add_values_epi32);
 add_values_epi32 = _mm256_min_epi32(add_values_epi32, max_epi32);
 
 return add_values_epi32;
}


//Resampling is done here per buffer
//void _resample_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
//{
// //TODO: Add cropping etc.
//
// //Choose best filter to use when downsampling
// //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
// int ver_filter = 0;
// int hor_filter = 0;
//
//  int src_height = param->src_height;
//  int src_width = param->src_width;
//  int trgt_height = param->trgt_height;//param->rnd_trgt_height;
//  int trgt_width = param->trgt_width;//param->rnd_trgt_width
//
// if (!is_upscaling) {
//  int crop_width = src_width - param->right_offset; //- param->left_offset;
//  int crop_height = src_height - param->bottom_offset; //- param->top_offset;
//
//  if (4 * crop_height > 15 * trgt_height)
//   ver_filter = 7;
//  else if (7 * crop_height > 20 * trgt_height)
//   ver_filter = 6;
//  else if (2 * crop_height > 5 * trgt_height)
//   ver_filter = 5;
//  else if (1 * crop_height > 2 * trgt_height)
//   ver_filter = 4;
//  else if (3 * crop_height > 5 * trgt_height)
//   ver_filter = 3;
//  else if (4 * crop_height > 5 * trgt_height)
//   ver_filter = 2;
//  else if (19 * crop_height > 20 * trgt_height)
//   ver_filter = 1;
//
//  if (4 * crop_width > 15 * trgt_width)
//   hor_filter = 7;
//  else if (7 * crop_width > 20 * trgt_width)
//   hor_filter = 6;
//  else if (2 * crop_width > 5 * trgt_width)
//   hor_filter = 5;
//  else if (1 * crop_width > 2 * trgt_width)
//   hor_filter = 4;
//  else if (3 * crop_width > 5 * trgt_width)
//   hor_filter = 3;
//  else if (4 * crop_width > 5 * trgt_width)
//   hor_filter = 2;
//  else if (19 * crop_width > 20 * trgt_width)
//   hor_filter = 1;
// }
//
// int shift_x = param->shift_x - 4;
// int shift_y = param->shift_y - 4;
// __m256i pointer, temp_mem, temp_filter, decrese;
// __m256i adder = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
// __m256i upscaling_adder = _mm256_set_epi32(8, 9, 10, 11, src_width, src_width, src_width, src_width);
// __m256i order = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
// __m128i smallest_epi16;
// int* start;
// int min;
//
// pic_data_t* tmp_row = buffer->tmp_row;
//
// // Horizontal downsampling
// for (int i = 0; i < src_height; i++) {
//  pic_data_t* src_row = &buffer->data[i * buffer->width];
//
//  for (int j = 0; j < trgt_width; j++) {
//   //Calculate reference position in src pic
//   int ref_pos_16 = (int)((unsigned int)(j * param->scale_x + param->add_x) >> shift_x) - param->delta_x;
//   int phase = ref_pos_16 & 15;
//   int ref_pos = ref_pos_16 >> 4;
//
//
//   //Choose filter
//   const int* filter;
//   int size = getFilter(&filter, is_upscaling, is_luma, phase, hor_filter);
//
//
//   pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, adder, 0, src_width - 1);
//   pointer = _mm256_permutevar8x32_epi32(pointer, order);
//
//   min = src_width - 1;
//   smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
//   smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
//   min = _mm_extract_epi16(smallest_epi16, 0);
//
//   tmp_row[j] = 0;
//   temp_mem = _mm256_load_si256((__m256i*)&(src_row[min]));
//   decrese = _mm256_set1_epi32(min);
//
//   pointer = _mm256_sub_epi32(pointer, decrese);
//
//   temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
//
//   temp_filter = _mm256_load_si256((__m256i*)&(filter[0]));
//   temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
//
//   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//
//   switch (size)
//   {
//   case 4:
//
//    tmp_row[j] = _mm256_extract_epi32(temp_mem, 0);
//    break;
//
//   case 8:
//    tmp_row[j] = _mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
//    break;
//
//
//   default:
//    tmp_row[j] = _mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
//
//    pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, upscaling_adder, 0, src_width - 1);
//    pointer = _mm256_permutevar8x32_epi32(pointer, order);
//
//    smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
//    smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
//    min = _mm_extract_epi16(smallest_epi16, 0);
//
//    temp_filter = _mm256_load_si256((__m256i*)&(filter[8]));
//    temp_mem = _mm256_load_si256((__m256i*)&(src_row[min]));
//    decrese = _mm256_set1_epi32(min);
//    pointer = _mm256_sub_epi32(pointer, decrese);
//    temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
//
//    temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
//
//
//    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//
//    tmp_row[j] += _mm256_extract_epi32(temp_mem, 0);
//    break;
//   }
//  }
//  //Copy tmp row to data
//  memcpy(src_row, tmp_row, sizeof(pic_data_t) * trgt_width);
// }
//
// pic_data_t* tmp_col = tmp_row; //rename for clarity
//
// // Vertical downsampling
// __m256i multiplier_epi32 = _mm256_set1_epi32(buffer->width);
// for (int i = 0; i < trgt_width; i++) {
//  pic_data_t* src_col = &buffer->data[i];
//  for (int j = 0; j < trgt_height; j++) {
//   //Calculate ref pos
//   int ref_pos_16 = (int)((unsigned int)(j * param->scale_y + param->add_y) >> shift_y) - param->delta_y;
//   int phase = ref_pos_16 & 15;
//   int ref_pos = ref_pos_16 >> 4;
//
//   //Choose filter
//   const int* filter;
//   int size = getFilter(&filter, is_upscaling, is_luma, phase, ver_filter);
//   /*
//   //Apply filter
//   tmp_col[j] = 0;
//   for (int k = 0; k < size; k++) {
//   int m = clip(ref_pos + k - (size >> 1) + 1, 0, src_height - 1);
//   tmp_col[j] += filter[k] * src_col[m * buffer->width];
//   }*/
//   //-------------------------------------------------------
//
//
//   pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, adder, 0, src_height - 1);
//
//   pointer = _mm256_permutevar8x32_epi32(pointer, order);
//   pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
//   start = (int*)&(pointer);
//
//   tmp_col[j] = 0;
//   temp_mem = _mm256_set_epi32(src_col[start[7]], src_col[start[6]], src_col[start[5]], src_col[start[4]], src_col[start[3]], src_col[start[2]], src_col[start[1]], src_col[start[0]]);
//
//   temp_filter = _mm256_load_si256((__m256i*)&(filter[0]));
//
//   temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
//
//   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//
//
//   switch (size)
//   {
//   case 4:
//
//    tmp_col[j] = _mm256_extract_epi32(temp_mem, 0);
//    break;
//
//   case 8:
//    tmp_col[j] = _mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
//    break;
//
//
//   default:
//    tmp_col[j] = _mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
//
//    pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, upscaling_adder, 0, src_height - 1);
//    pointer = _mm256_permutevar8x32_epi32(pointer, order);
//    pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
//
//    start = (int*)&(pointer);
//
//    temp_filter = _mm256_load_si256((__m256i*)&(filter[8]));
//
//    temp_mem = _mm256_set_epi32(src_col[start[7]], src_col[start[6]], src_col[start[5]], src_col[start[4]], src_col[start[3]], src_col[start[2]], src_col[start[1]], src_col[start[0]]);
//    temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
//
//
//    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
//
//    tmp_col[j] += _mm256_extract_epi32(temp_mem, 0);
//
//
//    break;
//   }
//   //---------------------------------------
//
//   //TODO: Why? Filter coefs summ up to 128 applied 2x 128*128= 2^14
//   //TODO: Why?  Filter coefs summ up to 64 applied 2x 64*64= 2^12
//   //Scale values back down
//   tmp_col[j] = is_upscaling ? (tmp_col[j] + 2048) >> 12 : (tmp_col[j] + 8192) >> 14;
//  }
//
//  //Clip and move to buffer data
//  for (int n = 0; n < trgt_height; n++) {
//   src_col[n * buffer->width] = clip(tmp_col[n], 0, 255);
//  }
// }
//}
//

//Resampling is done here per buffer
static void resample_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
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
 __m256i pointer, temp_mem, temp_filter, decrese;
 __m256i adder = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 __m256i upscaling_adder = _mm256_set_epi32(8, 9, 10, 11, src_width, src_width, src_width, src_width);
 __m256i order = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 __m128i smallest_epi16;
 int* start;
 int min;

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


   pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, adder, 0, src_width - 1);
   pointer = _mm256_permutevar8x32_epi32(pointer, order);

   min = src_width-1;
   smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
   smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
   min = _mm_extract_epi16(smallest_epi16, 0);
   
   tmp_row[j] = 0;
   //TODO: Gives a read access error when using the original load in vs2015 
   temp_mem = _mm256_maskload_epi32(&src_row[min], _mm256_set1_epi32(0xf0000000));//_mm256_load_si256((__m256i*)(&(src_row[min])));
   
   decrese = _mm256_set1_epi32(min);

   pointer = _mm256_sub_epi32(pointer, decrese);

   temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);

   temp_filter = _mm256_load_si256((__m256i*)&(filter[0])); //TODO: What if filter is not 256 bits?
   temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
   

   
   
   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);

   switch (size)
   {
   case 4:
    
    tmp_row[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0);//_mm256_extract_epi32(temp_mem, 0);
    break;

   case 8:
    tmp_row[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
    break;


   default:
    tmp_row[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);

    pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, upscaling_adder, 0, src_width - 1);
    pointer = _mm256_permutevar8x32_epi32(pointer, order);

    smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
    smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
    min = _mm_extract_epi16(smallest_epi16, 0);

    temp_filter = _mm256_load_si256((__m256i* )&(filter[8]));
    temp_mem = _mm256_load_si256((__m256i*)&(src_row[min]));

    decrese = _mm256_set1_epi32(min);
    pointer = _mm256_sub_epi32(pointer, decrese);
    temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);

    temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);   
    


    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);

    tmp_row[j] += _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0); //_mm256_extract_epi32(temp_mem, 0);
    break;
   }
  }
  //Copy tmp row to data
  memcpy(src_row, tmp_row, sizeof(pic_data_t) * trgt_width);
 }

 pic_data_t* tmp_col = tmp_row; //rename for clarity

  // Vertical resampling
 __m256i multiplier_epi32 = _mm256_set1_epi32(buffer->width);
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
   
   pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, adder, 0, src_height - 1);

   pointer = _mm256_permutevar8x32_epi32(pointer, order);
   pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
   start = (int*)&(pointer);

   tmp_col[j] = 0;
   temp_mem = _mm256_set_epi32(src_col[start[7]], src_col[start[6]], src_col[start[5]], src_col[start[4]], src_col[start[3]], src_col[start[2]], src_col[start[1]], src_col[start[0]]);

   temp_filter = _mm256_load_si256((__m256i*)&(filter[0]));

   temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);

   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
   temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);


   switch (size)
   {
   case 4:

    tmp_col[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0);//_mm256_extract_epi32(temp_mem, 0);
    break;

   case 8:
    tmp_col[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
    break;


   default:
    tmp_col[j] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);

    pointer = clip_add_avx2(ref_pos - (size >> 1) + 1, upscaling_adder, 0, src_height - 1);
    pointer = _mm256_permutevar8x32_epi32(pointer, order);
    pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);

    start = (int*)&(pointer);

    temp_filter = _mm256_load_si256((__m256i*)&(filter[8]));

    temp_mem = _mm256_set_epi32(src_col[start[7]], src_col[start[6]], src_col[start[5]], src_col[start[4]], src_col[start[3]], src_col[start[2]], src_col[start[1]], src_col[start[0]]);
    temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);


    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
    temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);

    tmp_col[j] += _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0); //_mm256_extract_epi32(temp_mem, 0);
    
    
    break;
   }
   //---------------------------------------

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


static void resampleBlockStep_avx2(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical)
{
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
  const int s_stride = is_vertical ? src_buffer->width : 1; //Multiplier to s_ind

  //Specify bounds for trgt buffer and filter
  //const int trgt_bound = is_vertical ? inner_bound : outer_bound;
  const int filter_bound = is_vertical ? outer_bound : inner_bound;

  const int o_step = is_vertical ? filter_size : 1; //Step size of outer loop. Adjust depending on how many values can be calculated concurrently with SIMD 
  const int i_step = is_vertical ? 1 : filter_size; //Step size of inner loop. Adjust depending on how many values can be calculated concurrently with SIMD

  __m256i base_mask = _mm256_set_epi32(0x80000000, 0xc0000000, 0xe0000000, 0xf0000000, 0xf8000000, 0xfc000000, 0xfe000000, 0xff000000); //Mask used to load/store. Shifting to the left allows masking out values from the end if width is not a multiple of step size etc.
  __m256i pointer, temp_mem, temp_filter, decrese;
  __m256i adder = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  __m256i downsampling_adder = _mm256_set_epi32(8, 9, 10, 11, src_size, src_size, src_size, src_size);
  __m256i order = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  __m128i smallest_epi16;
  __m256i multiplier_epi32 = _mm256_set1_epi32(s_stride); //Stride of src buffer
  //int* start;
  int min = 0;

  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < (block_y + block_height); y++) {

    pic_data_t* src = is_vertical ? src_buffer->data : &src_buffer->data[y * src_buffer->width];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width];

    //Outer loop:
    //  if is_vertical -> loop over k (filter inds)
    //  if !is_vertical -> loop over x (block width)
    for (int o_ind = outer_init; o_ind < outer_bound; o_ind += o_step) {

      const int t_ind = is_vertical ? y : o_ind; //trgt_buffer row/col index for cur resampling dir

      //Calculate reference position in src pic
      const int ref_pos_16 = (int)((unsigned int)(t_ind * scale + add) >> shift) - delta;
      const int phase = ref_pos_16 & 15;
      const int ref_pos = ref_pos_16 >> 4;

      //Inner loop:
      //  if is_vertical -> loop over x (block width)
      //  if !is_vertical -> loop over k (filter inds)-
      for (int i_ind = inner_init; i_ind < inner_bound; i_ind += i_step) {

        const int f_ind = is_vertical ? o_ind : i_ind; //Filter index
        const int f_step = is_vertical ? o_step : i_step; //Filter step aka how many filter coeff multiplys and accumulations done in one loop
        const int t_col = is_vertical ? i_ind : o_ind; //trgt_buffer column
        const int t_step = is_vertical ? i_step : o_step; //Target buffer step aka how many target buffer values are calculated in one loop

        const int mask_shift = SCALER_CLIP( f_ind + f_step - filter_bound, 0, 8); //lane can hold 8 integers. Shift determines how many values should be masked from the end

        //Choose filter
        //const int *filter;
        //const int f_size = getFilter(&filter, is_upscaling, is_luma, phase, filter_phase);

        //Set trgt buffer val to zero on first loop over filter
        if (f_ind == 0) {
          //Process step number of elements at the same time
          memset(&trgt_row[t_col + trgt_offset], 0, sizeof(pic_data_t)*t_step); //trgt_row[t_col + trgt_offset] = 0;
        }

        //Move src pointer to correct position (correct column in vertical resampling)
        pic_data_t *src_col = src + (is_vertical ? i_ind : 0);

        //const int s_ind = SCALER_CLIP(ref_pos + f_ind - (filter_size >> 1) + 1, 0, src_size - 1); //src_buffer row/col index for cur resampling dir

        //Get the source incides of all the elements that are processed 
        pointer = clip_add_avx2(ref_pos + f_ind - (filter_size >> 1), adder, 0, src_size - 1);
        pointer = _mm256_permutevar8x32_epi32(pointer, order);

        if (is_vertical) {
          //Get indices that form column
          pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
        } else {
          //Get the smallest indice in pointer
          min = src_size - 1;
          smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
          smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
          min = _mm_extract_epi16(smallest_epi16, 0);
        }

        //Load src values to mem
        //TODO: Use correct inds
        //TODO: would set be faster here? Would need to concider the mask though.
        temp_mem = is_vertical 
          ? _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), src_col, pointer, _mm256_slli_epi32(base_mask, mask_shift), 1) 
          : _mm256_maskload_epi32(&src_col[min], _mm256_slli_epi32(base_mask, mask_shift));
        
        if (!is_vertical) {
          //Sort indices in the correct order
          decrese = _mm256_set1_epi32(min);
          pointer = _mm256_sub_epi32(pointer, decrese);
          temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
        }

        //Load filter
        temp_filter = _mm256_maskload_epi32(&getFilterCoeff(filter, filter_size, phase, f_ind), _mm256_slli_epi32(base_mask, mask_shift));

        //Multiply source by the filter coeffs and sum
        temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
        temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
        temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
        
        //Prepare to store values back to trgt_row
        switch (f_step) {
        case 4:
          //Workaround for extract on earlier visual studios
          trgt_row[t_col + trgt_offset] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0);//_mm256_extract_epi32(temp_mem, 0);
          break;

        case 8:
          trgt_row[t_col + trgt_offset] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);
          break;

        case 12:
          //Do filter operations for the remaining coefficients
          trgt_row[t_col + trgt_offset] = _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 1), 0);//_mm256_extract_epi32(temp_mem, 0) + _mm256_extract_epi32(temp_mem, 4);

          pointer = clip_add_avx2(ref_pos + f_ind - (filter_size >> 1) + 1, downsampling_adder, 0, src_size - 1);
          pointer = _mm256_permutevar8x32_epi32(pointer, order);

          if (is_vertical){
            pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
          } else {
            smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
            smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
            min = _mm_extract_epi16(smallest_epi16, 0);
          }
          
          const int tmp_mask_shift = 4; //The code assumes 12-tap so only do the remaining 4 taps
          temp_mem = is_vertical
            ? _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), src_col, pointer, _mm256_slli_epi32(base_mask, tmp_mask_shift), 1)
            : _mm256_maskload_epi32(&src_col[min], _mm256_slli_epi32(base_mask, tmp_mask_shift));

          if (!is_vertical) {
            decrese = _mm256_set1_epi32(min);
            pointer = _mm256_sub_epi32(pointer, decrese);
            temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
          }

          temp_filter = _mm256_maskload_epi32(&getFilterCoeff(filter, filter_size, phase, 8), _mm256_slli_epi32(base_mask, tmp_mask_shift));

          temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
          temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);
          temp_mem = _mm256_hadd_epi32(temp_mem, temp_mem);

          trgt_row[t_col + trgt_offset] += _mm_extract_epi32(_mm256_extracti128_si256(temp_mem, 0), 0); //_mm256_extract_epi32(temp_mem, 0);
          break;

        default:
          //Should not go here
          assert(0);
          break;
        }//END switch

        //trgt_row[t_col + trgt_offset] += getFilterCoeff(filter, filter_size, phase, f_ind) * src_col[s_ind * s_stride + src_offset];

        //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
        if (is_vertical && o_ind == outer_bound - 1) {
          trgt_row[t_col + trgt_offset] = SCALER_CLIP(is_upscaling ? (trgt_row[t_col + trgt_offset] + 2048) >> 12 : (trgt_row[t_col + trgt_offset] + 8192) >> 14, 0, 255);
        }
      }
    }
  }
}

//Set the default resample function
resample_block_step_func *const kvz_default_block_step_resample_func_avx2 = &DEFAULT_RESAMPLE_BLOCK_STEP_FUNC_AVX2; resample_func *const kvz_default_resample_func_avx2 = &DEFAULT_RESAMPLE_FUNC_AVX2;