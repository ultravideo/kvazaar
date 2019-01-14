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
#define ALT1_RESAMPLE_BLOCK_STEP_FUNC_AVX2 resampleBlockStep_avx2_v2
#define ALT2_RESAMPLE_BLOCK_STEP_FUNC_AVX2 resampleBlockStep_avx2_v3
#define DEFAULT_RESAMPLE_FUNC_AVX2 resample_avx2
#define ALT_RESAMPLE_FUNC_AVX2 resample2resampleBlockStep_alt2_avx2

#define B11011000 0xD8 //0b11011000

//Consts for permute
#define HI_epi128 0x13
#define LOW_epi128 0x02

#define SELECT_LOW_4_BITS 0xF

#if defined(__clang__)
#define FALLTHROUGH [[clang::fallthrough]]
#elif defined(__GNUC__) || defined(__GNUG__)
#define FALLTHROUGH __attribute__((fallthrough))
#else
#define FALLTHROUGH
#endif

// Clip sum of add_val to each epi32 of lane
static __m256i clip_add_avx2(const int add_val, __m256i lane, const int min, const int max)
{
 __m256i min_epi32 = _mm256_set1_epi32(min);
 __m256i max_epi32 = _mm256_set1_epi32(max);

 __m256i add_values_epi32 = _mm256_set1_epi32(add_val);
 add_values_epi32 = _mm256_add_epi32(add_values_epi32, lane);

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
   //src_row is not aligned to a 32-bit boundary so using load gives an error
   temp_mem = _mm256_loadu_si256((__m256i*)(&(src_row[min])));//_mm256_maskload_epi32(&src_row[min], _mm256_set1_epi32(0xf0000000));//_mm256_load_si256((__m256i*)(&(src_row[min])));
   
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
    //src_row is not aligned to a 32-bit boundary so using load gives an error
    temp_mem = _mm256_loadu_si256((__m256i*)(&(src_row[min]))); //_mm256_maskload_epi32(&src_row[min], _mm256_set1_epi32(0xf0000000));//_mm256_load_si256((__m256i*)&(src_row[min]));

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

static int num_distinct_ordered(int *list, int n)
{
  if(n <= 0){
    return 0;
  }

  int count = 1;
  for (int i = 1; i < n; i++) {
    if(list[i-1] != list[i] ){
      count++;
    }
  }

  return count;
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

        const int mask_shift = SCALER_CLIP(f_ind + f_step - filter_bound, 0, 8); //lane can hold 8 integers. Shift determines how many values should be masked from the end

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
        pointer = clip_add_avx2(ref_pos + f_ind - (filter_size >> 1) + 1, adder, 0, src_size - 1);
        pointer = _mm256_permutevar8x32_epi32(pointer, order);

        const int src_mask_shift = 8 - num_distinct_ordered((int*)&pointer, 8 - mask_shift); //Shift determines how many values should be masked from the end of the src to not over index

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
        //Use src mask on load
        temp_mem = is_vertical 
          ? _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), src_col, pointer, _mm256_slli_epi32(base_mask, mask_shift), 4) 
          : _mm256_maskload_epi32(&src_col[min], _mm256_slli_epi32(base_mask, src_mask_shift));
        
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
            ? _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), src_col, pointer, _mm256_slli_epi32(base_mask, tmp_mask_shift), 4)
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
        if (is_vertical && o_ind + o_step >= outer_bound) {
          trgt_row[t_col + trgt_offset] = SCALER_CLIP(is_upscaling ? (trgt_row[t_col + trgt_offset] + 2048) >> 12 : (trgt_row[t_col + trgt_offset] + 8192) >> 14, 0, 255);
        }
      }
    }
  }
}

//AVX2 cumulator. Sum epi32 values in n (max 8) groups in __m256i arrays, with the array divided between m (1,2,4,8) groups, into one output __m256i array
//Only enough least significant input vectors are used to fufill the n groups requirement.
static __m256i _mm256_accumulate_nxm_epi32(__m256i v7, __m256i v6, __m256i v5, __m256i v4, __m256i v3, __m256i v2, __m256i v1, __m256i v0, const unsigned n, const unsigned m)
{
  /*
          v7 v6 v5 v4 v3 v2 v1 v0
          |  |  |  |  |  |  |  |
  shuffle \ /   \ /   \ /   \ /
          / \   / \   / \   / \
  add     -+-   -+-   -+-   -+-
           |     |     |     |
  shuffle  \     /     \     /
            \   /       \   /
             \ /         \ /
             / \         / \
  add        -+-         -+-
              |           |
              ----   ------
  shuffle         \ /
                  / \
  add             -+-
                   |
                  out
  */

  if( (n == 0) || (m == 0)){
    return _mm256_setzero_si256();
  }

  if( m == 8 ){
    return v0;
  }

  __m256i tmp00, tmp01;
  __m256i tmp10, tmp11;
  __m256i tmp20, tmp21;
  __m256i tmp30, tmp31;
  __m256i add3, add2, add1, add0;


  //Set unused vals to zero
  if (m == 1) { //Stage One
    switch (n) {
    case 1:
      v1 = _mm256_setzero_si256();
      FALLTHROUGH;
    case 2:
      add1 = _mm256_setzero_si256();
      break;
    case 3:
      v3 = _mm256_setzero_si256();
      break;
    case 5:
      v5 = _mm256_setzero_si256();
      FALLTHROUGH;
    case 6:
      add3 = _mm256_setzero_si256();
      break;
    case 7:
      v7 = _mm256_setzero_si256();
      break;

    default:
      break;
    }//END switch 
  } else if(m < 4) { //Stage Two
    switch (n) {
    case 1:
      v0 = _mm256_inserti128_si256(v0, _mm_setzero_si128(), 1);
      FALLTHROUGH;
    case 2:
      add1 = _mm256_setzero_si256();
      break;
    case 3:
      v1 = _mm256_inserti128_si256(v1, _mm_setzero_si128(), 1);
      break;
    case 5:
      v2 = _mm256_inserti128_si256(v2, _mm_setzero_si128(), 1);
      FALLTHROUGH;
    case 6:
      add3 = _mm256_setzero_si256();
      break;
    case 7:
      v3 = _mm256_inserti128_si256(v3, _mm_setzero_si128(), 1);
      break;

    default:
      break;
    }//END switch
  } else { //Final Stage
    switch (n) {

    case 7:
      v1 = _mm256_blend_epi32(v1, _mm256_setzero_si256(), 0xC0);
      break;
    case 6:
      v1 = _mm256_blend_epi32(v1, _mm256_setzero_si256(), 0xF0);
      break;
    case 5:
      v1 = _mm256_blend_epi32(v1, _mm256_setzero_si256(), 0xFC);
      break;
    case 3:
      v0 = _mm256_blend_epi32(v0, _mm256_setzero_si256(), 0xC0);
      break;
    case 2:
      v0 = _mm256_blend_epi32(v0, _mm256_setzero_si256(), 0xF0);
      break;
    case 1:
      v0 = _mm256_blend_epi32(v0, _mm256_setzero_si256(), 0xFC);
      break;

    default:
      break;
    }//END switch
  }

  //Swap orders
  const __m256i swap_tmp = _mm256_set_epi32(7, 5, 6, 4, 3, 1, 2, 0);
  const __m256i swap_final = _mm256_set_epi32(7, 3, 5, 1, 6, 2, 4, 0);

  //First stage
  if (m < 2) {
    switch (n) {
    default:
      //8 is max allowed n value
      FALLTHROUGH;
    case 8: FALLTHROUGH;
    case 7:
      tmp31 = _mm256_permute2x128_si256(v7, v6, HI_epi128); tmp30 = _mm256_permute2x128_si256(v7, v6, LOW_epi128);
      add3 = _mm256_add_epi32(tmp31, tmp30);
      FALLTHROUGH;
    case 6: FALLTHROUGH;
    case 5: 
      tmp21 = _mm256_permute2x128_si256(v5, v4, HI_epi128); tmp20 = _mm256_permute2x128_si256(v5, v4, LOW_epi128);
      add2 = _mm256_add_epi32(tmp21, tmp20);
      FALLTHROUGH;
    case 4: FALLTHROUGH;
    case 3: 
      tmp11 = _mm256_permute2x128_si256(v3, v2, HI_epi128); tmp10 = _mm256_permute2x128_si256(v3, v2, LOW_epi128);
      add1 = _mm256_add_epi32(tmp11, tmp10);
      FALLTHROUGH;
    case 2: FALLTHROUGH;
    case 1:
      tmp01 = _mm256_permute2x128_si256(v1, v0, HI_epi128); tmp00 = _mm256_permute2x128_si256(v1, v0, LOW_epi128);
      add0 = _mm256_add_epi32(tmp01, tmp00);
      break;
    }//END switch
    
  } else if (m < 4){
    switch (n) {
    default:
      //8 is max allowed n value
      FALLTHROUGH;
    case 8: FALLTHROUGH;
    case 7:
      add3 = v3;
      FALLTHROUGH;
    case 6: FALLTHROUGH;
    case 5: 
      add2 = v2;
      FALLTHROUGH;
    case 4: FALLTHROUGH;
    case 3: 
      add1 = v1;
      FALLTHROUGH;
    case 2: FALLTHROUGH;
    case 1: 
      add0 = v0;
      break;
    }//END switch
  }

  //Second stage
  if (m < 4) {
    if (n > 4) {
      tmp11 = _mm256_unpackhi_epi64(add2, add3); tmp10 = _mm256_unpacklo_epi64(add2, add3);
      add2 = _mm256_add_epi32(tmp11, tmp10);
    }

    tmp01 = _mm256_unpackhi_epi64(add0, add1); tmp00 = _mm256_unpacklo_epi64(add0, add1);
    add0 = _mm256_add_epi32(tmp01, tmp00);

  } else {
    if (n > 4) {
      add2 = _mm256_permute4x64_epi64(v1, B11011000);
    }
    add0 = _mm256_permute4x64_epi64(v0, B11011000);
  }

  //Final stage
  if (n > 4) {
    add1 = _mm256_permutevar8x32_epi32(add2, swap_tmp);
  } else {
    add1 = _mm256_setzero_si256();
  }
  add0 = _mm256_permutevar8x32_epi32(add0, swap_tmp);

  tmp01 = _mm256_unpackhi_epi32(add0, add1); tmp00 = _mm256_unpacklo_epi32(add0, add1);

  add0 = _mm256_add_epi32(tmp01, tmp00);

  return _mm256_permutevar8x32_epi32(add0, swap_final);
}


//AVX2 cumulator. Sum epi32 values in 8 __m256i arrays into one output __m256i array
static __m256i _mm256_accumulate_8_epi32(__m256i v7, __m256i v6, __m256i v5, __m256i v4, __m256i v3, __m256i v2, __m256i v1, __m256i v0)
{  /*
          v7 v6 v5 v4 v3 v2 v1 v0
          |  |  |  |  |  |  |  |
  shuffle \ /   \ /   \ /   \ /
          / \   / \   / \   / \
  add     -+-   -+-   -+-   -+-   
           |     |     |     |
  shuffle  \     /     \     /
            \   /       \   /
             \ /         \ /
             / \         / \
  add        -+-         -+-
              |           |
              ----   ------
  shuffle         \ /
                  / \
  add             -+-
                   |
                  out
   */
  __m256i tmp00, tmp01;
  __m256i tmp10, tmp11;
  __m256i tmp20, tmp21;
  __m256i tmp30, tmp31;
  __m256i add3, add2, add1, add0;
  
  //Swap orders
  const __m256i swap_tmp = _mm256_set_epi32(7, 5, 6, 4, 3, 1, 2, 0);
  const __m256i swap_final = _mm256_set_epi32(7, 3, 5, 1, 6, 2, 4, 0);

  //First stage
  tmp31 = _mm256_permute2x128_si256(v7, v6, HI_epi128); tmp30 = _mm256_permute2x128_si256(v7, v6, LOW_epi128);
  tmp21 = _mm256_permute2x128_si256(v5, v4, HI_epi128); tmp20 = _mm256_permute2x128_si256(v5, v4, LOW_epi128);
  tmp11 = _mm256_permute2x128_si256(v3, v2, HI_epi128); tmp10 = _mm256_permute2x128_si256(v3, v2, LOW_epi128);
  tmp01 = _mm256_permute2x128_si256(v1, v0, HI_epi128); tmp00 = _mm256_permute2x128_si256(v1, v0, LOW_epi128);

  add3 = _mm256_add_epi32(tmp31, tmp30);
  add2 = _mm256_add_epi32(tmp21, tmp20);
  add1 = _mm256_add_epi32(tmp11, tmp10);
  add0 = _mm256_add_epi32(tmp01, tmp00);

  //Second stage
  tmp11 = _mm256_unpackhi_epi64(add2, add3); tmp10 = _mm256_unpacklo_epi64(add2, add3);
  tmp01 = _mm256_unpackhi_epi64(add0, add1); tmp00 = _mm256_unpacklo_epi64(add0, add1);

  add1 = _mm256_add_epi32(tmp11, tmp10);
  add0 = _mm256_add_epi32(tmp01, tmp00);

  //Final stage
  add1 = _mm256_permutevar8x32_epi32(add1, swap_tmp);
  add0 = _mm256_permutevar8x32_epi32(add0, swap_tmp); 

  tmp01 = _mm256_unpackhi_epi32(add0, add1); tmp00 = _mm256_unpacklo_epi32(add0, add1);

  add0 = _mm256_add_epi32(tmp01, tmp00);

  return _mm256_permutevar8x32_epi32(add0, swap_final);
}

//Avx load for loading n (max 8) epi32 values from memory
static __m256i _mm256_loadu_n_epi32(const int *src, const unsigned n)
{
  __m256i dst = _mm256_setzero_si256();

  switch (n) {
  case 0:
    break;

  case 1:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, src[0]);
    dst = _mm256_set_epi64x(0, 0, 0, (int64_t)src[0]);
    break;

  case 2:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, 0, src[1], src[0]);
    dst = _mm256_set_epi64x(0, 0, 0, (((int64_t)src[1]) << 32) | (int64_t)src[0]);
    break;

  case 3:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, src[2], src[1], src[0]);
    dst = _mm256_set_epi64x(0, 0, (int64_t)src[2], (((int64_t)src[1]) << 32) | (int64_t)src[0]);
    break;

  case 4:
    dst = _mm256_inserti128_si256(dst, _mm_loadu_si128((__m128i*)src), 0);
    break;

  case 5:
    //dst = _mm256_set_epi32(0, 0, 0, src[4], src[3], src[2], src[1], src[0]);
    dst = _mm256_set_epi64x(0, (int64_t)src[4], (((int64_t)src[3]) << 32) | (int64_t)src[2], (((int64_t)src[1]) << 32) | (int64_t)src[0]);
    break;

  case 6:
    //dst = _mm256_set_epi32(0, 0, src[5], src[4], src[3], src[2], src[1], src[0]);
    dst = _mm256_set_epi64x(0, (((int64_t)src[5]) << 32) | (int64_t)src[4], (((int64_t)src[3]) << 32) | (int64_t)src[2], (((int64_t)src[1]) << 32) | (int64_t)src[0]);
    break;

  case 7:
    //dst = _mm256_set_epi32(0, src[6], src[5], src[4], src[3], src[2], src[1], src[0]);
    dst = _mm256_set_epi64x((int64_t)src[6], (((int64_t)src[5])<<32) | (int64_t)src[4], (((int64_t)src[3]) << 32) | (int64_t)src[2], (((int64_t)src[1]) << 32) | (int64_t)src[0]);
    break;

  default:
    //Not a valid number of values to load. Only max 8 values can be loaded 
    FALLTHROUGH;
  case 8:
    dst = _mm256_loadu_si256((__m256i*)src);
    break;

  } //END Switch

  return dst;
}

//Read n (max 8) epi32i values from src specified by idx
static __m256i _mm256_gather_n_epi32(const int *src, const unsigned idx[8], const unsigned n)
{
  __m256i dst = _mm256_setzero_si256();

  switch(n){
  
  default:
    //Only eight values can be loaded
    //Fall through
  case 8:
    //dst = _mm256_set_epi32(src[idx[7]], src[idx[6]], src[idx[5]], src[idx[4]], src[idx[3]], src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x( (((int64_t)src[idx[7]])<<32) | (int64_t)src[idx[6]], (((int64_t)src[idx[5]])<<32) | (int64_t)src[idx[4]], (((int64_t)src[idx[3]])<<32) | (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 7:
    //dst = _mm256_set_epi32(0, src[idx[6]], src[idx[5]], src[idx[4]], src[idx[3]], src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x((int64_t)src[idx[6]], (((int64_t)src[idx[5]])<<32) | (int64_t)src[idx[4]], (((int64_t)src[idx[3]])<<32) | (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 6:
    //dst = _mm256_set_epi32(0, 0, src[idx[5]], src[idx[4]], src[idx[3]], src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x(0, (((int64_t)src[idx[5]])<<32) | (int64_t)src[idx[4]], (((int64_t)src[idx[3]])<<32) | (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 5:
    //dst = _mm256_set_epi32(0, 0, 0, src[idx[4]], src[idx[3]], src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x(0, (int64_t)src[idx[4]], (((int64_t)src[idx[3]]) << 32) | (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 4:
    //dst = _mm256_set_epi32(0, 0, 0, 0, src[idx[3]], src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x(0, 0, (((int64_t)src[idx[3]])<<32) | (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 3:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, src[idx[2]], src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x(0, 0, (int64_t)src[idx[2]], (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 2:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, 0, src[idx[1]], src[idx[0]]);
    dst = _mm256_set_epi64x(0, 0, 0, (((int64_t)src[idx[1]])<<32) | (int64_t)src[idx[0]]);
    break;
    
  case 1:
    //dst = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, src[idx[0]]);
    dst = _mm256_set_epi64x(0, 0, 0, (int64_t)src[idx[0]]);
    break;
    
  case 0:
    break;

  }//END Switch

  return dst;
}

//Avx store n epi32 from src to dst
static void _mm256_storeu_n_epi32(int *dst, __m256i src, const unsigned n)
{
  switch (n) {
  case 0:
    break;

  case 1:
    dst[0] = _mm_extract_epi32(_mm256_castsi256_si128(src), 0);
    break;

  case 2: {
    int64_t v1v0 = _mm_extract_epi64(_mm256_castsi256_si128(src), 0);

    dst[0] = (int)v1v0;
    dst[1] = (int)(v1v0 >> 32);

    break;
  }

  case 3: {
    __m128i tmp = _mm256_castsi256_si128(src);
    int64_t v1v0 = _mm_extract_epi64(tmp, 0);

    dst[0] = (int)v1v0;
    dst[1] = (int)(v1v0 >> 32);
    dst[2] = _mm_extract_epi32(tmp, 2);

    break;
  }

  case 4:
    _mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(src));
    break;

  case 5: {
    _mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(src));

    dst[4] = _mm_extract_epi32(_mm256_extracti128_si256(src, 1), 0);

    break;
  }

  case 6: {
    _mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(src));

    int64_t v1v0 = _mm_extract_epi64(_mm256_extracti128_si256(src, 1), 0);

    dst[4] = (int)v1v0;
    dst[5] = (int)(v1v0 >> 32);

    break;
  } 

  case 7: {
    _mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(src));

    __m128i tmp = _mm256_extracti128_si256(src, 1);
    int64_t v1v0 = _mm_extract_epi64(tmp, 0);

    dst[4] = (int)(v1v0);
    dst[5] = (int)(v1v0 >> 32);
    dst[6] = _mm_extract_epi32(tmp, 2);
    
    break;
  }

  default:
    //Only max 8 values can be stored
    FALLTHROUGH;
  case 8:
    _mm256_storeu_si256((__m256i*)dst, src);
    break;
  }
}

static void resampleBlockStep_avx2_v2(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical)
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
  const int trgt_bound = is_vertical ? inner_bound : outer_bound;
  const int filter_bound = is_vertical ? outer_bound : inner_bound;

  //Calculate outer and inner step so as to maximize lane/register usage:
  //  The accumulation can be done for 8 pixels at the same time
  //  If filter_size is 12, need to limit f_step to 8
  //  If filter_size is 4, can use only 4 ymm to hold eight pixels' filters
  const int t_step = 8; //Target buffer step aka how many target buffer values are calculated in one loop
  const int f_step = SCALER_MIN(filter_size, 8); //Filter step aka how many filter coeff multiplys and accumulations done in one loop
  //const int fm = 8 >> (f_step >> 1); //How many filter inds can be fit in one ymm

  const int o_step = is_vertical ? f_step : t_step; //Step size of outer loop. Adjust depending on how many values can be calculated concurrently with SIMD 
  const int i_step = is_vertical ? t_step : f_step; //Step size of inner loop. Adjust depending on how many values can be calculated concurrently with SIMD

  //const __m256i zero = _mm256_setzero_si256(); //Zero vector
  const __m256i scale_round = is_upscaling ? _mm256_set1_epi32(2048) : _mm256_set1_epi32(8192); //Rounding constant for normalizing pixel values to the correct range
  const int scale_shift = is_upscaling ? 12 : 14; //Amount of shift in the final pixel value normalization

  __m256i pointer, temp_trgt_epi32, decrese, filter_res_epi32;
  __m256i temp_mem[8], temp_filter[8];
  //const __m256i adder = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  const __m256i adderr = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  //const __m256i order = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  __m128i smallest_epi16;
  const __m256i multiplier_epi32 = _mm256_set1_epi32(s_stride); //Stride of src buffer
  int min = 0;

  __m256i t_ind_epi32, ref_pos_16_epi32, phase_epi32, ref_pos_epi32;
  const __m256i scale_epi32 = _mm256_set1_epi32(scale);
  const __m256i add_epi32 = _mm256_set1_epi32(add);
  const __m256i delta_epi32 = _mm256_set1_epi32(delta);
  const __m256i phase_mask = _mm256_set1_epi32(SELECT_LOW_4_BITS);

  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < (block_y + block_height); y++) {

    pic_data_t* src = is_vertical ? src_buffer->data : &src_buffer->data[y * src_buffer->width];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width];

    //Outer loop:
    //  if is_vertical -> loop over k (filter inds)
    //  if !is_vertical -> loop over x (block width)
    for (int o_ind = outer_init; o_ind < outer_bound; o_ind += o_step) {

      //const int t_ind = is_vertical ? y : o_ind; //trgt_buffer row/col index for cur resampling dir
      t_ind_epi32 = is_vertical ? _mm256_set1_epi32(y) : clip_add_avx2(o_ind, adderr, 0, outer_bound - 1);

      //Calculate reference position in src pic
      //const int ref_pos_16 = (int)((unsigned int)(t_ind * scale + add) >> shift) - delta;
      //const int phase = ref_pos_16 & 15;
      //const int ref_pos = ref_pos_16 >> 4;
      ref_pos_16_epi32 = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(t_ind_epi32, scale_epi32), add_epi32), shift), delta_epi32);
      phase_epi32 = _mm256_and_si256(ref_pos_16_epi32, phase_mask);
      ref_pos_epi32 = _mm256_srai_epi32(ref_pos_16_epi32, 4);

      const int *phase = (int*)&phase_epi32;
      const int *ref_pos = (int*)&ref_pos_epi32;

      //Inner loop:
      //  if is_vertical -> loop over x (block width)
      //  if !is_vertical -> loop over k (filter inds)-
      for (int i_ind = inner_init; i_ind < inner_bound; i_ind += i_step) {

        const int f_ind = is_vertical ? o_ind : i_ind; //Filter index
        const int t_col = is_vertical ? i_ind : o_ind; //trgt_buffer column

        //lane can hold 8 integers. f/t_num determines how many elements can be processed this loop (without going out of bounds)
        const unsigned f_num = SCALER_CLIP(filter_bound - f_ind, 0, f_step);
        const unsigned t_num = SCALER_CLIP(trgt_bound - t_col, 0, t_step);
        
        //Define a special case when filter_num is 4 that puts two "loops" into the same temp_mem[ind]
        const unsigned fm = f_num == 4 ? 2 : 1; //How many filter inds can be fit in one ymm

        //Get filter pixels for each ref_pos and do multiplication
        for (unsigned i = 0; i < t_num; i++) {
          const unsigned ind = i >> (fm >> 1); //Index to the temp avx2 vector arrays

          //Move src pointer to correct position (correct column in vertical resampling)
          pic_data_t *src_col = src + (is_vertical ? i_ind + i : 0);

          //Get the source incides of all the elements that are processed 
          pointer = clip_add_avx2(ref_pos[i] + f_ind - (filter_size >> 1) + 1, adderr, 0, src_size - 1);
          //pointer = _mm256_permutevar8x32_epi32(pointer, order);

          const unsigned src_num = num_distinct_ordered((int*)&pointer, f_num);

          if (is_vertical) {
            //Get indices that form a column
            pointer = _mm256_mullo_epi32(pointer, multiplier_epi32);
          }
          else {
            //Get the smallest indice in pointer
            min = src_size - 1;
            smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
            smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
            min = _mm_extract_epi16(smallest_epi16, 0);
          }

          //Load src values to mem
          if (fm == 1 || (i % 2) == 0) {
            temp_mem[ind] = is_vertical
              ? _mm256_gather_n_epi32(src_col, (unsigned*)&pointer, f_num)
              : _mm256_loadu_n_epi32(&src_col[min], src_num);
          } else {
            //Filter less than 8 elements at a time so can fit more "loops" in the same temp_mem[ind]
            temp_mem[ind] = _mm256_inserti128_si256(temp_mem[ind], _mm256_extracti128_si256(is_vertical
              ? _mm256_gather_n_epi32(src_col, (unsigned*)&pointer, f_num)
              : _mm256_loadu_n_epi32(&src_col[min], src_num), 0), 1);
          }

          if (!is_vertical) {
            //Sort indices in the correct order
            
            if (fm == 1 || (i % 2) == 0) {
              decrese = _mm256_set1_epi32(min);
              pointer = _mm256_sub_epi32(pointer, decrese);
              temp_mem[ind] = _mm256_permutevar8x32_epi32(temp_mem[ind], pointer);
            } else {
              //Only permute high 128bits
              decrese = _mm256_set1_epi32(min - 4);
              pointer = _mm256_inserti128_si256(pointer, _mm256_castsi256_si128(_mm256_sub_epi32(pointer, decrese)), 1);
              temp_mem[ind] = _mm256_blend_epi32(temp_mem[ind], _mm256_permutevar8x32_epi32(temp_mem[ind], pointer), 0xF0);
            }
          }

          //Load filter
          if (fm == 1 || (i % 2) == 0) {
            temp_filter[ind] = _mm256_loadu_n_epi32(&getFilterCoeff(filter, filter_size, phase[i], f_ind), f_num);
          } else {
            temp_filter[ind] = _mm256_inserti128_si256(temp_filter[ind], _mm256_castsi256_si128(_mm256_loadu_n_epi32(&getFilterCoeff(filter, filter_size, phase[i], f_ind), f_num)), 1);
          }

          //Multiply source by the filter coeffs and sum
          //Only do multiplication when fm number of "loops" are set or no more future "loops"
          if (fm == 1 || (i % 2) == 1 || i + 1 >= t_num ) {
            temp_mem[ind] = _mm256_mullo_epi32(temp_mem[ind], temp_filter[ind]);
          }
        }

        filter_res_epi32 = t_num == 8 && fm == 1
                             ? _mm256_accumulate_8_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0])
                             : _mm256_accumulate_nxm_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0], t_num, fm);

        //Sum filtered pixel values back to trgt_row so need to load the existing values (except for first pass)
        if (f_ind != 0) {
          temp_trgt_epi32 = _mm256_loadu_n_epi32(&trgt_row[t_col + trgt_offset], t_num);
          filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, temp_trgt_epi32);
        }
        
        //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
        if (is_vertical && o_ind + o_step >= outer_bound) {
          //trgt_row[t_col + trgt_offset] = SCALER_CLIP(is_upscaling ? (trgt_row[t_col + trgt_offset] + 2048) >> 12 : (trgt_row[t_col + trgt_offset] + 8192) >> 14, 0, 255);
          filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, scale_round);
          filter_res_epi32 = _mm256_srai_epi32(filter_res_epi32, scale_shift);
          filter_res_epi32 = clip_add_avx2(0, filter_res_epi32, 0, 255);
        }

        //Write back the new values for the current t_num pixels
        _mm256_storeu_n_epi32(&trgt_row[t_col + trgt_offset], filter_res_epi32, t_num);        
      }
    }
  }
}

static void resampleBlockStep_avx2_v3(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical)
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
  
  //Calculate outer and inner step so as to maximize lane/register usage:
  //  The accumulation can be done for 8 pixels at the same time
  //  If filter_size is 12, need to limit f_step to 8
  //  If filter_size is 4, can use 4 vector registers to hold eight pixels' filters
  const int t_step = 8; //Target buffer step aka how many target buffer values are calculated in one loop
  const int f_step = SCALER_MIN(filter_size, 8); //Filter step aka how many filter coeff multiplys and accumulations done in one loop
  
  const unsigned num_filter_parts = is_vertical ? 1 : (filter_size + 7) >> 3; //Number of loops needed to perform filtering in max 8 coeff chunks

  const int x_bound = block_x + block_width;
  const int y_bound = block_y + block_height;
  const unsigned i_bound = is_vertical ? filter_size : t_step;

  const unsigned s_stride = src_buffer->width;

  const __m256i scale_round = is_upscaling ? _mm256_set1_epi32(2048) : _mm256_set1_epi32(8192); //Rounding constant for normalizing pixel values to the correct range
  const int scale_shift = is_upscaling ? 12 : 14; //Amount of shift in the final pixel value normalization

  //Only filter size of max 12 supported
  assert(filter_size <= 12);

  const __m256i adderr = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

  const __m256i scale_epi32 = _mm256_set1_epi32(scale);
  const __m256i add_epi32 = _mm256_set1_epi32(add);
  const __m256i delta_epi32 = _mm256_set1_epi32(delta);
  const __m256i phase_mask = _mm256_set1_epi32(SELECT_LOW_4_BITS);

  __m256i temp_mem[8], temp_filter[8];
  __m256i filter_res_epi32;
  __m256i phase_epi32, ref_pos_epi32;


  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < y_bound; y++) {

    pic_data_t* src = is_vertical ? src_buffer->data + src_offset : &src_buffer->data[y * src_buffer->width + src_offset];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width + trgt_offset];
    
    //Outer loop:
    // loop over x (target block width)
    for (int x = block_x; x < x_bound; x += t_step) {

      const unsigned *phase;
      const int *ref_pos;
      int ref_pos_tmp;
      unsigned phase_tmp;

      const unsigned t_num = SCALER_CLIP(x_bound - x, 0, t_step);

      //Calculate reference position in src pic
      if (!is_vertical) {
        const __m256i t_ind_epi32 = clip_add_avx2(x, adderr, 0, x_bound - 1);
      
        const __m256i ref_pos_16_epi32 = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(t_ind_epi32, scale_epi32), add_epi32), shift), delta_epi32);
        phase_epi32 = _mm256_and_si256(ref_pos_16_epi32, phase_mask);
        ref_pos_epi32 = _mm256_srai_epi32(ref_pos_16_epi32, 4);

        phase = (unsigned*)&phase_epi32;
        ref_pos = (int*)&ref_pos_epi32;
      } else {
        const int ref_pos_16 = (int)((unsigned int)(y * scale + add) >> shift) - delta;
        phase_tmp = ref_pos_16 & 15; phase = &phase_tmp;
        ref_pos_tmp = ref_pos_16 >> 4; ref_pos = &ref_pos_tmp;
      }                                               
      
      //Loop over filter segments if filter is longer than f_step 
      for (unsigned filter_part = 0; filter_part < num_filter_parts; filter_part++) {
        
        const int f_ind = filter_part * f_step; //Filter index

        //lane can hold 8 integers. f/t_num determines how many elements can be processed this loop (without going out of bounds)
        const unsigned f_num = SCALER_CLIP(filter_size - f_ind, 0, f_step);
        const int fm = f_num == 4 ? 2 : 1; //How many filter inds can be fit in one ymm

        //Inner loop:
        // if is_vertical -> loop over k (filter inds)
        // if !is_vertical -> loop over x (in t_step)
        for (unsigned i = 0; i < i_bound; i++) {

          if (!is_vertical) {
            //Define a special case when filter_num is 4 that puts two "loops" into the same temp_mem[ind]

            const int ind = i >> (fm >> 1); //Index to the temp avx2 vector arrays

            //Get the source incides of all the elements that are processed 
            const __m256i pointer = clip_add_avx2(ref_pos[i] + f_ind - (filter_size >> 1) + 1, adderr, 0, src_size - 1);

            const int src_num = num_distinct_ordered((int*)&pointer, f_num);
            
            const int min = _mm_extract_epi32(_mm256_castsi256_si128(pointer), 0);

            if (fm == 1 || (i % 2) == 0) {
              //Load src values to mem
              temp_mem[ind] = _mm256_loadu_n_epi32(&src[min], src_num);

              //Get correct pixel values pased on pointer
              const __m256i decrese = _mm256_set1_epi32(min);
              const __m256i perm = _mm256_sub_epi32(pointer, decrese);
              temp_mem[ind] = _mm256_permutevar8x32_epi32(temp_mem[ind], perm);

              //Load filter
              temp_filter[ind] = _mm256_loadu_n_epi32(&getFilterCoeff(filter, filter_size, phase[i], f_ind), f_num);

            }
            else {
              //Filter less than 8 elements at a time so can fit more "loops" in the same temp_mem[ind]
              temp_mem[ind] = _mm256_inserti128_si256(temp_mem[ind], _mm256_extracti128_si256(_mm256_loadu_n_epi32(&src[min], src_num), 0), 1);

              //Only permute high 128bits
              const __m256i decrese = _mm256_set1_epi32(min - 4);
              const __m256i perm = _mm256_inserti128_si256(pointer, _mm256_castsi256_si128(_mm256_sub_epi32(pointer, decrese)), 1);
              temp_mem[ind] = _mm256_blend_epi32(temp_mem[ind], _mm256_permutevar8x32_epi32(temp_mem[ind], perm), 0xF0);

              temp_filter[ind] = _mm256_inserti128_si256(temp_filter[ind], _mm256_castsi256_si128(_mm256_loadu_n_epi32(&getFilterCoeff(filter, filter_size, phase[i], f_ind), f_num)), 1);
            }

            //Multiply source by the filter coeffs and sum
            //Only do multiplication when fm number of "loops" are set or no more future "loops"
            if (fm == 1 || (i % 2) == 1 || i + 1 >= t_num) {
              temp_mem[ind] = _mm256_mullo_epi32(temp_mem[ind], temp_filter[ind]);
            }
          }
          else {
            //Get src row corresponding to cur filter index i
            const int s_ind = SCALER_CLIP(ref_pos[0] + (int)i - (filter_size >> 1) + 1, 0, src_size - 1);
            *temp_mem = _mm256_loadu_n_epi32(&src[s_ind * s_stride + x], t_num);

            *temp_filter = _mm256_set1_epi32(getFilterCoeff(filter, filter_size, *phase, i));

            *temp_mem = _mm256_mullo_epi32(*temp_mem, *temp_filter);

            if (i == 0) {
              filter_res_epi32 = *temp_mem;
            }
            else {
              filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, *temp_mem);
            }
          }
        }

        if (!is_vertical) {
          //Accumulate multiplication step results to the final filtered values
          if (f_ind == 0) {
            filter_res_epi32 = t_num == 8 && fm == 1
              ? _mm256_accumulate_8_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0])
              : _mm256_accumulate_nxm_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0], t_num, fm);
          } else {
            filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, t_num == 8 && fm == 1
              ? _mm256_accumulate_8_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0])
              : _mm256_accumulate_nxm_epi32(temp_mem[7], temp_mem[6], temp_mem[5], temp_mem[4], temp_mem[3], temp_mem[2], temp_mem[1], temp_mem[0], t_num, fm));
          }
        }
      }

      //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
      if (is_vertical) {
        filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, scale_round);
        filter_res_epi32 = _mm256_srai_epi32(filter_res_epi32, scale_shift);
        filter_res_epi32 = clip_add_avx2(0, filter_res_epi32, 0, 255);
      }

      //Write back the new values for the current t_num pixels
      _mm256_storeu_n_epi32(&trgt_row[x], filter_res_epi32, t_num);

    }
  }
}

static void resample2resampleBlockStep_alt2_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
{
  pic_buffer_t* tmp = kvz_newPictureBuffer(param->trgt_width + param->trgt_padding_x, param->src_height + param->src_padding_y, 0);

  //Vertical resampling
  kvz_alt2_block_step_resample_func_avx2(buffer, tmp, 0, 0, 0, 0, param->trgt_width + param->trgt_padding_x, param->src_height + param->src_padding_y, param, is_upscaling, is_luma, 0);

  //Horizontal resampling
  kvz_alt2_block_step_resample_func_avx2(tmp, buffer, 0, 0, 0, 0, param->trgt_width + param->trgt_padding_x, param->trgt_height + param->trgt_padding_y, param, is_upscaling, is_luma, 1);

  kvz_deallocatePictureBuffer(tmp);
}



//Takes in six vectors of data and filter coeffs. Apply filter and return results for eight pixels at a time.
//Data layout:
//  data0[X]:=|D3X D1X|D2X D0X|
//  data1[X]:=|D7X D5X|D6X D4X|
//  where DYX:={d_(Y,X+0), d_(Y,X+1), d_(Y,X+2), d_(Y,X+3)} and d_(x,y) is the yth sample for the xth pixel given as a 16-bit value
//
//  filter0[X]:=|F3X F1X|F2X F0X|
//  filter1[X]:=|F7X F5X|F6X F4X|
//  where FYX:={f_(Y,X+3), f_(Y,X+2), f_(Y,X+1), f_(Y,X+0)} and f_(x,y) is the yth filter coeff for the xth pixel given as a 16-bit value
//
//Returns:
//  |r_7 r_6 r_5 r_4|r_3 r_2 r_1 r_0|
//  where r_x = is the filterd pixel value of the xth pixel given as a 32-bit value
//
//Only dataX[i], where 0 < i < filter_rounds, is used. Data and filter values should be given as 
static __m256i apply_filter_4x4_dual(const __m256i* data0, const __m256i* data1, const __m256i* filter0, const __m256i* filter1, const unsigned filter_rounds)
{
  
  __m256i res;
  const __m256i perm = _mm256_set_epi32(7, 5, 3, 1, 6, 2, 4, 0);

  //Filtering done four taps at a time for four pixels, so repeate until desired number of taps performed
  __m256i tmp0 = _mm256_madd_epi16(data0[0], filter0[0]);
  __m256i tmp1 = _mm256_madd_epi16(data1[0], filter1[0]);

  for (size_t i = 1; i < filter_rounds; i++)
  {
    //Do multiplication and sum with results for previous taps
    tmp0 = _mm256_add_epi32(tmp0, _mm256_madd_epi16(data0[i], filter0[i]));
    tmp1 = _mm256_add_epi32(tmp1, _mm256_madd_epi16(data1[i], filter1[i]));
  }

  //Need to do one more round of additions to get final values
  //  Blend to get |dd bb|cc aa|, |hh ff|gg ee| => |hd fb|gc ea|, |dh bf|cg ae| for the sum (need suffle to get the latter vector to the correct order)
  res = _mm256_add_epi32(_mm256_blend_epi32(tmp0, tmp1, /*1010 1010*/0xAA), _mm256_shuffle_epi32(_mm256_blend_epi32(tmp0, tmp1, /*0101 0101*/0x55), /*1011 0001*/0xB1));

  //Finally permute to get the order right
  return _mm256_permutevar8x32_epi32(res, perm);
}

static void resampleBlockStep_avx2_v3(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical)
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

  //Check that scaling ratio is not too large for pre-loading to fail
  //Need to have four pixel positions worth of samples at any time, so the reference/sample pixel position distance between x_i and x_i+3 cannot exceed 32 (the max number of continous samples that are pre-loaded)
  //if ((((int)((unsigned int)(4 * scale + add) >> shift) - delta) >> 4) + (filter_size >> 1) + 1 > 32)
  //{
  //  //Fallback on another resampling function
  //  resampleBlockStep_avx2_v3(src_buffer, trgt_buffer, src_offset, trgt_offset, block_x, block_y, block_width, block_height, param, is_upscaling, is_luma, is_vertical);
  //  return;
  //}

  //Only filter size of max 12 supported
  assert(filter_size <= 12);

  //Calculate outer and inner step so as to maximize lane/register usage:
  //  The accumulation can be done for 8 pixels at the same time
  const int t_step = 8; //Target buffer step aka how many target buffer values are calculated in one loop
  const int f_step = 4;//SCALER_MIN(filter_size, 8); //Filter step aka how many filter coeff multiplys and accumulations done in one loop

  const unsigned filter_parts = filter_size / f_step; //How many f_step sized chunks the filter size contains
  const unsigned num_filter_parts = is_vertical ? 0 : filter_parts;//(filter_size + 7) >> 3; //Number of loops needed to perform filtering

  const int x_bound = block_x + block_width;
  const int y_bound = block_y + block_height;
  const unsigned i_bound = is_vertical ? filter_size : 1;

  const unsigned s_stride = src_buffer->width;

  const __m256i scale_round = is_upscaling ? _mm256_set1_epi32(2048) : _mm256_set1_epi32(8192); //Rounding constant for normalizing pixel values to the correct range
  const int scale_shift = is_upscaling ? 12 : 14; //Amount of shift in the final pixel value normalization

  
  const __m256i adderr = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

  const __m256i scale_epi32 = _mm256_set1_epi32(scale);
  const __m256i add_epi32 = _mm256_set1_epi32(add);
  const __m256i delta_epi32 = _mm256_set1_epi32(delta);
  const __m256i phase_mask = _mm256_set1_epi32(SELECT_LOW_4_BITS);

  __m256i filters_epi16[12]; //Contain preloaded filters. Data layouts: |Ph(X+3) Ph(X+1)|Ph(X+2) Ph(X)| or |Ph(X+1)_1 Ph(X)_1|Ph(X+1)_0 Ph(X)_0| or |Ph(X+1)_0 Ph(X)_1|Ph(X)_2 Ph(X)_0| and |Ph(X+2)_1 Ph(X+1)_2|Ph(X+2)_0 Ph(X+1)_1| and |Ph(X+3)_2 Ph(X+3)_0|Ph(X+3)_1 Ph(X+2)_2|
  __m256i temp_mem[12], temp_filter[12];
  __m256i filter_res_epi32;
  __m256i phase_epi32, ref_pos_epi32;
  __m256i data0[3], data1[3], filter0[3], filter1[3];
  
  const unsigned phase_map[4][4] = { { 0, 2, 1, 3 }, { 3, 0, 2, 1 }, { 1, 3, 0, 2 }, { 2, 1, 3, 0} }; //Phase map for 4 tap filter and 12 tap filter
  const unsigned filter_map[4][3] = { { 0, 0, 0 }, { 0, 1, 1 }, { 1, 1, 2 }, { 2, 2, 2 } }; //Filter array map for 12 tap

  //Pre-load filters
  if (filter_size <= 4) {
    
    filters_epi16[0] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[0]), _mm256_loadu_si256((__m256i*)&filter[8]));
    filters_epi16[1] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[16]), _mm256_loadu_si256((__m256i*)&filter[24]));
    filters_epi16[2] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[32]), _mm256_loadu_si256((__m256i*)&filter[40]));
    filters_epi16[3] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[48]), _mm256_loadu_si256((__m256i*)&filter[56]));
  }
  if (4 < filter_size && filter_size <= 8) {
    
    filters_epi16[4] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[64]), _mm256_loadu_si256((__m256i*)&filter[72]));
    filters_epi16[5] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[80]), _mm256_loadu_si256((__m256i*)&filter[88]));
    filters_epi16[6] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[96]), _mm256_loadu_si256((__m256i*)&filter[104]));
    filters_epi16[7] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[112]), _mm256_loadu_si256((__m256i*)&filter[120]));
  }
  if( 8 < filter_size ) {

    filters_epi16[8] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[128]), _mm256_loadu_si256((__m256i*)&filter[136]));
    filters_epi16[9] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[144]), _mm256_loadu_si256((__m256i*)&filter[152]));
    filters_epi16[10] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[160]), _mm256_loadu_si256((__m256i*)&filter[168]));
    filters_epi16[11] = _mm256_packs_epi32(_mm256_loadu_si256((__m256i*)&filter[176]), _mm256_loadu_si256((__m256i*)&filter[184]));
  }

  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < y_bound; y++) {

    pic_data_t* src = is_vertical ? src_buffer->data + src_offset : &src_buffer->data[y * src_buffer->width + src_offset];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width + trgt_offset];

    //Outer loop:
    // loop over x (target block width)
    for (int x = block_x; x < x_bound; x += t_step) {

      const unsigned *phase;
      const int *ref_pos;
      int ref_pos_tmp;
      unsigned phase_tmp;

      const unsigned t_num = SCALER_CLIP(x_bound - x, 0, t_step);

      //Calculate reference position in src pic
      if (!is_vertical) {
        const __m256i t_ind_epi32 = clip_add_avx2(x, adderr, 0, x_bound - 1);

        const __m256i ref_pos_16_epi32 = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(t_ind_epi32, scale_epi32), add_epi32), shift), delta_epi32);
        phase_epi32 = _mm256_and_si256(ref_pos_16_epi32, phase_mask);
        ref_pos_epi32 = _mm256_srai_epi32(ref_pos_16_epi32, 4);

        phase = (unsigned*)&phase_epi32;
        ref_pos = (int*)&ref_pos_epi32;
      }
      else {
        const int ref_pos_16 = (int)((unsigned int)(y * scale + add) >> shift) - delta;
        phase_tmp = ref_pos_16 & 15; phase = &phase_tmp;
        ref_pos_tmp = ref_pos_16 >> 4; ref_pos = &ref_pos_tmp;
      }

      //Loop over filter segments and load needed data 
      for (unsigned filter_part = 0; filter_part < num_filter_parts; filter_part++) {

        const int f_ind = filter_part * f_step; //Filter index
        
        //Load filter
        if (filter_size <= 4) {
          filter0[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[phase[3] >> 2], phase_map[0][phase[3] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[1] >> 2], phase_map[0][phase[1] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[2] >> 2], phase_map[0][phase[2] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[0] >> 2], phase_map[0][phase[0] % 4]));
          filter1[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[phase[7] >> 2], phase_map[0][phase[7] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[5] >> 2], phase_map[0][phase[5] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[6] >> 2], phase_map[0][phase[6] % 4]),
                                                   _mm256_extract_epi64(filters_epi16[phase[4] >> 2], phase_map[0][phase[4] % 4]));
        } else if (filter_size <= 8) {
          filter0[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[phase[3] >> 1], (phase[3] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[1] >> 1], (phase[1] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[2] >> 1], (phase[2] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[0] >> 1], (phase[0] % 2) + filter_part));
          filter1[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[phase[7] >> 1], (phase[7] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[5] >> 1], (phase[5] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[6] >> 1], (phase[6] % 2) + filter_part),
                                                   _mm256_extract_epi64(filters_epi16[phase[4] >> 1], (phase[4] % 2) + filter_part));
        } else {
          filter0[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[(phase[3] >> 2) + filter_map[phase[3] % 4][filter_part]], phase_map[phase[3] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[1] >> 2) + filter_map[phase[1] % 4][filter_part]], phase_map[phase[1] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[2] >> 2) + filter_map[phase[2] % 4][filter_part]], phase_map[phase[2] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[0] >> 2) + filter_map[phase[0] % 4][filter_part]], phase_map[phase[0] % 4][filter_part]));
          filter1[filter_part] = _mm256_set_epi64x(_mm256_extract_epi64(filters_epi16[(phase[7] >> 2) + filter_map[phase[7] % 4][filter_part]], phase_map[phase[7] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[5] >> 2) + filter_map[phase[5] % 4][filter_part]], phase_map[phase[5] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[6] >> 2) + filter_map[phase[6] % 4][filter_part]], phase_map[phase[6] % 4][filter_part]),
                                                   _mm256_extract_epi64(filters_epi16[(phase[4] >> 2) + filter_map[phase[4] % 4][filter_part]], phase_map[phase[4] % 4][filter_part]));
        }
        
        //Load src samples in four pixel chunks into registers
        temp_mem[f_ind + 0] = _mm256_loadu2_m128i((__m128i*)&src[ref_pos[1] + f_ind], (__m128i*)&src[ref_pos[0] + f_ind]);
        temp_mem[f_ind + 1] = _mm256_loadu2_m128i((__m128i*)&src[ref_pos[3] + f_ind], (__m128i*)&src[ref_pos[2] + f_ind]);
        temp_mem[f_ind + 2] = _mm256_loadu2_m128i((__m128i*)&src[ref_pos[5] + f_ind], (__m128i*)&src[ref_pos[4] + f_ind]);
        temp_mem[f_ind + 3] = _mm256_loadu2_m128i((__m128i*)&src[ref_pos[7] + f_ind], (__m128i*)&src[ref_pos[6] + f_ind]);

        //Pack samples into 16-bit values. Data order will be |P3|P1|P2|P0| etc.
        data0[filter_part] = _mm256_packus_epi32(temp_mem[f_ind + 0], temp_mem[f_ind + 1]);
        data1[filter_part] = _mm256_packus_epi32(temp_mem[f_ind + 2], temp_mem[f_ind + 3]);

      }
      
      //Inner loop:
      // if is_vertical -> loop over k (filter inds)
      // if !is_vertical -> loop over x (in t_step)
      for (unsigned i = 0; i < i_bound; i++) {

        if (!is_vertical) {

          filter_res_epi32 = apply_filter_4x4_dual(data0, data1, filter0, filter1, num_filter_parts);

        } else {
          //Get src row corresponding to cur filter index i
          const int s_ind = SCALER_CLIP(ref_pos[0] + (int)i - (filter_size >> 1) + 1, 0, src_size - 1);
          *temp_mem = _mm256_loadu_n_epi32(&src[s_ind * s_stride + x], t_num);

          *temp_filter = _mm256_set1_epi32(getFilterCoeff(filter, filter_size, *phase, i));

          *temp_mem = _mm256_mullo_epi32(*temp_mem, *temp_filter);

          if (i == 0) {
            filter_res_epi32 = *temp_mem;
          } else {
            filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, *temp_mem);
          }
        }
      }

      //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
      if (is_vertical) {
        filter_res_epi32 = _mm256_add_epi32(filter_res_epi32, scale_round);
        filter_res_epi32 = _mm256_srai_epi32(filter_res_epi32, scale_shift);
        filter_res_epi32 = clip_add_avx2(0, filter_res_epi32, 0, 255);
      }

      //Write back the new values for the current t_num pixels
      _mm256_storeu_n_epi32(&trgt_row[x], filter_res_epi32, t_num);

    }
  }

}

//static __m256i _mm256_gather_n_epi32_v2(const int *src, unsigned idx[8], unsigned n)
//{
//  __m256i dst = _mm256_setzero_si256();
//
//  switch (n) {
//
//  default:
//    //Only eight values can be loaded
//    //Fall through
//  case 8:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 1), src[idx[7]], 3), 1);
//    //Fall through
//  case 7:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 1), src[idx[6]], 2), 1);
//    //Fall through
//  case 6:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 1), src[idx[5]], 1), 1);
//    //Fall through
//  case 5:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 1), src[idx[4]], 0), 1);
//    //Fall through
//  case 4:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 0), src[idx[3]], 3), 0);
//    //Fall through
//  case 3:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 0), src[idx[2]], 2), 0);
//    //Fall through
//  case 2:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 0), src[idx[1]], 1), 0);
//    //Fall through
//  case 1:
//    dst = _mm256_inserti128_si256(dst, _mm_insert_epi32(_mm256_extracti128_si256(dst, 0), src[idx[0]], 0), 0);
//    //Fall through
//  case 0:
//    break;
//
//  }//END Switch
//
//  return dst;
//}

//int test_avx()
//{
//  unsigned in[18] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
//  unsigned in2[8] = {2, 3, 5, 7, 9, 11, 13, 17};
//  unsigned out[8] = {0,0,0,0,0,0,0,0};
//  unsigned out2[4] = {0,0,0,0};
//  unsigned out3 = 0;
//  __m256i z = _mm256_setzero_si256();
//  __m256i v1 = _mm256_gather_n_epi32(in, in2, 8);
//  __m256i v2 = _mm256_loadu_n_epi32(in, 8);
//  _mm256_storeu_n_epi32(out, v2, 8);
//  _mm256_storeu_n_epi32(out2, v1, 4);
//
//  __m256i v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 2, 1);
//  v3 = _mm256_accumulate_8_epi32(z, z, z, z, z, z, v1, v2);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 8, 1);
//
//  v3 = _mm256_accumulate_nxm_epi32(v2, v2, v1, v1, v2, v1, v1, v2, 8, 2);
//
//  v1 = _mm256_gather_n_epi32(in, in2+4, 4);
//  v2 = _mm256_loadu_n_epi32(out2, 4);
//
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 4, 2);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 1, 2);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 3, 2);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 2, 2);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 8, 4);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 7, 4);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 6, 4);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 2, 4);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 3, 4);
//  v3 = _mm256_accumulate_nxm_epi32(z, z, z, z, z, z, v1, v2, 1, 4);
//
//  
//  return _mm_extract_epi32( _mm256_extracti128_si256(v1, 0), 0) + _mm_extract_epi32(_mm256_extracti128_si256(v2, 1), 0) + out[3];
//}

//Set the default resample function
resample_block_step_func *const kvz_default_block_step_resample_func_avx2 = &DEFAULT_RESAMPLE_BLOCK_STEP_FUNC_AVX2;
resample_block_step_func *const kvz_alt1_block_step_resample_func_avx2 = &ALT1_RESAMPLE_BLOCK_STEP_FUNC_AVX2;
resample_block_step_func *const kvz_alt2_block_step_resample_func_avx2 = &ALT2_RESAMPLE_BLOCK_STEP_FUNC_AVX2;
resample_func *const kvz_default_resample_func_avx2 = &DEFAULT_RESAMPLE_FUNC_AVX2;
resample_func *const kvz_alt_resample_func_avx2 = &ALT_RESAMPLE_FUNC_AVX2;
