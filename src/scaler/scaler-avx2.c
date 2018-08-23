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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>

#define DEFAULT_RESAMPLE_BLOCK_STEP_FUNC_AVX2 resampleBlockStep_avx2

#define B11011000 0xD8 //0b11011000

// Used only to clip values
static __m256i clip_avx2(int ref_pos, int src_width, int size, __m256i adder)
{
 __m256i zeros = _mm256_setzero_si256();
 __m256i max_epi32 = _mm256_set1_epi32(src_width - 1);

 int temp = ref_pos - (size >> 1) + 1;

 __m256i values_epi32 = _mm256_set1_epi32(temp);
 values_epi32 = _mm256_add_epi32(values_epi32, adder);

 // Compare integers so value is in range 0 - (src-width-1)
 values_epi32 = _mm256_max_epi32(zeros, values_epi32);
 values_epi32 = _mm256_min_epi32(values_epi32, max_epi32);
 
 return values_epi32;
}

//Used for clipping values
static int clip(int val, int min, int max)
{

 if (val <= min)
  return min;
 if (val >= max)
  return max;

 return val;
}

//Define filters for scaling operations
//Values from SHM
static const int downFilter[8][16][12] = {
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

static const int lumaUpFilter[16][8] = {
 { 0, 0, 0, 64, 0, 0, 0, 0 },
 { 0, 1, -3, 63, 4, -2, 1, 0 },
 { -1, 2, -5, 62, 8, -3, 1, 0 },
 { -1, 3, -8, 60, 13, -4, 1, 0 },
 { -1, 4, -10, 58, 17, -5, 1, 0 },
 { -1, 4, -11, 52, 26, -8, 3, -1 },
 { -1, 3, -9, 47, 31, -10, 4, -1 },
 { -1, 4, -11, 45, 34, -10, 4, -1 },
 { -1, 4, -11, 40, 40, -11, 4, -1 },
 { -1, 4, -10, 34, 45, -11, 4, -1 },
 { -1, 4, -10, 31, 47, -9, 3, -1 },
 { -1, 3, -8, 26, 52, -11, 4, -1 },
 { 0, 1, -5, 17, 58, -10, 4, -1 },
 { 0, 1, -4, 13, 60, -8, 3, -1 },
 { 0, 1, -3, 8, 62, -5, 2, -1 },
 { 0, 1, -2, 4, 63, -3, 1, 0 }
};

static const int chromaUpFilter[16][4] = {
 { 0, 64, 0, 0 },
 { -2, 62, 4, 0 },
 { -2, 58, 10, -2 },
 { -4, 56, 14, -2 },
 { -4, 54, 16, -2 },
 { -6, 52, 20, -2 },
 { -6, 46, 28, -4 },
 { -4, 42, 30, -4 },
 { -4, 36, 36, -4 },
 { -4, 30, 42, -4 },
 { -4, 28, 46, -6 },
 { -2, 20, 52, -6 },
 { -2, 16, 54, -4 },
 { -2, 14, 56, -4 },
 { -2, 10, 58, -2 },
 { 0, 4, 62, -2 }
};

static const int downFilter1D[8][16 * 12] = {
  // ratio <= 20/19
  {
    0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, -6, 127, 7, -2, 0, 0, 0, 0,
    0, 0, 0, 3, -12, 125, 16, -5, 1, 0, 0, 0,
    0, 0, 0, 4, -16, 120, 26, -7, 1, 0, 0, 0,
    0, 0, 0, 5, -18, 114, 36, -10, 1, 0, 0, 0,
    0, 0, 0, 5, -20, 107, 46, -12, 2, 0, 0, 0,
    0, 0, 0, 5, -21, 99, 57, -15, 3, 0, 0, 0,
    0, 0, 0, 5, -20, 89, 68, -18, 4, 0, 0, 0,
    0, 0, 0, 4, -19, 79, 79, -19, 4, 0, 0, 0,
    0, 0, 0, 4, -18, 68, 89, -20, 5, 0, 0, 0,
    0, 0, 0, 3, -15, 57, 99, -21, 5, 0, 0, 0,
    0, 0, 0, 2, -12, 46, 107, -20, 5, 0, 0, 0,
    0, 0, 0, 1, -10, 36, 114, -18, 5, 0, 0, 0,
    0, 0, 0, 1, -7, 26, 120, -16, 4, 0, 0, 0,
    0, 0, 0, 1, -5, 16, 125, -12, 3, 0, 0, 0,
    0, 0, 0, 0, -2, 7, 127, -6, 2, 0, 0, 0,
  },
  // 20/19 < ratio <= 5/4
  {
    0, 2, 0, -14, 33, 86, 33, -14, 0, 2, 0, 0,
    0, 1, 1, -14, 29, 85, 38, -13, -1, 2, 0, 0,
    0, 1, 2, -14, 24, 84, 43, -12, -2, 2, 0, 0,
    0, 1, 2, -13, 19, 83, 48, -11, -3, 2, 0, 0,
    0, 0, 3, -13, 15, 81, 53, -10, -4, 3, 0, 0,
    0, 0, 3, -12, 11, 79, 57, -8, -5, 3, 0, 0,
    0, 0, 3, -11, 7, 76, 62, -5, -7, 3, 0, 0,
    0, 0, 3, -10, 3, 73, 65, -2, -7, 3, 0, 0,
    0, 0, 3, -9, 0, 70, 70, 0, -9, 3, 0, 0,
    0, 0, 3, -7, -2, 65, 73, 3, -10, 3, 0, 0,
    0, 0, 3, -7, -5, 62, 76, 7, -11, 3, 0, 0,
    0, 0, 3, -5, -8, 57, 79, 11, -12, 3, 0, 0,
    0, 0, 3, -4, -10, 53, 81, 15, -13, 3, 0, 0,
    0, 0, 2, -3, -11, 48, 83, 19, -13, 2, 1, 0,
    0, 0, 2, -2, -12, 43, 84, 24, -14, 2, 1, 0,
    0, 0, 2, -1, -13, 38, 85, 29, -14, 1, 1, 0,
  },
  // 5/4 < ratio <= 5/3
  {
    0, 5, -6, -10, 37, 76, 37, -10, -6, 5, 0, 0,
    0, 5, -4, -11, 33, 76, 40, -9, -7, 5, 0, 0,
    -1, 5, -3, -12, 29, 75, 45, -7, -8, 5, 0, 0,
    -1, 4, -2, -13, 25, 75, 48, -5, -9, 5, 1, 0,
    -1, 4, -1, -13, 22, 73, 52, -3, -10, 4, 1, 0,
    -1, 4, 0, -13, 18, 72, 55, -1, -11, 4, 2, -1,
    -1, 4, 1, -13, 14, 70, 59, 2, -12, 3, 2, -1,
    -1, 3, 1, -13, 11, 68, 62, 5, -12, 3, 2, -1,
    -1, 3, 2, -13, 8, 65, 65, 8, -13, 2, 3, -1,
    -1, 2, 3, -12, 5, 62, 68, 11, -13, 1, 3, -1,
    -1, 2, 3, -12, 2, 59, 70, 14, -13, 1, 4, -1,
    -1, 2, 4, -11, -1, 55, 72, 18, -13, 0, 4, -1,
    0, 1, 4, -10, -3, 52, 73, 22, -13, -1, 4, -1,
    0, 1, 5, -9, -5, 48, 75, 25, -13, -2, 4, -1,
    0, 0, 5, -8, -7, 45, 75, 29, -12, -3, 5, -1,
    0, 0, 5, -7, -9, 40, 76, 33, -11, -4, 5, 0,
  },
  // 5/3 < ratio <= 2
  {
    2, -3, -9, 6, 39, 58, 39, 6, -9, -3, 2, 0,
    2, -3, -9, 4, 38, 58, 43, 7, -9, -4, 1, 0,
    2, -2, -9, 2, 35, 58, 44, 9, -8, -4, 1, 0,
    1, -2, -9, 1, 34, 58, 46, 11, -8, -5, 1, 0,
    1, -1, -8, -1, 31, 57, 47, 13, -7, -5, 1, 0,
    1, -1, -8, -2, 29, 56, 49, 15, -7, -6, 1, 1,
    1, 0, -8, -3, 26, 55, 51, 17, -7, -6, 1, 1,
    1, 0, -7, -4, 24, 54, 52, 19, -6, -7, 1, 1,
    1, 0, -7, -5, 22, 53, 53, 22, -5, -7, 0, 1,
    1, 1, -7, -6, 19, 52, 54, 24, -4, -7, 0, 1,
    1, 1, -6, -7, 17, 51, 55, 26, -3, -8, 0, 1,
    1, 1, -6, -7, 15, 49, 56, 29, -2, -8, -1, 1,
    0, 1, -5, -7, 13, 47, 57, 31, -1, -8, -1, 1,
    0, 1, -5, -8, 11, 46, 58, 34, 1, -9, -2, 1,
    0, 1, -4, -8, 9, 44, 58, 35, 2, -9, -2, 2,
    0, 1, -4, -9, 7, 43, 58, 38, 4, -9, -3, 2,
  },
  // 2 < ratio <= 5/2
  {
    -2, -7, 0, 17, 35, 43, 35, 17, 0, -7, -5, 2,
    -2, -7, -1, 16, 34, 43, 36, 18, 1, -7, -5, 2,
    -1, -7, -1, 14, 33, 43, 36, 19, 1, -6, -5, 2,
    -1, -7, -2, 13, 32, 42, 37, 20, 3, -6, -5, 2,
    0, -7, -3, 12, 31, 42, 38, 21, 3, -6, -5, 2,
    0, -7, -3, 11, 30, 42, 39, 23, 4, -6, -6, 1,
    0, -7, -4, 10, 29, 42, 40, 24, 5, -6, -6, 1,
    1, -7, -4, 9, 27, 41, 40, 25, 6, -5, -6, 1,
    1, -6, -5, 7, 26, 41, 41, 26, 7, -5, -6, 1,
    1, -6, -5, 6, 25, 40, 41, 27, 9, -4, -7, 1,
    1, -6, -6, 5, 24, 40, 42, 29, 10, -4, -7, 0,
    1, -6, -6, 4, 23, 39, 42, 30, 11, -3, -7, 0,
    2, -5, -6, 3, 21, 38, 42, 31, 12, -3, -7, 0,
    2, -5, -6, 3, 20, 37, 42, 32, 13, -2, -7, -1,
    2, -5, -6, 1, 19, 36, 43, 33, 14, -1, -7, -1,
    2, -5, -7, 1, 18, 36, 43, 34, 16, -1, -7, -2,
  },
  // 5/2 < ratio <= 20/7
  {
    -6, -3, 5, 19, 31, 36, 31, 19, 5, -3, -6, 0,
    -6, -4, 4, 18, 31, 37, 32, 20, 6, -3, -6, -1,
    -6, -4, 4, 17, 30, 36, 33, 21, 7, -3, -6, -1,
    -5, -5, 3, 16, 30, 36, 33, 22, 8, -2, -6, -2,
    -5, -5, 2, 15, 29, 36, 34, 23, 9, -2, -6, -2,
    -5, -5, 2, 15, 28, 36, 34, 24, 10, -2, -6, -3,
    -4, -5, 1, 14, 27, 36, 35, 24, 10, -1, -6, -3,
    -4, -5, 0, 13, 26, 35, 35, 25, 11, 0, -5, -3,
    -4, -6, 0, 12, 26, 36, 36, 26, 12, 0, -6, -4,
    -3, -5, 0, 11, 25, 35, 35, 26, 13, 0, -5, -4,
    -3, -6, -1, 10, 24, 35, 36, 27, 14, 1, -5, -4,
    -3, -6, -2, 10, 24, 34, 36, 28, 15, 2, -5, -5,
    -2, -6, -2, 9, 23, 34, 36, 29, 15, 2, -5, -5,
    -2, -6, -2, 8, 22, 33, 36, 30, 16, 3, -5, -5,
    -1, -6, -3, 7, 21, 33, 36, 30, 17, 4, -4, -6,
    -1, -6, -3, 6, 20, 32, 37, 31, 18, 4, -4, -6,
  },
  // 20/7 < ratio <= 15/4
  {
    -9, 0, 9, 20, 28, 32, 28, 20, 9, 0, -9, 0,
    -9, 0, 8, 19, 28, 32, 29, 20, 10, 0, -4, -5,
    -9, -1, 8, 18, 28, 32, 29, 21, 10, 1, -4, -5,
    -9, -1, 7, 18, 27, 32, 30, 22, 11, 1, -4, -6,
    -8, -2, 6, 17, 27, 32, 30, 22, 12, 2, -4, -6,
    -8, -2, 6, 16, 26, 32, 31, 23, 12, 2, -4, -6,
    -8, -2, 5, 16, 26, 31, 31, 23, 13, 3, -3, -7,
    -8, -3, 5, 15, 25, 31, 31, 24, 14, 4, -3, -7,
    -7, -3, 4, 14, 25, 31, 31, 25, 14, 4, -3, -7,
    -7, -3, 4, 14, 24, 31, 31, 25, 15, 5, -3, -8,
    -7, -3, 3, 13, 23, 31, 31, 26, 16, 5, -2, -8,
    -6, -4, 2, 12, 23, 31, 32, 26, 16, 6, -2, -8,
    -6, -4, 2, 12, 22, 30, 32, 27, 17, 6, -2, -8,
    -6, -4, 1, 11, 22, 30, 32, 27, 18, 7, -1, -9,
    -5, -4, 1, 10, 21, 29, 32, 28, 18, 8, -1, -9,
    -5, -4, 0, 10, 20, 29, 32, 28, 19, 8, 0, -9,
  },
  // ratio > 15/4
  {
    -8, 7, 13, 18, 22, 24, 22, 18, 13, 7, 2, -10,
    -8, 7, 13, 18, 22, 23, 22, 19, 13, 7, 2, -10,
    -8, 6, 12, 18, 22, 23, 22, 19, 14, 8, 2, -10,
    -9, 6, 12, 17, 22, 23, 23, 19, 14, 8, 3, -10,
    -9, 6, 12, 17, 21, 23, 23, 19, 14, 9, 3, -10,
    -9, 5, 11, 17, 21, 23, 23, 20, 15, 9, 3, -10,
    -9, 5, 11, 16, 21, 23, 23, 20, 15, 9, 4, -10,
    -9, 5, 10, 16, 21, 23, 23, 20, 15, 10, 4, -10,
    -10, 5, 10, 16, 20, 23, 23, 20, 16, 10, 5, -10,
    -10, 4, 10, 15, 20, 23, 23, 21, 16, 10, 5, -9,
    -10, 4, 9, 15, 20, 23, 23, 21, 16, 11, 5, -9,
    -10, 3, 9, 15, 20, 23, 23, 21, 17, 11, 5, -9,
    -10, 3, 9, 14, 19, 23, 23, 21, 17, 12, 6, -9,
    -10, 3, 8, 14, 19, 23, 23, 22, 17, 12, 6, -9,
    -10, 2, 8, 14, 19, 22, 23, 22, 18, 12, 6, -8,
    -10, 2, 7, 13, 19, 22, 23, 22, 18, 13, 7, -8,
  }
};

static const int lumaUpFilter1D[16 * 8] = {
  0, 0, 0, 64, 0, 0, 0, 0,
  0, 1, -3, 63, 4, -2, 1, 0,
  -1, 2, -5, 62, 8, -3, 1, 0,
  -1, 3, -8, 60, 13, -4, 1, 0,
  -1, 4, -10, 58, 17, -5, 1, 0,
  -1, 4, -11, 52, 26, -8, 3, -1,
  -1, 3, -9, 47, 31, -10, 4, -1,
  -1, 4, -11, 45, 34, -10, 4, -1,
  -1, 4, -11, 40, 40, -11, 4, -1,
  -1, 4, -10, 34, 45, -11, 4, -1,
  -1, 4, -10, 31, 47, -9, 3, -1,
  -1, 3, -8, 26, 52, -11, 4, -1,
  0, 1, -5, 17, 58, -10, 4, -1,
  0, 1, -4, 13, 60, -8, 3, -1,
  0, 1, -3, 8, 62, -5, 2, -1,
  0, 1, -2, 4, 63, -3, 1, 0
};

static const int chromaUpFilter1D[16 * 4] = {
  0, 64, 0, 0,
  -2, 62, 4, 0,
  -2, 58, 10, -2,
  -4, 56, 14, -2,
  -4, 54, 16, -2,
  -6, 52, 20, -2,
  -6, 46, 28, -4,
  -4, 42, 30, -4,
  -4, 36, 36, -4,
  -4, 30, 42, -4,
  -4, 28, 46, -6,
  -2, 20, 52, -6,
  -2, 16, 54, -4,
  -2, 14, 56, -4,
  -2, 10, 58, -2,
  0, 4, 62, -2
};

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
 *filter = downFilter[filter_ind][phase];
 return (sizeof(downFilter[0][0]) / sizeof(downFilter[0][0][0]));
}

//Helper function for choosing the correct filter
//Returns the size of the filter and the filter param is set to the correct filter
static int prepareFilter(const int** const filter, int is_upsampling, int is_luma, int filter_ind)
{
  if (is_upsampling) {
    //Upsampling so use 8- or 4-tap filters
    if (is_luma) {
      *filter = lumaUpFilter1D;
      return sizeof(lumaUpFilter[0]) / sizeof(lumaUpFilter[0][0]);
    }

    *filter = chromaUpFilter1D;
    return sizeof(chromaUpFilter[0]) / sizeof(chromaUpFilter[0][0]);
  }

  //Downsampling so use 12-tap filter
  *filter = downFilter1D[filter_ind];
  return (sizeof(downFilter[0][0]) / sizeof(downFilter[0][0][0]));
}

#define getFilterCoeff(filter, stride, phase, ind) ((filter)[(ind)+(phase)*(stride)])

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
//   pointer = clip_avx2(ref_pos, src_width, size, adder);
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
//    pointer = clip_avx2(ref_pos, src_width, size, upscaling_adder);
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
//   pointer = clip_avx2(ref_pos, src_height, size, adder);
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
//    pointer = clip_avx2(ref_pos, src_height, size, upscaling_adder);
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

//Resampling is done here per buffer
//void resample_avx2(const pic_buffer_t* const buffer, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma)
//{
// //TODO: Add cropping etc.
//
// //Choose best filter to use when downsampling
// //Need to use rounded values (to the closest multiple of 2,4,16 etc.)?
// int ver_filter = 0;
// int hor_filter = 0;
//
// int src_width = param->src_width + param->src_padding_x;
// int src_height = param->src_height + param->src_padding_y;
// int trgt_width = param->rnd_trgt_width;
// int trgt_height = param->rnd_trgt_height;
//
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
// // Horizontal resampling
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
//   pointer = clip_avx2(ref_pos, src_width, size, adder);
//   pointer = _mm256_permutevar8x32_epi32(pointer, order);
//
//   min = src_width-1;
//   smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
//   smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
//   min = _mm_extract_epi16(smallest_epi16, 0);
//   
//   tmp_row[j] = 0;
//   temp_mem = _mm256_load_si256((__m256i*)&(src_row[min]));
//   
//   decrese = _mm256_set1_epi32(min);
//
//   pointer = _mm256_sub_epi32(pointer, decrese);
//
//   temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
//
//   temp_filter = _mm256_load_si256((__m256i*)&(filter[0]));
//   temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);
//   
//
//   
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
//    pointer = clip_avx2(ref_pos, src_width, size, upscaling_adder);
//    pointer = _mm256_permutevar8x32_epi32(pointer, order);
//
//    smallest_epi16 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_packus_epi32(pointer, pointer), B11011000));
//    smallest_epi16 = _mm_minpos_epu16(smallest_epi16);
//    min = _mm_extract_epi16(smallest_epi16, 0);
//
//    temp_filter = _mm256_load_si256((__m256i* )&(filter[8]));
//    temp_mem = _mm256_load_si256((__m256i*)&(src_row[min]));
//
//    decrese = _mm256_set1_epi32(min);
//    pointer = _mm256_sub_epi32(pointer, decrese);
//    temp_mem = _mm256_permutevar8x32_epi32(temp_mem, pointer);
//
//    temp_mem = _mm256_mullo_epi32(temp_mem, temp_filter);   
//    
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
//  // Vertical resampling
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
//   
//   pointer = clip_avx2(ref_pos, src_height, size, adder);
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
//    pointer = clip_avx2(ref_pos, src_height, size, upscaling_adder);
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
#define SCALER_MIN(x,y) (((x) < (y)) ? (x) : (y))
#define SCALER_MAX(x,y) (((x) > (y)) ? (x) : (y))
#define SCALER_CLIP(val, mn, mx) (SCALER_MIN(mx,SCALER_MAX(mn,val)))
void resampleBlockStep_avx2(const pic_buffer_t* const src_buffer, const pic_buffer_t *const trgt_buffer, const int src_offset, const int trgt_offset, const int block_x, const int block_y, const int block_width, const int block_height, const scaling_parameter_t* const param, const int is_upscaling, const int is_luma, const int is_vertical)
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

  //Do resampling (vertical/horizontal) of the specified block into trgt_buffer using src_buffer
  for (int y = block_y; y < (block_y + block_height); y++) {

    pic_data_t* src = is_vertical ? src_buffer->data : &src_buffer->data[y * src_buffer->width];
    pic_data_t* trgt_row = &trgt_buffer->data[y * trgt_buffer->width];

    //Outer loop:
    //  if is_vertical -> loop over k (filter inds)
    //  if !is_vertical -> loop over x (block width)
    for (int o_ind = outer_init; o_ind < outer_bound; o_ind++) {

      const int t_ind = is_vertical ? y : o_ind; //trgt_buffer row/col index for cur resampling dir

                                                 //Inner loop:
                                                 //  if is_vertical -> loop over x (block width)
                                                 //  if !is_vertical -> loop over k (filter inds)-
      for (int i_ind = inner_init; i_ind < inner_bound; i_ind++) {

        const int f_ind = is_vertical ? o_ind : i_ind; //Filter index
        const int t_col = is_vertical ? i_ind : o_ind; //trgt_buffer column

                                                       //Calculate reference position in src pic
        int ref_pos_16 = (int)((unsigned int)(t_ind * scale + add) >> shift) - delta;
        int phase = ref_pos_16 & 15;
        int ref_pos = ref_pos_16 >> 4;

        //Choose filter
        //const int *filter;
        //const int f_size = getFilter(&filter, is_upscaling, is_luma, phase, filter_phase);

        //Set trgt buffer val to zero on first loop over filter
        if (f_ind == 0) {
          trgt_row[t_col + trgt_offset] = 0;
        }

        const int s_ind = SCALER_CLIP(ref_pos + f_ind - (filter_size >> 1) + 1, 0, src_size - 1); //src_buffer row/col index for cur resampling dir

                                                                                                  //Move src pointer to correct position (correct column in vertical resampling)
        pic_data_t *src_col = src + (is_vertical ? i_ind : 0);
        trgt_row[t_col + trgt_offset] += getFilterCoeff(filter, filter_size, phase, f_ind) * src_col[s_ind * s_stride + src_offset];

        //Scale values in trgt buffer to the correct range. Only done in the final loop over o_ind (block width)
        if (is_vertical && o_ind == outer_bound - 1) {
          trgt_row[t_col + trgt_offset] = SCALER_CLIP(is_upscaling ? (trgt_row[t_col + trgt_offset] + 2048) >> 12 : (trgt_row[t_col + trgt_offset] + 8192) >> 14, 0, 255);
        }
      }
    }
  }
}

//Set the default resample function
resample_block_step_func *const kvz_default_block_step_resample_func_avx2 = &DEFAULT_RESAMPLE_BLOCK_STEP_FUNC_AVX2;