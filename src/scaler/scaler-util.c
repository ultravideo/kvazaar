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

#include "scaler-util.h"

//Used for clipping values
int kvz_clip_scaler(int val, int min, int max)
{
  if (val <= min)
    return min;
  if (val >= max)
    return max;

  return val;
}

//Helper function for choosing the correct filter
//Returns the size of the filter and the filter param is set to the correct filter
int kvz_getFilter(const int** const filter, int is_upsampling, int is_luma, int phase, int filter_ind)
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
unsigned kvz_prepareFilter(const int** const filter, int is_upsampling, int is_luma, int filter_ind)
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


#define FILTER_SELECT_TEMPLATE(postfix,...) \
  if (is_upsampling) {\
    if (is_luma) {\
      *filter = (void *)lumaUpFilter1D##postfix;\
      return sizeof(lumaUpFilter[0]) / sizeof(lumaUpFilter[0][0]);\
    }\
    *filter = (void *)chromaUpFilter1D##postfix;\
    return sizeof(chromaUpFilter[0]) / sizeof(chromaUpFilter[0][0]);\
  }\
  *filter = (void *)downFilter1D##postfix[filter_ind];\
  return (sizeof(downFilter[0][0]) / sizeof(downFilter[0][0][0]));\

//Helper function for choosing the correct filter with the given depth
//Returns the size of the filter and the filter param is set to the correct filter
unsigned kvz_prepareFilterDepth(const void** const filter, int is_upsampling, int is_luma, int filter_ind, int depth)
{
  switch (depth)
  {
    case sizeof(int) :
    {
      FILTER_SELECT_TEMPLATE(,);
    }

    case sizeof(short) :
    {
      FILTER_SELECT_TEMPLATE(_16bit);
    }

    case sizeof(char) :
    {
      FILTER_SELECT_TEMPLATE(_8bit);
    }
  }

  //No valid depth
  *filter = (void *)0;
  return 0;
}