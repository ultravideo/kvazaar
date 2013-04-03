/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file picture.c
    \brief Functions to handle pictures
    \author Marko Viitanen
    \date 2012-06
    
  This file contains all the needed functions to handle pictures

*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "picture.h"

/** \defgroup picture_group Picture handler group
 *  This group contains all picture related stuff
 *  @{
 */


/*!
    \brief Allocate memory for picture_list
    \param size initial array size
    \return picture_list pointer, NULL on failure
*/
  picture_list *picture_list_init(int size)
  {
    picture_list *list = (picture_list *)malloc(sizeof(picture_list));
    list->size = size;
    if(size > 0)
    {
      list->pics = (picture**)malloc(sizeof(picture*)*size);
    }

    list->used_size = 0;

    return list;
  }

  /*!
    \brief Resize picture_list array
    \param list picture_list pointer
    \param size new array size
    \return 1 on success, 0 on failure
  */
  int picture_list_resize(picture_list *list, int size)
  {
    unsigned int i;
    picture** old_pics = NULL;

    //No need to do anything when resizing to same size
    if(size == list->size)
    {
      return 1;
    }

    //Save the old list
    if(list->used_size > 0)
    {
      old_pics = list->pics;
    }

    //allocate space for the new list
    list->pics = (picture**)malloc(sizeof(picture*)*size);

    //Copy everthing from the old list to the new if needed.
    if(old_pics != NULL)
    {
      for(i = 0; i < list->used_size; i++)
      {
        list->pics[i] = old_pics[i];
      }

      free(old_pics);
    }

    return 1;
  }

  /*!
    \brief Free memory allocated to the picture_list
    \param list picture_list pointer
    \return 1 on success, 0 on failure
  */
  int picture_list_destroy(picture_list *list)
  {
    unsigned int i;
    if(list->used_size > 0)
    {
      for(i = 0; i < list->used_size; i++)
      {
        picture_destroy(list->pics[i]);
      }
    }
    
    if(list->size > 0)
    {
      free(list->pics);
    }
    free(list);
    return 1;
  }

  /*!
    \brief Free memory allocated to picture
    \param pic picture pointer
    \return 1 on success, 0 on failure
  */
  int picture_destroy(picture *pic)
  {
    free(pic->uData);
    free(pic->vData);
    free(pic->yData);
    return 1;
  }


  /** @} */ // end of group1


#include <math.h>
#define PSNRMAX (255.0*255.0)

//Calculates image PSNR value
double imagePSNR(uint8_t *frame1, uint8_t *frame2, uint32_t x, uint32_t y)
{   
  uint64_t MSE=0;
  uint64_t MSEtemp=0;
  double psnr=0.0;
  int32_t index;

  //Calculate MSE
  for(index=x*y-1;index>=0;index--)
  {
    MSEtemp=frame1[index]-frame2[index];
    MSE+=MSEtemp*MSEtemp;
  }
  MSE/=x*y;

  //Avoid division by zero
  if(MSE==0) return 99.0;

  //The PSNR
  psnr=10*log10(PSNRMAX/MSE);

  //Thats it.
  return psnr;
}

//Sum of Absolute Difference for block
uint32_t SAD(uint8_t *block,uint8_t* block2, uint32_t x, uint32_t y)
{
  uint32_t i;
  uint32_t sum=0;
  for(i=0;i<x*y;i+=4)
  {
    sum+=abs((int32_t)block[i]-(int32_t)block2[i]);
    sum+=abs((int32_t)block[i+1]-(int32_t)block2[i+1]);
    sum+=abs((int32_t)block[i+2]-(int32_t)block2[i+2]);
    sum+=abs((int32_t)block[i+3]-(int32_t)block2[i+3]);
  }

  return sum;    
}

uint32_t SAD64x64(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2)
{
  int32_t i,ii,y,x;
  uint32_t sum=0;
  for(y=0;y<64;y++)
  {
    i = y*stride1; 
    ii = y*stride2;
    for(x = 0; x < 64;x++)
    {
      sum+=abs((int32_t)block[i+x]-(int32_t)block2[ii+x]);
    }

  }

  return sum;    
}

uint32_t SAD32x32(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2)
{
  int32_t i,ii,y;
  uint32_t sum=0;
  for(y=0;y<32;y++)
  {
    i = y*stride1; 
    ii = y*stride2;
    sum+=abs((int32_t)block[i]-(int32_t)block2[ii]);
    sum+=abs((int32_t)block[i+1]-(int32_t)block2[ii+1]);
    sum+=abs((int32_t)block[i+2]-(int32_t)block2[ii+2]);
    sum+=abs((int32_t)block[i+3]-(int32_t)block2[ii+3]);
    sum+=abs((int32_t)block[i+4]-(int32_t)block2[ii+4]);
    sum+=abs((int32_t)block[i+5]-(int32_t)block2[ii+5]);
    sum+=abs((int32_t)block[i+6]-(int32_t)block2[ii+6]);
    sum+=abs((int32_t)block[i+7]-(int32_t)block2[ii+7]);
    sum+=abs((int32_t)block[i+8]-(int32_t)block2[ii+8]);
    sum+=abs((int32_t)block[i+9]-(int32_t)block2[ii+9]);
    sum+=abs((int32_t)block[i+10]-(int32_t)block2[ii+10]);
    sum+=abs((int32_t)block[i+11]-(int32_t)block2[ii+11]);
    sum+=abs((int32_t)block[i+12]-(int32_t)block2[ii+12]);
    sum+=abs((int32_t)block[i+13]-(int32_t)block2[ii+13]);
    sum+=abs((int32_t)block[i+14]-(int32_t)block2[ii+14]);
    sum+=abs((int32_t)block[i+15]-(int32_t)block2[ii+15]);
    sum+=abs((int32_t)block[i+16]-(int32_t)block2[ii+16]);
    sum+=abs((int32_t)block[i+17]-(int32_t)block2[ii+17]);
    sum+=abs((int32_t)block[i+18]-(int32_t)block2[ii+18]);
    sum+=abs((int32_t)block[i+19]-(int32_t)block2[ii+19]);
    sum+=abs((int32_t)block[i+20]-(int32_t)block2[ii+20]);
    sum+=abs((int32_t)block[i+21]-(int32_t)block2[ii+21]);
    sum+=abs((int32_t)block[i+22]-(int32_t)block2[ii+22]);
    sum+=abs((int32_t)block[i+23]-(int32_t)block2[ii+23]);
    sum+=abs((int32_t)block[i+24]-(int32_t)block2[ii+24]);
    sum+=abs((int32_t)block[i+25]-(int32_t)block2[ii+25]);
    sum+=abs((int32_t)block[i+26]-(int32_t)block2[ii+26]);
    sum+=abs((int32_t)block[i+27]-(int32_t)block2[ii+27]);
    sum+=abs((int32_t)block[i+28]-(int32_t)block2[ii+28]);
    sum+=abs((int32_t)block[i+29]-(int32_t)block2[ii+29]);
    sum+=abs((int32_t)block[i+30]-(int32_t)block2[ii+30]);
    sum+=abs((int32_t)block[i+31]-(int32_t)block2[ii+31]);
  }

  return sum;    
}


uint32_t SAD16x16(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2)
{
  int32_t i,ii,y;
  uint32_t sum=0;
  for(y=0;y<16;y++)
  {
    i = y*stride1; 
    ii = y*stride2;
    sum+=abs((int32_t)block[i]-(int32_t)block2[ii]);
    sum+=abs((int32_t)block[i+1]-(int32_t)block2[ii+1]);
    sum+=abs((int32_t)block[i+2]-(int32_t)block2[ii+2]);
    sum+=abs((int32_t)block[i+3]-(int32_t)block2[ii+3]);
    sum+=abs((int32_t)block[i+4]-(int32_t)block2[ii+4]);
    sum+=abs((int32_t)block[i+5]-(int32_t)block2[ii+5]);
    sum+=abs((int32_t)block[i+6]-(int32_t)block2[ii+6]);
    sum+=abs((int32_t)block[i+7]-(int32_t)block2[ii+7]);
    sum+=abs((int32_t)block[i+8]-(int32_t)block2[ii+8]);
    sum+=abs((int32_t)block[i+9]-(int32_t)block2[ii+9]);
    sum+=abs((int32_t)block[i+10]-(int32_t)block2[ii+10]);
    sum+=abs((int32_t)block[i+11]-(int32_t)block2[ii+11]);
    sum+=abs((int32_t)block[i+12]-(int32_t)block2[ii+12]);
    sum+=abs((int32_t)block[i+13]-(int32_t)block2[ii+13]);
    sum+=abs((int32_t)block[i+14]-(int32_t)block2[ii+14]);
    sum+=abs((int32_t)block[i+15]-(int32_t)block2[ii+15]);
  }

  return sum;    
}


uint32_t SAD8x8(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2)
{
  int32_t i,ii,y;
  uint32_t sum=0;
  for(y=0;y<8;y++)
  {
    i = y*stride1; 
    ii = y*stride2;
    sum+=abs((int32_t)block[i]-(int32_t)block2[ii]);
    sum+=abs((int32_t)block[i+1]-(int32_t)block2[ii+1]);
    sum+=abs((int32_t)block[i+2]-(int32_t)block2[ii+2]);
    sum+=abs((int32_t)block[i+3]-(int32_t)block2[ii+3]);
    sum+=abs((int32_t)block[i+4]-(int32_t)block2[ii+4]);
    sum+=abs((int32_t)block[i+5]-(int32_t)block2[ii+5]);
    sum+=abs((int32_t)block[i+6]-(int32_t)block2[ii+6]);
    sum+=abs((int32_t)block[i+7]-(int32_t)block2[ii+7]);
  }

  return sum;    
}

uint32_t SAD4x4(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2)
{
  int32_t i,ii,y;
  uint32_t sum=0;
  for(y=0;y<4;y++)
  {
    i = y*stride1; 
    ii = y*stride2;
    sum+=abs((int32_t)block[i]-(int32_t)block2[ii]);
    sum+=abs((int32_t)block[i+1]-(int32_t)block2[ii+1]);
    sum+=abs((int32_t)block[i+2]-(int32_t)block2[ii+2]);
    sum+=abs((int32_t)block[i+3]-(int32_t)block2[ii+3]);
  }

  return sum;
}