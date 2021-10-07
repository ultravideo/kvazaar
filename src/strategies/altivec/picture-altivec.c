/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "strategies/altivec/picture-altivec.h"

#if COMPILE_POWERPC_ALTIVEC
#undef bool
#include <altivec.h>
#define bool _Bool
#include <stdlib.h>

#include "kvazaar.h"
#include "strategyselector.h"


static unsigned reg_sad_altivec(const kvz_pixel * const data1, const kvz_pixel * const data2,
                        const int width, const int height, const unsigned stride1, const unsigned stride2)
{
  vector unsigned int vsad = {0,0,0,0}, vzero = {0,0,0,0}; 
  vector signed int sumdiffs;
  int tmpsad, sad = 0;
  
  int y, x;
  
  for (y = 0; y < height; ++y) {
    vector unsigned char perm1, perm2;
    
    perm1 = vec_lvsl(0, &data1[y * stride1]);
    perm2 = vec_lvsl(0, &data2[y * stride2]);
    
    for (x = 0; x <= width-16; x+=16) {
      vector unsigned char t1, t2, t3, t4, t5;
      vector unsigned char *current, *previous;
      
      current = (vector unsigned char *) &data1[y * stride1 + x];
      previous = (vector unsigned char *) &data2[y * stride2 + x];
      
      t1  = vec_perm(current[0], current[1], perm1 );  /* align current vector  */ 
      t2  = vec_perm(previous[0], previous[1], perm2 );/* align previous vector */ 
      t3  = vec_max(t1, t2 );      /* find largest of two           */ 
      t4  = vec_min(t1, t2 );      /* find smaller of two           */ 
      t5  = vec_sub(t3, t4);       /* find absolute difference      */ 
      vsad = vec_sum4s(t5, vsad);    /* accumulate sum of differences */
    }

    for (; x < width; ++x) {
      sad += abs(data1[y * stride1 + x] - data2[y * stride2 + x]);
    }
  }
  
  sumdiffs = vec_sums((vector signed int) vsad, (vector signed int) vzero);
  /* copy vector sum into unaligned result */
  sumdiffs = vec_splat( sumdiffs, 3);
  vec_ste( sumdiffs, 0, &tmpsad );
  sad += tmpsad;
  
  return sad;
}

#endif //COMPILE_POWERPC_ALTIVEC

int kvz_strategy_register_picture_altivec(void* opaque, uint8_t bitdepth) {
  bool success = true;
#if COMPILE_POWERPC_ALTIVEC
  if(bitdepth == 8){
    success &= kvz_strategyselector_register(opaque, "reg_sad", "altivec", 10, &reg_sad_altivec);
  }
#endif
  return success;
}
