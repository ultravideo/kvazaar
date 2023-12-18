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

#include "strategies/strategies-picture.h"

#include "strategies/altivec/picture-altivec.h"
#include "strategies/avx2/picture-avx2.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/sse2/picture-sse2.h"
#include "strategies/sse41/picture-sse41.h"
#include "strategyselector.h"


// Define function pointers.
reg_sad_func * kvz_reg_sad = 0;

cost_pixel_nxn_func * kvz_sad_4x4 = 0;
cost_pixel_nxn_func * kvz_sad_8x8 = 0;
cost_pixel_nxn_func * kvz_sad_16x16 = 0;
cost_pixel_nxn_func * kvz_sad_32x32 = 0;
cost_pixel_nxn_func * kvz_sad_64x64 = 0;

cost_pixel_nxn_func * kvz_satd_4x4 = 0;
cost_pixel_nxn_func * kvz_satd_8x8 = 0;
cost_pixel_nxn_func * kvz_satd_16x16 = 0;
cost_pixel_nxn_func * kvz_satd_32x32 = 0;
cost_pixel_nxn_func * kvz_satd_64x64 = 0;

cost_pixel_nxn_multi_func * kvz_sad_4x4_dual = 0;
cost_pixel_nxn_multi_func * kvz_sad_8x8_dual = 0;
cost_pixel_nxn_multi_func * kvz_sad_16x16_dual = 0;
cost_pixel_nxn_multi_func * kvz_sad_32x32_dual = 0;
cost_pixel_nxn_multi_func * kvz_sad_64x64_dual = 0;

cost_pixel_nxn_multi_func * kvz_satd_4x4_dual = 0;
cost_pixel_nxn_multi_func * kvz_satd_8x8_dual = 0;
cost_pixel_nxn_multi_func * kvz_satd_16x16_dual = 0;
cost_pixel_nxn_multi_func * kvz_satd_32x32_dual = 0;
cost_pixel_nxn_multi_func * kvz_satd_64x64_dual = 0;

cost_pixel_any_size_func * kvz_satd_any_size = 0;
cost_pixel_any_size_multi_func * kvz_satd_any_size_quad = 0;

pixels_calc_ssd_func * kvz_pixels_calc_ssd = 0;

inter_recon_bipred_func * kvz_bipred_average = 0;

get_optimized_sad_func *kvz_get_optimized_sad = 0;
ver_sad_func *kvz_ver_sad = 0;
hor_sad_func *kvz_hor_sad = 0;

pixel_var_func *kvz_pixel_var = 0;


int kvz_strategy_register_picture(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= kvz_strategy_register_picture_generic(opaque, bitdepth);

  if (kvz_g_hardware_flags.intel_flags.sse2) {
    success &= kvz_strategy_register_picture_sse2(opaque, bitdepth);
  }
  if (kvz_g_hardware_flags.intel_flags.sse41) {
    success &= kvz_strategy_register_picture_sse41(opaque, bitdepth);
  }
  if (kvz_g_hardware_flags.intel_flags.avx2) {
    success &= kvz_strategy_register_picture_avx2(opaque, bitdepth);
  }
  if (kvz_g_hardware_flags.powerpc_flags.altivec) {
    success &= kvz_strategy_register_picture_altivec(opaque, bitdepth);
  }

  return success;
}


/**
* \brief  Get a function that calculates SATD for NxN block.
*
* \param n  Width of the region for which SATD is calculated.
*
* \returns  Pointer to cost_16bit_nxn_func.
*/
cost_pixel_nxn_func * kvz_pixels_get_satd_func(unsigned n)
{
  switch (n) {
  case 4:
    return kvz_satd_4x4;
  case 8:
    return kvz_satd_8x8;
  case 16:
    return kvz_satd_16x16;
  case 32:
    return kvz_satd_32x32;
  case 64:
    return kvz_satd_64x64;
  default:
    return NULL;
  }
}


/**
* \brief  Get a function that calculates SAD for NxN block.
*
* \param n  Width of the region for which SAD is calculated.
*
* \returns  Pointer to cost_16bit_nxn_func.
*/
cost_pixel_nxn_func * kvz_pixels_get_sad_func(unsigned n)
{
  switch (n) {
  case 4:
    return kvz_sad_4x4;
  case 8:
    return kvz_sad_8x8;
  case 16:
    return kvz_sad_16x16;
  case 32:
    return kvz_sad_32x32;
  case 64:
    return kvz_sad_64x64;
  default:
    return NULL;
  }
}

/**
* \brief  Get a function that calculates SATDs for 2 NxN blocks.
*
* \param n  Width of the region for which SATD is calculated.
*
* \returns  Pointer to cost_pixel_nxn_multi_func.
*/
cost_pixel_nxn_multi_func * kvz_pixels_get_satd_dual_func(unsigned n)
{
  switch (n) {
  case 4:
    return kvz_satd_4x4_dual;
  case 8:
    return kvz_satd_8x8_dual;
  case 16:
    return kvz_satd_16x16_dual;
  case 32:
    return kvz_satd_32x32_dual;
  case 64:
    return kvz_satd_64x64_dual;
  default:
    return NULL;
  }
}


/**
* \brief  Get a function that calculates SADs for 2 NxN blocks.
*
* \param n  Width of the region for which SAD is calculated.
*
* \returns  Pointer to cost_pixel_nxn_multi_func.
*/
cost_pixel_nxn_multi_func * kvz_pixels_get_sad_dual_func(unsigned n)
{
  switch (n) {
  case 4:
    return kvz_sad_4x4_dual;
  case 8:
    return kvz_sad_8x8_dual;
  case 16:
    return kvz_sad_16x16_dual;
  case 32:
    return kvz_sad_32x32_dual;
  case 64:
    return kvz_sad_64x64_dual;
  default:
    return NULL;
  }
}
