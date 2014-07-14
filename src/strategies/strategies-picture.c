#include "strategies-picture.h"
#include "strategyselector.h"

// Define function pointers.
reg_sad_func * reg_sad = 0;

cost_pixel_nxn_func * sad_8bit_4x4 = 0;
cost_pixel_nxn_func * sad_8bit_8x8 = 0;
cost_pixel_nxn_func * sad_8bit_16x16 = 0;
cost_pixel_nxn_func * sad_8bit_32x32 = 0;
cost_pixel_nxn_func * sad_8bit_64x64 = 0;

cost_pixel_nxn_func * satd_8bit_4x4 = 0;
cost_pixel_nxn_func * satd_8bit_8x8 = 0;
cost_pixel_nxn_func * satd_8bit_16x16 = 0;
cost_pixel_nxn_func * satd_8bit_32x32 = 0;
cost_pixel_nxn_func * satd_8bit_64x64 = 0;


// Headers for platform optimizations.
#include "generic/picture-generic.h"
#include "sse2/picture-sse2.h"
#include "sse41/picture-sse41.h"
#include "altivec/picture-altivec.h"
#include "picture/picture-avx.c"


int strategy_register_picture(void* opaque) {
  bool success = true;

  success &= strategy_register_picture_generic(opaque);

  if (g_hardware_flags.intel_flags.sse2) {
    success &= strategy_register_picture_sse2(opaque);
    if (g_hardware_flags.intel_flags.avx) {
    success &= strategy_register_picture_avx(opaque);
    }
  }
  if (g_hardware_flags.intel_flags.sse41) {
    success &= strategy_register_picture_sse41(opaque);
  }
  if (g_hardware_flags.powerpc_flags.altivec) {
    success &= strategy_register_picture_altivec(opaque);
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
cost_pixel_nxn_func * pixels_get_satd_func(unsigned n)
{
  switch (n) {
  case 4:
    return satd_8bit_4x4;
  case 8:
    return satd_8bit_8x8;
  case 16:
    return satd_8bit_16x16;
  case 32:
    return satd_8bit_32x32;
  case 64:
    return satd_8bit_64x64;
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
cost_pixel_nxn_func * pixels_get_sad_func(unsigned n)
{
  switch (n) {
  case 4:
    return sad_8bit_4x4;
  case 8:
    return sad_8bit_8x8;
  case 16:
    return sad_8bit_16x16;
  case 32:
    return sad_8bit_32x32;
  case 64:
    return sad_8bit_64x64;
  default:
    return NULL;
  }
}
