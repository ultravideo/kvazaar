#include "strategies-picture.h"

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


// Include inline functions.
#include "picture/picture-generic.c"
#if COMPILE_INTEL_SSE2
#include "picture/picture-sse2.c"
#endif
#if COMPILE_INTEL_SSE2 && COMPILE_INTEL_SSE41
#include "picture/picture-sse41.c"
#endif
#if COMPILE_INTEL_AVX
#include "picture/picture-avx.c"
#endif
#if COMPILE_POWERPC_ALTIVEC
#include "picture/picture-altivec.c"
#endif


int strategy_register_picture(void* opaque) {
  if (!strategy_register_picture_generic(opaque)) return 0;
  
#if COMPILE_INTEL
  if (g_hardware_flags.intel_flags.sse2) {
#if COMPILE_INTEL_SSE2
    if (!strategy_register_picture_sse2(opaque)) return 0;
#endif
    if (g_hardware_flags.intel_flags.sse41) {
#if COMPILE_INTEL_SSE2 && COMPILE_INTEL_SSE41
      if (!strategy_register_picture_sse41(opaque)) return 0;
#endif
    }
    if (g_hardware_flags.intel_flags.avx) {
#if COMPILE_INTEL_AVX
      if (!strategy_register_picture_avx(opaque)) return 0;
#endif
    }
  }
#endif //COMPILE_INTEL

#if COMPILE_POWERPC
  if (g_hardware_flags.powerpc_flags.altivec) {
#if COMPILE_POWERPC_ALTIVEC
    if (!strategy_register_picture_altivec(opaque)) return 0;
#endif //COMPILE_POWERPC_ALTIVEC
  }
#endif //COMPILE_POWERPC
  return 1;
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
