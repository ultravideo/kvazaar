#include "strategies/strategies-dct.h"
#include "strategyselector.h"

// Define function pointers.
dct_func * fast_forward_dst_4x4 = 0;

dct_func * dct_4x4 = 0;
dct_func * dct_8x8 = 0;
dct_func * dct_16x16 = 0;
dct_func * dct_32x32 = 0;

dct_func * fast_inverse_dst_4x4 = 0;

dct_func * idct_4x4 = 0;
dct_func * idct_8x8= 0;
dct_func * idct_16x16 = 0;
dct_func * idct_32x32 = 0;


// Headers for platform optimizations.
#include "generic/dct-generic.h"
#include "avx2/dct-avx2.h"


int strategy_register_dct(void* opaque) {
  bool success = true;

  success &= strategy_register_dct_generic(opaque);

  if (g_hardware_flags.intel_flags.avx2) {
    success &= strategy_register_dct_avx2(opaque);
  }

  return success;
}


/**
* \brief  Get a function that calculates SAD for NxN block.
*
* \param n  Width of the region for which SAD is calculated.
*
* \returns  Pointer to cost_16bit_nxn_func.
*/
dct_func * get_dct_func(int8_t width, int32_t mode)
{
  switch (width) {
  case 4:
    switch (mode){
    case 65535:
      return dct_4x4;
    default:
      return fast_forward_dst_4x4;
  }
  case 8:
    return dct_8x8;
  case 16:
    return dct_16x16;
  case 32:
    return dct_32x32;
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
dct_func * get_idct_func(int8_t width, int32_t mode)
{
  switch (width) {
  case 4:
    switch (mode){
    case 65535:
      return idct_4x4;
    default:
      return fast_inverse_dst_4x4;
    }
  case 8:
    return idct_8x8;
  case 16:
    return idct_16x16;
  case 32:
    return idct_32x32;
  default:
    return NULL;
  }
}
