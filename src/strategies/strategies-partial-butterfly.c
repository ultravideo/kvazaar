#include "strategies/strategies-partial-butterfly.h"
#include "strategyselector.h"

// Define function pointers.
partial_butterfly_func * partial_butterfly_4 = 0;
partial_butterfly_func * partial_butterfly_8 = 0;
partial_butterfly_func * partial_butterfly_16 = 0;
partial_butterfly_func * partial_butterfly_32 = 0;

partial_butterfly_func * partial_butterfly_inverse_4 = 0;
partial_butterfly_func * partial_butterfly_inverse_8 = 0;
partial_butterfly_func * partial_butterfly_inverse_16 = 0;
partial_butterfly_func * partial_butterfly_inverse_32 = 0;


// Headers for platform optimizations.
#include "generic/partial-butterfly-generic.h"
#include "avx2/partial-butterfly-avx2.h"


int strategy_register_partial_butterfly(void* opaque) {
  bool success = true;

  success &= strategy_register_partial_butterfly_generic(opaque);

  if (g_hardware_flags.intel_flags.avx2) {
    success &= strategy_register_partial_butterfly_avx2(opaque);
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
partial_butterfly_func * get_partial_butterfly_func(unsigned n)
{
  switch (n) {
  case 4:
    return partial_butterfly_4;
  case 8:
    return partial_butterfly_8;
  case 16:
    return partial_butterfly_16;
  case 32:
    return partial_butterfly_32;
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
partial_butterfly_func * get_partial_butterfly_inverse_func(unsigned n)
{
  switch (n) {
  case 4:
    return partial_butterfly_inverse_4;
  case 8:
    return partial_butterfly_inverse_8;
  case 16:
    return partial_butterfly_inverse_16;
  case 32:
    return partial_butterfly_inverse_32;
  default:
    return NULL;
  }
}
