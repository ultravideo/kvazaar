#include "strategies-ipol.h"
#include "strategyselector.h"

// Define function pointers.
ipol_func *filter_inter_quarterpel_luma;
ipol_func *filter_inter_halfpel_chroma;
ipol_func *filter_inter_octpel_chroma;
epol_func *extend_borders;

// Headers for platform optimizations.
#include "generic/ipol-generic.h"
#include "avx2/ipol-avx2.h"


int strategy_register_ipol(void* opaque) {
  bool success = true;

  success &= strategy_register_ipol_generic(opaque);

  if (g_hardware_flags.intel_flags.avx2) {
    success &= strategy_register_ipol_avx2(opaque);
  }
  return success;
}