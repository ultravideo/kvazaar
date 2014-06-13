#include "strategies-picture.h"

#include "picture/picture-generic.c"
#if COMPILE_INTEL_SSE2
#include "picture/picture-sse2.c"
#endif
#if COMPILE_INTEL_SSE2 && COMPILE_INTEL_SSE41
#include "picture/picture-sse41.c"
#endif
#if COMPILE_POWERPC_ALTIVEC
#include "picture/picture-altivec.c"
#endif


reg_sad_func reg_sad;


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
