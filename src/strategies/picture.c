#include "picture-generic.c"
#if COMPILE_INTEL_SSE2
#include "picture-sse2.c"
#endif
#if COMPILE_INTEL_SSE2 && COMPILE_INTEL_SSE41
#include "picture-sse41.c"
#endif

unsigned (*reg_sad)(const pixel * const data1, const pixel * const data2,
                    const int width, const int height, const unsigned stride1, const unsigned stride2);


static int strategy_register_picture(void* opaque) {
  if (!strategy_register_picture_generic(opaque)) return 0;
  
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
  return 1;
}
