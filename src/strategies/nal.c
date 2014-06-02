#include "nal.h"
#include "nal-generic.c"

void (*array_checksum)(const pixel* data,
                       const int height, const int width,
                       const int stride,
                       unsigned char checksum_out[SEI_HASH_MAX_LENGTH]);


static int strategy_register_nal(void* opaque) {
  if (!strategy_register_nal_generic(opaque)) return 0;
  
  return 1;
}
