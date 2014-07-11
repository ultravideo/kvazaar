#include "strategies-nal.h"

#include "generic/nal-generic.h"

void (*array_checksum)(const pixel* data,
                       const int height, const int width,
                       const int stride,
                       unsigned char checksum_out[SEI_HASH_MAX_LENGTH]);


int strategy_register_nal(void* opaque) {
  bool success = true;

  success &= strategy_register_nal_generic(opaque);
  
  return success;
}
