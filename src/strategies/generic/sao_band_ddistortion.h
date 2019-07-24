#ifndef SAO_BAND_DDISTORTION_H_
#define SAO_BAND_DDISTORTION_H_

// #include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "sao.h"

static int sao_band_ddistortion_generic(const encoder_state_t * const state,
                                        const kvz_pixel *orig_data,
                                        const kvz_pixel *rec_data,
                                        int block_width,
                                        int block_height,
                                        int band_pos,
                                        const int sao_bands[4])
{
  int y, x;
  int shift = state->encoder_control->bitdepth-5;
  int sum = 0;
  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      const int32_t curr_pos = y * block_width + x;

      kvz_pixel rec  =  rec_data[curr_pos];
      kvz_pixel orig = orig_data[curr_pos];

      int32_t band = (rec >> shift) - band_pos;
      int32_t offset = 0;
      if (band >= 0 && band <= 3) {
        offset = sao_bands[band];
      }
      // Offset is applied to reconstruction, so it is subtracted from diff.

      int32_t diff  = orig - rec;
      int32_t delta = diff - offset;

      int32_t dmask = (offset == 0) ? -1 : 0;
      diff  &= ~dmask;
      delta &= ~dmask;

      sum += delta * delta - diff * diff;
    }
  }

  return sum;
}

#endif
