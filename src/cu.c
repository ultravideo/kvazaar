/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include <string.h>
#include <stdlib.h>

#include "cu.h"
#include "threads.h"


/**
 * \brief Number of PUs in a CU.
 *
 * Indexed by part_mode_t values.
 */
const uint8_t kvz_part_mode_num_parts[] = {
  1, // 2Nx2N
  2, // 2NxN
  2, // Nx2N
  4, // NxN
  2, // 2NxnU
  2, // 2NxnD
  2, // nLx2N
  2, // nRx2N
};

/**
 * \brief PU offsets.
 *
 * Indexed by [part mode][PU number][axis].
 *
 * Units are 1/4 of the width of the CU.
 */
const uint8_t kvz_part_mode_offsets[][4][2] = {
  { {0, 0}                         }, // 2Nx2N
  { {0, 0}, {0, 2}                 }, // 2NxN
  { {0, 0}, {2, 0}                 }, // Nx2N
  { {0, 0}, {2, 0}, {0, 2}, {2, 2} }, // NxN
  { {0, 0}, {0, 1}                 }, // 2NxnU
  { {0, 0}, {0, 3}                 }, // 2NxnD
  { {0, 0}, {1, 0}                 }, // nLx2N
  { {0, 0}, {3, 0}                 }, // nRx2N
};

/**
 * \brief PU sizes.
 *
 * Indexed by [part mode][PU number][axis].
 *
 * Units are 1/4 of the width of the CU.
 */
const uint8_t kvz_part_mode_sizes[][4][2] = {
  { {4, 4}                         }, // 2Nx2N
  { {4, 2}, {4, 2}                 }, // 2NxN
  { {2, 4}, {2, 4}                 }, // Nx2N
  { {2, 2}, {2, 2}, {2, 2}, {2, 2} }, // NxN
  { {4, 1}, {4, 3}                 }, // 2NxnU
  { {4, 3}, {4, 1}                 }, // 2NxnD
  { {1, 4}, {3, 4}                 }, // nLx2N
  { {3, 4}, {1, 4}                 }, // nRx2N
};


#define BLIT_COEFF_CASE(n) case n:\
  for (y = 0; y < n; ++y) {\
    memcpy(&dst[y*dst_stride], &orig[y*orig_stride], n * sizeof(coeff_t));\
  }\
  break;

void kvz_coefficients_blit(const coeff_t * const orig, coeff_t * const dst,
                         const unsigned width, const unsigned height,
                         const unsigned orig_stride, const unsigned dst_stride)
{
  unsigned y;
  
  int nxn_width = (width == height) ? width : 0;
  switch (nxn_width) {
    BLIT_COEFF_CASE(4)
    BLIT_COEFF_CASE(8)
    BLIT_COEFF_CASE(16)
    BLIT_COEFF_CASE(32)
    BLIT_COEFF_CASE(64)
  default:
    for (y = 0; y < height; ++y) {
      memcpy(&dst[y*dst_stride], &orig[y*orig_stride], width * sizeof(coeff_t));
    }
    break;
  }
}

unsigned kvz_coefficients_calc_abs(const coeff_t *const buf, const int buf_stride,
                        const int width)
{
  int sum = 0;
  int y, x;

  for (y = 0; y < width; ++y) {
    for (x = 0; x < width; ++x) {
      sum += abs(buf[x + y * buf_stride]);
    }
  }

  return sum;
}

/**
 * \brief Allocate a CU array.
 *
 * \param width   width of the array in luma pixels
 * \param height  height of the array in luma pixels
 */
cu_array_t * kvz_cu_array_alloc(const int width, const int height) {
  cu_array_t *cua = MALLOC(cu_array_t, 1);

  const int width_scu  = (width  + 7) >> 3;
  const int height_scu = (height + 7) >> 3;
  assert(width_scu * 8 >= width);
  assert(height_scu * 8 >= height);
  const unsigned cu_array_size = width_scu * height_scu;
  cua->data = calloc(cu_array_size, sizeof(cu_info_t));
  cua->width = width_scu << 3;
  cua->height = height_scu << 3;
  cua->refcount = 1;

  return cua;
}


int kvz_cu_array_free(cu_array_t * const cua)
{
  int32_t new_refcount;
  if (!cua) return 1;

  new_refcount = KVZ_ATOMIC_DEC(&(cua->refcount));
  //Still we have some references, do nothing
  if (new_refcount > 0) return 1;

  FREE_POINTER(cua->data);
  free(cua);

  return 1;
}


/**
 * \brief Copy part of a cu array to another cu array.
 *
 * All values are in luma pixels.
 *
 * \param dst     destination array
 * \param dst_x   x-coordinate of the left edge of the copied area in dst
 * \param dst_y   y-coordinate of the top edge of the copied area in dst
 * \param src     source array
 * \param src_x   x-coordinate of the left edge of the copied area in src
 * \param src_y   y-coordinate of the top edge of the copied area in src
 * \param width   width of the area to copy
 * \param height  height of the area to copy
 */
void kvz_cu_array_copy(cu_array_t* dst,       int dst_x, int dst_y,
                       const cu_array_t* src, int src_x, int src_y,
                       int width, int height)
{
  // Convert values from pixel coordinates to array indices.
  int src_stride = src->width >> 3;
  int dst_stride = dst->width >> 3;
  const cu_info_t* src_ptr = &src->data[(src_x >> 3) + (src_y >> 3) * src_stride];
  cu_info_t* dst_ptr       = &dst->data[(dst_x >> 3) + (dst_y >> 3) * dst_stride];

  // Number of bytes to copy per row.
  const size_t row_size = sizeof(cu_info_t) * (width >> 3);

  width = MIN(width,   MIN(src->width  - src_x, dst->width  - dst_x));
  height = MIN(height, MIN(src->height - src_y, dst->height - dst_y));

  assert(src_x + width  <= src->width);
  assert(src_y + height <= src->height);
  assert(dst_x + width  <= dst->width);
  assert(dst_y + height <= dst->height);

  for (int i = 0; i < (height >> 3); ++i) {
    memcpy(dst_ptr, src_ptr, row_size);
    src_ptr += src_stride;
    dst_ptr += dst_stride;
  }
}
