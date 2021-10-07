/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
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


cu_info_t* kvz_cu_array_at(cu_array_t *cua, unsigned x_px, unsigned y_px)
{
  return (cu_info_t*) kvz_cu_array_at_const(cua, x_px, y_px);
}


const cu_info_t* kvz_cu_array_at_const(const cu_array_t *cua, unsigned x_px, unsigned y_px)
{
  assert(x_px < cua->width);
  assert(y_px < cua->height);
  return &(cua)->data[(x_px >> 2) + (y_px >> 2) * ((cua)->stride >> 2)];
}


/**
 * \brief Allocate a CU array.
 *
 * \param width   width of the array in luma pixels
 * \param height  height of the array in luma pixels
 */
cu_array_t * kvz_cu_array_alloc(const int width, const int height)
{
  cu_array_t *cua = MALLOC(cu_array_t, 1);

  // Round up to a multiple of LCU width and divide by cell width.
  const int width_scu  = CEILDIV(width,  LCU_WIDTH) * LCU_WIDTH / SCU_WIDTH;
  const int height_scu = CEILDIV(height, LCU_WIDTH) * LCU_WIDTH / SCU_WIDTH;
  const unsigned cu_array_size = width_scu * height_scu;

  cua->base     = NULL;
  cua->data     = calloc(cu_array_size, sizeof(cu_info_t));
  cua->width    = width_scu  * SCU_WIDTH;
  cua->height   = height_scu * SCU_WIDTH;
  cua->stride   = cua->width;
  cua->refcount = 1;

  return cua;
}


cu_array_t * kvz_cu_subarray(cu_array_t *base,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height)
{
  assert(x_offset + width <= base->width);
  assert(y_offset + height <= base->height);

  if (x_offset == 0 &&
      y_offset == 0 &&
      width == base->width &&
      height == base->height)
  {
    return kvz_cu_array_copy_ref(base);
  }

  cu_array_t *cua = MALLOC(cu_array_t, 1);

  // Find the real base array.
  cu_array_t *real_base = base;
  while (real_base->base) {
    real_base = real_base->base;
  }
  cua->base     = kvz_cu_array_copy_ref(real_base);
  cua->data     = kvz_cu_array_at(base, x_offset, y_offset);
  cua->width    = width;
  cua->height   = height;
  cua->stride   = base->stride;
  cua->refcount = 1;

  return cua;
}

void kvz_cu_array_free(cu_array_t **cua_ptr)
{
  cu_array_t *cua = *cua_ptr;
  if (cua == NULL) return;
  *cua_ptr = NULL;

  int new_refcount = KVZ_ATOMIC_DEC(&cua->refcount);
  if (new_refcount > 0) {
    // Still we have some references, do nothing.
    return;
  }

  assert(new_refcount == 0);

  if (!cua->base) {
    FREE_POINTER(cua->data);
  } else {
    kvz_cu_array_free(&cua->base);
    cua->data = NULL;
  }

  FREE_POINTER(cua);
}


/**
 * \brief Get a new pointer to a cu array.
 *
 * Increment reference count and return the cu array.
 */
cu_array_t * kvz_cu_array_copy_ref(cu_array_t* cua)
{
  int32_t new_refcount = KVZ_ATOMIC_INC(&cua->refcount);
  // The caller should have had another reference and we added one
  // reference so refcount should be at least 2.
  assert(new_refcount >= 2);
  return cua;
}


/**
 * \brief Copy an lcu to a cu array.
 *
 * All values are in luma pixels.
 *
 * \param dst     destination array
 * \param dst_x   x-coordinate of the left edge of the copied area in dst
 * \param dst_y   y-coordinate of the top edge of the copied area in dst
 * \param src     source lcu
 */
void kvz_cu_array_copy_from_lcu(cu_array_t* dst, int dst_x, int dst_y, const lcu_t *src)
{
  const int dst_stride = dst->stride >> 2;
  for (int y = 0; y < LCU_WIDTH; y += SCU_WIDTH) {
    for (int x = 0; x < LCU_WIDTH; x += SCU_WIDTH) {
      const cu_info_t *from_cu = LCU_GET_CU_AT_PX(src, x, y);
      const int x_scu = (dst_x + x) >> 2;
      const int y_scu = (dst_y + y) >> 2;
      cu_info_t *to_cu = &dst->data[x_scu + y_scu * dst_stride];
      memcpy(to_cu,                  from_cu, sizeof(*to_cu));
    }
  }
}
