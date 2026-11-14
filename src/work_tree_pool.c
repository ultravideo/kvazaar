/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2026, project contributors
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


#include "work_tree_pool.h"

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

//! Alignment of the work tree and of each of its buffers.
#define WORK_TREE_ALIGNMENT 64
#define WORK_TREE_ALIGN_UP(size) \
  (((size) + WORK_TREE_ALIGNMENT - 1) & ~((size_t)WORK_TREE_ALIGNMENT - 1))

//! Number of slots in one word of the bit mask.
#define MASK_BITS 64

struct kvz_work_tree_pool {
  //! Protects used.
  pthread_mutex_t lock;
  enum kvz_chroma_format chroma_format;
  int num_slots;
  //! Bit mask of the slots that are in use, MASK_BITS slots per word.
  uint64_t *used;
  //! The work trees of the slots.
  lcu_t **slots;
};


/**
 * \brief Allocate a work tree with buffers for the given chroma format.
 *
 * The lcu_t structs and all of their buffers are in a single allocation.
 * The buffers are not initialized, the search does that for every LCU.
 *
 * \return the work tree (KVZ_WORK_TREE_SIZE lcu_t structs), or NULL
 */
static lcu_t * work_tree_alloc(const enum kvz_chroma_format csp)
{
  const size_t chroma_size     = KVZ_LCU_CHROMA_SIZE(csp);
  const size_t chroma_rec_size = KVZ_LCU_CHROMA_REC_SIZE(csp);

  const size_t lcus_bytes         = WORK_TREE_ALIGN_UP(KVZ_WORK_TREE_SIZE * sizeof(lcu_t));
  const size_t luma_px_bytes      = WORK_TREE_ALIGN_UP(LCU_LUMA_SIZE * sizeof(kvz_pixel));
  const size_t chroma_px_bytes    = WORK_TREE_ALIGN_UP(chroma_size * sizeof(kvz_pixel));
  const size_t chroma_rec_bytes   = WORK_TREE_ALIGN_UP(chroma_rec_size * sizeof(kvz_pixel));
  const size_t luma_coeff_bytes   = WORK_TREE_ALIGN_UP(LCU_LUMA_SIZE * sizeof(coeff_t));
  const size_t chroma_coeff_bytes = WORK_TREE_ALIGN_UP(chroma_size * sizeof(coeff_t));

  // One ref buffer shared by all the depths, and rec and coeff buffers for
  // every depth.
  const size_t ref_bytes   = luma_px_bytes + 2 * chroma_px_bytes;
  const size_t depth_bytes = luma_px_bytes + 2 * chroma_rec_bytes +
                             luma_coeff_bytes + 2 * chroma_coeff_bytes;
  const size_t total_bytes = lcus_bytes + ref_bytes + KVZ_WORK_TREE_SIZE * depth_bytes;

  // Reserve room for the alignment and for a pointer to the start of the
  // allocation right before the work tree.
  void *const memory = malloc(total_bytes + WORK_TREE_ALIGNMENT + sizeof(void*));
  if (!memory) return NULL;

  uint8_t *ptr = (uint8_t*)WORK_TREE_ALIGN_UP((uintptr_t)memory + sizeof(void*));
  ((void**)ptr)[-1] = memory;

  lcu_t *const work_tree = (lcu_t*)ptr;
  ptr += lcus_bytes;

  const bool has_chroma = chroma_size > 0;
  lcu_yuv_t ref;
  ref.chroma_format = csp;
  ref.y = (kvz_pixel*)ptr;                              ptr += luma_px_bytes;
  ref.u = has_chroma ? (kvz_pixel*)ptr : NULL;          ptr += chroma_px_bytes;
  ref.v = has_chroma ? (kvz_pixel*)ptr : NULL;          ptr += chroma_px_bytes;

  for (int depth = 0; depth < KVZ_WORK_TREE_SIZE; ++depth) {
    lcu_t *const lcu = &work_tree[depth];
    lcu->ref = ref;

    lcu->rec.chroma_format = csp;
    lcu->rec.y = (kvz_pixel*)ptr;                       ptr += luma_px_bytes;
    lcu->rec.u = has_chroma ? (kvz_pixel*)ptr : NULL;   ptr += chroma_rec_bytes;
    lcu->rec.v = has_chroma ? (kvz_pixel*)ptr : NULL;   ptr += chroma_rec_bytes;

    lcu->coeff.y = (coeff_t*)ptr;                       ptr += luma_coeff_bytes;
    lcu->coeff.u = has_chroma ? (coeff_t*)ptr : NULL;   ptr += chroma_coeff_bytes;
    lcu->coeff.v = has_chroma ? (coeff_t*)ptr : NULL;   ptr += chroma_coeff_bytes;
  }

  return work_tree;
}


static void work_tree_free(lcu_t *work_tree)
{
  if (work_tree) {
    free(((void**)work_tree)[-1]);
  }
}


/**
 * \brief Allocate a pool of work trees.
 *
 * \param chroma_format   chroma format of the encoded video
 * \param num_slots       maximum number of LCU searches running at the same
 *                        time, i.e. the number of threads
 * \return the pool, or NULL on failure
 */
kvz_work_tree_pool_t * kvz_work_tree_pool_alloc(const enum kvz_chroma_format chroma_format,
                                                const int num_slots)
{
  kvz_work_tree_pool_t *pool = calloc(1, sizeof(kvz_work_tree_pool_t));
  if (!pool) return NULL;

  pool->chroma_format = chroma_format;
  pool->num_slots = MAX(1, num_slots);

  const int num_words = (pool->num_slots + MASK_BITS - 1) / MASK_BITS;
  pool->used  = calloc(num_words, sizeof(uint64_t));
  pool->slots = calloc(pool->num_slots, sizeof(lcu_t*));
  if (!pool->used || !pool->slots || pthread_mutex_init(&pool->lock, NULL) != 0) {
    free(pool->used);
    free(pool->slots);
    free(pool);
    return NULL;
  }

  for (int i = 0; i < pool->num_slots; ++i) {
    pool->slots[i] = work_tree_alloc(chroma_format);
    if (!pool->slots[i]) {
      kvz_work_tree_pool_free(pool);
      return NULL;
    }
  }

  return pool;
}


void kvz_work_tree_pool_free(kvz_work_tree_pool_t *pool)
{
  if (!pool) return;

  for (int i = 0; i < pool->num_slots; ++i) {
    work_tree_free(pool->slots[i]);
  }
  pthread_mutex_destroy(&pool->lock);
  free(pool->used);
  free(pool->slots);
  free(pool);
}


/**
 * \brief Take a work tree from the pool.
 *
 * The contents of the work tree are undefined. If all of the slots are in
 * use, a new work tree is allocated, which kvz_work_tree_pool_put frees.
 *
 * \return KVZ_WORK_TREE_SIZE lcu_t structs, or NULL if allocation fails
 */
lcu_t * kvz_work_tree_pool_get(kvz_work_tree_pool_t *pool)
{
  lcu_t *work_tree = NULL;

  pthread_mutex_lock(&pool->lock);
  for (int i = 0; i < pool->num_slots; ++i) {
    const uint64_t bit = (uint64_t)1 << (i % MASK_BITS);
    if (!(pool->used[i / MASK_BITS] & bit)) {
      pool->used[i / MASK_BITS] |= bit;
      work_tree = pool->slots[i];
      break;
    }
  }
  pthread_mutex_unlock(&pool->lock);

  if (!work_tree) {
    // More searches running at the same time than there are slots. This
    // should not happen, but it is not an error either.
    work_tree = work_tree_alloc(pool->chroma_format);
  }
  return work_tree;
}


/**
 * \brief Return a work tree taken with kvz_work_tree_pool_get to the pool.
 */
void kvz_work_tree_pool_put(kvz_work_tree_pool_t *pool, lcu_t *work_tree)
{
  if (!work_tree) return;

  pthread_mutex_lock(&pool->lock);
  for (int i = 0; i < pool->num_slots; ++i) {
    if (pool->slots[i] == work_tree) {
      assert(pool->used[i / MASK_BITS] & ((uint64_t)1 << (i % MASK_BITS)));
      pool->used[i / MASK_BITS] &= ~((uint64_t)1 << (i % MASK_BITS));
      pthread_mutex_unlock(&pool->lock);
      return;
    }
  }
  pthread_mutex_unlock(&pool->lock);

  // Not one of the slots, so it was allocated by kvz_work_tree_pool_get.
  work_tree_free(work_tree);
}
