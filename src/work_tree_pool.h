#ifndef WORK_TREE_POOL_H_
#define WORK_TREE_POOL_H_
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

/**
 * \ingroup Compression
 * \file
 * \brief Thread safe pool of search work trees.
 *
 * The search of an LCU needs MAX_PU_DEPTH + 1 lcu_t structs (the work
 * tree). Their pixel and coefficient buffers are allocated for the exact
 * chroma format in use, so they are much smaller for 4:2:0 than they
 * would be if they were sized for 4:4:4. Allocating them for every LCU
 * would be slow, so the encoder keeps a pool of work trees with one slot
 * for each thread that may run an LCU search at the same time. A bit mask
 * tells which slots are in use.
 *
 * The buffers of a work tree are set up once, when the pool is allocated:
 * - Every depth has its own reconstruction and coefficient buffers.
 * - The source pixels (ref) are not modified by the search, so all of the
 *   depths share one ref buffer.
 */

#include "global.h" // IWYU pragma: keep
#include "cu.h"
#include "kvazaar.h"

typedef struct kvz_work_tree_pool kvz_work_tree_pool_t;

//! Number of lcu_t structs in a work tree.
#define KVZ_WORK_TREE_SIZE (MAX_PU_DEPTH + 1)

kvz_work_tree_pool_t * kvz_work_tree_pool_alloc(enum kvz_chroma_format chroma_format,
                                                int num_slots);
void kvz_work_tree_pool_free(kvz_work_tree_pool_t *pool);

lcu_t * kvz_work_tree_pool_get(kvz_work_tree_pool_t *pool);
void kvz_work_tree_pool_put(kvz_work_tree_pool_t *pool, lcu_t *work_tree);

#endif // WORK_TREE_POOL_H_
