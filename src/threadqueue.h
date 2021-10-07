#ifndef THREADQUEUE_H_
#define THREADQUEUE_H_
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
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Threading
 * \file
 * Container for worker tasks.
 */

#include "global.h" // IWYU pragma: keep

#include <pthread.h>

typedef struct threadqueue_job_t threadqueue_job_t;
typedef struct threadqueue_queue_t threadqueue_queue_t;

threadqueue_queue_t * kvz_threadqueue_init(int thread_count);

threadqueue_job_t * kvz_threadqueue_job_create(void (*fptr)(void *arg), void *arg);
int kvz_threadqueue_submit(threadqueue_queue_t * threadqueue, threadqueue_job_t *job);

int kvz_threadqueue_job_dep_add(threadqueue_job_t *job, threadqueue_job_t *dependency);

threadqueue_job_t *kvz_threadqueue_copy_ref(threadqueue_job_t *job);

void kvz_threadqueue_free_job(threadqueue_job_t **job_ptr);

int kvz_threadqueue_waitfor(threadqueue_queue_t * threadqueue, threadqueue_job_t * job);
int kvz_threadqueue_stop(threadqueue_queue_t * threadqueue);
void kvz_threadqueue_free(threadqueue_queue_t * threadqueue);

#endif // THREADQUEUE_H_
