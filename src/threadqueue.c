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

#include "threadqueue.h"

#include <errno.h> // ETIMEDOUT
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "threads.h"


typedef struct {
  threadqueue_queue_t * threadqueue;
  int worker_id;
} threadqueue_worker_spec;

#define THREADQUEUE_LIST_REALLOC_SIZE 32

//#define PTHREAD_COND_SIGNAL(c) fprintf(stderr, "%s:%d pthread_cond_signal(%s=%p)\n", __FUNCTION__, __LINE__, #c, c); if (pthread_cond_signal((c)) != 0) { fprintf(stderr, "pthread_cond_signal(%s=%p) failed!\n", #c, c); assert(0); return 0; }
//#define PTHREAD_COND_BROADCAST(c) fprintf(stderr, "%s:%d pthread_cond_broadcast(%s=%p)\n", __FUNCTION__, __LINE__, #c, c); if (pthread_cond_broadcast((c)) != 0) { fprintf(stderr, "pthread_cond_broadcast(%s=%p) failed!\n", #c, c); assert(0); return 0; }
//#define PTHREAD_COND_WAIT(c,l) fprintf(stderr, "%s:%d pthread_cond_wait(%s=%p, %s=%p)\n", __FUNCTION__, __LINE__, #c, c, #l, l); if (pthread_cond_wait((c),(l)) != 0) { fprintf(stderr, "pthread_cond_wait(%s=%p, %s=%p) failed!\n", #c, c, #l, l); assert(0); return 0; } else {fprintf(stderr, "%s:%d pthread_cond_wait(%s=%p, %s=%p) (done)\n", __FUNCTION__, __LINE__, #c, c, #l, l);}
//#define PTHREAD_LOCK(l) fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p) (try)\n", __FUNCTION__, __LINE__, #l, l); if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s=%p) failed!\n", #l, l); assert(0); return 0; } else {fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}
//#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s=%p) failed!\n", #l, l); assert(0); return 0; }  else {fprintf(stderr, "%s:%d pthread_mutex_unlock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}


#define PTHREAD_COND_SIGNAL(c) if (pthread_cond_signal((c)) != 0) { fprintf(stderr, "pthread_cond_signal(%s=%p) failed!\n", #c, c); assert(0); return 0; }
#define PTHREAD_COND_BROADCAST(c) if (pthread_cond_broadcast((c)) != 0) { fprintf(stderr, "pthread_cond_broadcast(%s=%p) failed!\n", #c, c); assert(0); return 0; }

#ifndef _PTHREAD_DUMP
#define PTHREAD_COND_WAIT(c,l) if (pthread_cond_wait((c),(l)) != 0) { fprintf(stderr, "pthread_cond_wait(%s=%p, %s=%p) failed!\n", #c, c, #l, l); assert(0); return 0; }
#define PTHREAD_LOCK(l) if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s) failed!\n", #l); assert(0); return 0; }
#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s) failed!\n", #l); assert(0); return 0; }

#else  //PTHREAD_DUMP
#define PTHREAD_LOCK(l) do { \
  PERFORMANCE_MEASURE_START(); \
  if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s) failed!\n", #l); assert(0); return 0; } \
  PERFORMANCE_MEASURE_END(NULL, "pthread_mutex_lock(%s=%p)@%s:%d",#l,l,__FUNCTION__, __LINE__); \
} while (0);

#define PTHREAD_UNLOCK(l) do { \
  PERFORMANCE_MEASURE_START(); \
  if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s) failed!\n", #l); assert(0); return 0; } \
  PERFORMANCE_MEASURE_END(NULL, "pthread_mutex_unlock(%s=%p)@%s:%d",#l,l,__FUNCTION__, __LINE__); \
} while (0);

#define PTHREAD_COND_WAIT(c,l) do { \
  PERFORMANCE_MEASURE_START(); \
  if (pthread_cond_wait((c),(l)) != 0) { fprintf(stderr, "pthread_cond_wait(%s=%p, %s=%p) failed!\n", #c, c, #l, l); assert(0); return 0;} \
  PERFORMANCE_MEASURE_END(NULL, "pthread_cond_wait(%s=%p, %s=%p)@%s:%d",#c, c, #l, l,__FUNCTION__, __LINE__); \
} while (0);
#endif //PTHREAD_DUMP

static void* threadqueue_worker(void* threadqueue_worker_spec_opaque)
{
  threadqueue_worker_spec * const threadqueue_worker_spec = threadqueue_worker_spec_opaque;
  threadqueue_queue_t * const threadqueue = threadqueue_worker_spec->threadqueue;
  threadqueue_job_t * next_job = NULL;

#ifdef KVZ_DEBUG
  KVZ_GET_TIME(&threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id]);
#endif //KVZ_DEBUG

  for(;;) {
    threadqueue_job_t * job = NULL;

    PTHREAD_LOCK(&threadqueue->lock);

    while(!threadqueue->stop && threadqueue->queue_waiting_execution == 0 && !next_job) {
      // Wait until there is something to do in the queue.
      PTHREAD_COND_WAIT(&threadqueue->cond, &threadqueue->lock);
    }

    if(threadqueue->stop) {
      if (next_job) {
        // Put a job we had already reserved back into the queue.
        // FIXME: This lock should be unnecessary, as nobody else is allowed
        // to touch this job when it's running.
        PTHREAD_LOCK(&next_job->lock);
        next_job->state = THREADQUEUE_JOB_STATE_QUEUED;
        PTHREAD_UNLOCK(&next_job->lock);
      }
      break;
    }

    //Find a task (should be fast enough)
    job = NULL;
    if (next_job) {
      assert(next_job->ndepends == 0);
      job = next_job;
    } else {
      //FIXME: if not using OWF, the first is better than the second, otherwise we should use the second order
      //for (i = threadqueue->queue_count - 1; i >= threadqueue->queue_start; --i) {
      //for (i = threadqueue->queue_start; i < threadqueue->queue_count; ++i) {

      for (int i = (threadqueue->fifo ? threadqueue->queue_start : threadqueue->queue_count - 1);
           (threadqueue->fifo ? i < threadqueue->queue_count : i >= threadqueue->queue_start); 
           (threadqueue->fifo ? ++i : --i)) {
        threadqueue_job_t * const i_job = threadqueue->queue[i];

        if (i_job->state == THREADQUEUE_JOB_STATE_QUEUED && i_job->ndepends == 0) {
          // Once we found the job with no dependancies, lock it and change
          // its state to running, so nobody else can claim it.
          PTHREAD_LOCK(&i_job->lock);
          if (i_job->state == THREADQUEUE_JOB_STATE_QUEUED && i_job->ndepends == 0) {
            job = i_job;
            job->state = THREADQUEUE_JOB_STATE_RUNNING;
          }
          PTHREAD_UNLOCK(&i_job->lock);
          if (job) break;
        }
      }
    }

    if (!job) {
      // We have no job. Probably because more threads were woken up than
      // there were jobs to do.
      PTHREAD_UNLOCK(&threadqueue->lock);
    } else {
      // We have a job with ndepends==0 and its state is running.
      assert(job->state == THREADQUEUE_JOB_STATE_RUNNING);

      // Advance queue_start to skip all the running jobs.
      while (threadqueue->queue_start < threadqueue->queue_count &&
             threadqueue->queue[threadqueue->queue_start]->state != THREADQUEUE_JOB_STATE_QUEUED)
      {
        threadqueue->queue_start++;
      }
      
      if (!next_job) {
        --threadqueue->queue_waiting_execution;
        ++threadqueue->queue_running;
      }

      PTHREAD_UNLOCK(&threadqueue->lock);

#ifdef KVZ_DEBUG
      job->debug_worker_id = threadqueue_worker_spec->worker_id;
      KVZ_GET_TIME(&job->debug_clock_start);
#endif //KVZ_DEBUG

      job->fptr(job->arg);

#ifdef KVZ_DEBUG
      job->debug_worker_id = threadqueue_worker_spec->worker_id;
      KVZ_GET_TIME(&job->debug_clock_stop);
#endif //KVZ_DEBUG

      // FIXME: This lock should be unnecessary, as nobody else is allowed
      // to touch this job when it's running.
      PTHREAD_LOCK(&job->lock);
      assert(job->state == THREADQUEUE_JOB_STATE_RUNNING);

      job->state = THREADQUEUE_JOB_STATE_DONE;

      next_job = NULL;

      int queue_waiting_dependency_decr = 0;
      int queue_waiting_execution_incr = 0;

      // Go throught all the jobs that depend on this one, decresing their ndepends.
      for (int i = 0; i < job->rdepends_count; ++i) {
        threadqueue_job_t * const depjob = job->rdepends[i];
        // Note that we lock the dependency AFTER locking the source. This avoids a deadlock in dep_add.
        PTHREAD_LOCK(&depjob->lock);

        assert(depjob->state == THREADQUEUE_JOB_STATE_QUEUED);
        assert(depjob->ndepends > 0);
        --depjob->ndepends;

        // Count how many jobs can now start executing so we know how many
        // threads to wake up.
        if (depjob->ndepends == 0) {
          if (!next_job) {
            // Avoid having to find a new job for this worker through the
            // queue by taking one of the jobs that depended on current job.
            next_job = depjob;
            depjob->state = THREADQUEUE_JOB_STATE_RUNNING;
          } else {
            ++queue_waiting_execution_incr;
          }
          ++queue_waiting_dependency_decr;
        }

        PTHREAD_UNLOCK(&depjob->lock);
      }
      
      PTHREAD_UNLOCK(&job->lock);

      PTHREAD_LOCK(&threadqueue->lock);

      assert(threadqueue->queue_waiting_dependency >= queue_waiting_dependency_decr);

      // This thread will 
      if (!next_job) {
        // We didn't find a new job, so this thread will have to go wait.
        threadqueue->queue_running--;
      }
      threadqueue->queue_waiting_dependency -= queue_waiting_dependency_decr;
      threadqueue->queue_waiting_execution += queue_waiting_execution_incr;

      // Wake up enough threads to take care of the tasks now lacking dependancies.
      for (int i = 0; i < queue_waiting_execution_incr; ++i) {
        PTHREAD_COND_SIGNAL(&threadqueue->cond);
      }

      // Signal main thread that a job has been completed.
      pthread_cond_signal(&threadqueue->cb_cond);

      PTHREAD_UNLOCK(&threadqueue->lock);
    }
  }

  // We got out of the loop because threadqueue->stop == 1. The queue is locked.
  assert(threadqueue->stop);
  --threadqueue->threads_running;
  
#ifdef KVZ_DEBUG
  KVZ_GET_TIME(&threadqueue->debug_clock_thread_end[threadqueue_worker_spec->worker_id]);
  
  fprintf(threadqueue->debug_log, "\t%d\t-\t%lf\t+%lf\t-\tthread\n", threadqueue_worker_spec->worker_id, KVZ_CLOCK_T_AS_DOUBLE(threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id]), KVZ_CLOCK_T_DIFF(threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id], threadqueue->debug_clock_thread_end[threadqueue_worker_spec->worker_id]));
#endif //KVZ_DEBUG
  
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  free(threadqueue_worker_spec_opaque);
  
  pthread_exit(NULL);
  
  return NULL;
}


int kvz_threadqueue_init(threadqueue_queue_t * const threadqueue, int thread_count, int fifo) {
  int i;
  if (pthread_mutex_init(&threadqueue->lock, NULL) != 0) {
    fprintf(stderr, "pthread_mutex_init failed!\n");
    assert(0);
    return 0;
  }
  if (pthread_cond_init(&threadqueue->cond, NULL) != 0) {
    fprintf(stderr, "pthread_cond_init failed!\n");
    assert(0);
    return 0;
  }
  
  if (pthread_cond_init(&threadqueue->cb_cond, NULL) != 0) {
    fprintf(stderr, "pthread_cond_init failed!\n");
    assert(0);
    return 0;
  }
  
  threadqueue->stop = 0;
  threadqueue->fifo = !!fifo;
  threadqueue->threads_running = 0;
  threadqueue->threads_count = thread_count;
  
  threadqueue->threads = MALLOC(pthread_t, thread_count);
  if (!threadqueue->threads) {
    fprintf(stderr, "Could not malloc threadqueue->threads!\n");
    return 0;
  }
#ifdef KVZ_DEBUG
  threadqueue->debug_clock_thread_start = MALLOC(KVZ_CLOCK_T, thread_count);
  assert(threadqueue->debug_clock_thread_start);
  threadqueue->debug_clock_thread_end = MALLOC(KVZ_CLOCK_T, thread_count);
  assert(threadqueue->debug_clock_thread_end);
  threadqueue->debug_log = fopen("threadqueue.log", "w");
#endif //KVZ_DEBUG
  
  threadqueue->queue = NULL;
  threadqueue->queue_size = 0;
  threadqueue->queue_count = 0;
  threadqueue->queue_start = 0;
  threadqueue->queue_waiting_execution = 0;
  threadqueue->queue_waiting_dependency = 0;
  threadqueue->queue_running = 0;
  
  //Lock the queue before creating threads, to ensure they all have correct information
  PTHREAD_LOCK(&threadqueue->lock);
  
  for(i = 0; i < thread_count; i++) {
    threadqueue_worker_spec *tqws = MALLOC(threadqueue_worker_spec, 1);
    if (tqws) {
      tqws->threadqueue = threadqueue;
      tqws->worker_id = i;
      if(pthread_create(&(threadqueue->threads[i]), NULL, threadqueue_worker, (void*)tqws) != 0) {
          fprintf(stderr, "pthread_create failed!\n");
          assert(0);
          return 0;
      }
      threadqueue->threads_running++;
    } else {
      fprintf(stderr, "Could not allocate threadqueue_worker_spec structure!\n");
      PTHREAD_UNLOCK(&threadqueue->lock);
      return 0;
    }
  }
  
  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
}

/**
 * \brief Get a new pointer to a job.
 *
 * Increment reference count and return the job.
 */
threadqueue_job_t *kvz_threadqueue_copy_ref(threadqueue_job_t *job)
{
  // The caller should have had another reference.
  assert(job->refcount > 0);
  KVZ_ATOMIC_INC(&job->refcount);
  return job;
}


/**
 * \brief Free a job.
 *
 * Decrement reference count of the job. If no references exist any more,
 * deallocate associated memory and destroy mutexes.
 *
 * Sets the job pointer to NULL.
 */
void kvz_threadqueue_free_job(threadqueue_job_t **job_ptr)
{
  threadqueue_job_t *job = *job_ptr;
  if (job == NULL) return;
  *job_ptr = NULL;

  int new_refcount = KVZ_ATOMIC_DEC(&job->refcount);
  if (new_refcount > 0) {
    // There are still references so we don't free the data yet.
    return;
  }

  assert(new_refcount >= 0);

#ifdef KVZ_DEBUG
  FREE_POINTER(job->debug_description);
#endif

  FREE_POINTER(job->rdepends);
  pthread_mutex_destroy(&job->lock);
  FREE_POINTER(job);
}

static void threadqueue_free_jobs(threadqueue_queue_t * const threadqueue) {
  for (int i = 0; i < threadqueue->queue_count; ++i) {
    kvz_threadqueue_free_job(&threadqueue->queue[i]);
  }
  threadqueue->queue_count = 0;
  threadqueue->queue_start = 0;
}

int kvz_threadqueue_finalize(threadqueue_queue_t * const threadqueue) {
  //Flush the queue
  if (!kvz_threadqueue_flush(threadqueue)) {
    fprintf(stderr, "Unable to flush threadqueue!\n");
    return 0;
  }
  
  //Lock threadqueue
  PTHREAD_LOCK(&threadqueue->lock);
  
  //Free job memory
  threadqueue_free_jobs(threadqueue);
  
  if (threadqueue->stop) {
    fprintf(stderr, "threadqueue already stopping\n");
    
    if (pthread_mutex_unlock(&threadqueue->lock) != 0) {
      fprintf(stderr, "pthread_mutex_unlock failed!\n");
      assert(0);
      return 0;
    }
    assert(0); //We should get here...
    return 0;
  }
  
  threadqueue->stop = 1;
  
  if (pthread_cond_broadcast(&(threadqueue->cond)) != 0) {
    fprintf(stderr, "pthread_cond_broadcast failed!\n");
    PTHREAD_UNLOCK(&threadqueue->lock);
    assert(0);
    return 0;
  }
  //Unlock it now, since all jobs have to stpo
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  //Join threads
  for(int i = 0; i < threadqueue->threads_count; i++) {
    if(pthread_join(threadqueue->threads[i], NULL) != 0) {
      fprintf(stderr, "pthread_join failed!\n");
      return 0;
    }
  }
  
#ifdef KVZ_DEBUG
  FREE_POINTER(threadqueue->debug_clock_thread_start);
  FREE_POINTER(threadqueue->debug_clock_thread_end);
  fclose(threadqueue->debug_log);
#endif
  
  //Free allocated stuff
  FREE_POINTER(threadqueue->queue);
  threadqueue->queue_count = 0;
  threadqueue->queue_size = 0;
  threadqueue->queue_start = 0;
  
  FREE_POINTER(threadqueue->threads);
  threadqueue->threads_count = 0;
  
  if (pthread_mutex_destroy(&threadqueue->lock) != 0) {
    fprintf(stderr, "pthread_mutex_destroy failed!\n");
    assert(0);
    return 0;
  }
  if (pthread_cond_destroy(&threadqueue->cond) != 0) {
    fprintf(stderr, "pthread_cond_destroy failed!\n");
    assert(0);
    return 0;
  }
  
  if (pthread_cond_destroy(&threadqueue->cb_cond) != 0) {
    fprintf(stderr, "pthread_cond_destroy failed!\n");
    assert(0);
    return 0;
  }
  
  return 1;
}

int kvz_threadqueue_flush(threadqueue_queue_t * const threadqueue) {
  int notdone = 1;
  
  //Lock the queue
  PTHREAD_LOCK(&threadqueue->lock);
  
  do {
    notdone = threadqueue->queue_waiting_execution + threadqueue->queue_waiting_dependency + threadqueue->queue_running;

    if (notdone > 0) {
      PTHREAD_COND_BROADCAST(&(threadqueue->cond));
      PTHREAD_COND_WAIT(&threadqueue->cb_cond, &threadqueue->lock);
    }
  } while (notdone > 0);
  
  threadqueue_free_jobs(threadqueue);

  assert(threadqueue->queue_waiting_dependency == 0 && threadqueue->queue_waiting_execution == 0 && threadqueue->queue_running == 0);

  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
}

int kvz_threadqueue_waitfor(threadqueue_queue_t * const threadqueue, threadqueue_job_t * const job) {
  int job_done = 0;
  
  //NULL job is clearly OK :-)
  if (!job) return 1;
  
  //Lock the queue
  PTHREAD_LOCK(&threadqueue->lock);
  do {
    
    PTHREAD_LOCK(&job->lock);
    job_done = (job->state == THREADQUEUE_JOB_STATE_DONE);
    PTHREAD_UNLOCK(&job->lock);
    
    if (!job_done) {
      PTHREAD_COND_BROADCAST(&(threadqueue->cond));
      PTHREAD_COND_WAIT(&threadqueue->cb_cond, &threadqueue->lock);
    }
  } while (!job_done);

  // Free jobs submitted before this job.
  int i = 0;
  for (; i < threadqueue->queue_count; ++i) {
    if (threadqueue->queue[i] == job) break;
    kvz_threadqueue_free_job(&threadqueue->queue[i]);
  }
  // Move remaining jobs to the beginning of the array.
  if (i > 0) {
    threadqueue->queue_count -= i;
    threadqueue->queue_start = 0;
    memmove(threadqueue->queue, &threadqueue->queue[i], threadqueue->queue_count * sizeof(*threadqueue->queue));
    FILL_ARRAY(&threadqueue->queue[threadqueue->queue_count], 0, i);
  }

  PTHREAD_UNLOCK(&threadqueue->lock);
  
  return 1;
}

threadqueue_job_t * kvz_threadqueue_submit(threadqueue_queue_t * const threadqueue, void (*fptr)(void *arg), void *arg, int wait, const char* const debug_description) {
  threadqueue_job_t *job;
  //No lock here... this should be constant
  if (threadqueue->threads_count == 0) {
    //FIXME: This should be improved in order to handle dependencies
    PERFORMANCE_MEASURE_START(KVZ_PERF_JOB);
    fptr(arg);
    PERFORMANCE_MEASURE_END(KVZ_PERF_JOB, threadqueue, "%s", debug_description);
    return NULL;
  }
  
  assert(wait == 0 || wait == 1);
  
  job = MALLOC(threadqueue_job_t, 1);
  job->refcount = 1;
  
#ifdef KVZ_DEBUG
  if (debug_description) {
    size_t desc_len = MIN(255, strlen(debug_description));
    char* desc;
    
    //Copy description
    desc = MALLOC(char, desc_len + 1);
    assert(desc);
    memcpy(desc, debug_description, desc_len);
    desc[desc_len] = 0;
    
    job->debug_description = desc;
  } else {
    char* desc;
    desc = MALLOC(char, 255);
    sprintf(desc, "(*%p)(%p)", fptr, arg);
    
    job->debug_description = desc;
  }
  KVZ_GET_TIME(&job->debug_clock_enqueue);
#endif //KVZ_DEBUG
  
  if (!job) {
    fprintf(stderr, "Could not alloc job!\n");
    assert(0);
    return NULL;
  }
  
  job->fptr = fptr;
  job->arg = arg;
  if (pthread_mutex_init(&job->lock, NULL) != 0) {
    fprintf(stderr, "pthread_mutex_init(job) failed!\n");
    assert(0);
    return NULL;
  }
  job->ndepends = wait;
  job->rdepends = NULL;
  job->rdepends_count = 0;
  job->rdepends_size = 0;
  job->state = THREADQUEUE_JOB_STATE_QUEUED;
  
  PTHREAD_LOCK(&threadqueue->lock);
  
  //Add the reverse dependency
  if (threadqueue->queue_count >= threadqueue->queue_size) {
    threadqueue->queue = realloc(threadqueue->queue, sizeof(threadqueue_job_t *) * (threadqueue->queue_size + THREADQUEUE_LIST_REALLOC_SIZE));
    if (!threadqueue->queue) {
      fprintf(stderr, "Could not realloc queue!\n");
      assert(0);
      return NULL;
    }
    threadqueue->queue_size += THREADQUEUE_LIST_REALLOC_SIZE;
  }
  threadqueue->queue[threadqueue->queue_count++] = job;
  job->refcount++;
  
  if (job->ndepends == 0) {
    ++threadqueue->queue_waiting_execution;
    //Hope a thread can do it...
    PTHREAD_COND_SIGNAL(&(threadqueue->cond));
  } else {
    ++threadqueue->queue_waiting_dependency;
  }
  
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  return job;
}

int kvz_threadqueue_job_dep_add(threadqueue_job_t *job, threadqueue_job_t *depends_on) {
  //If we are not using threads, job are NULL pointers, so we can skip that
  if (!job && !depends_on) return 1;
  
  assert(job && depends_on);
  
  //Lock first the job, and then the dependency
  PTHREAD_LOCK(&job->lock);
  PTHREAD_LOCK(&depends_on->lock);
  
  if (depends_on->state != THREADQUEUE_JOB_STATE_DONE) {
    job->ndepends++;
  }
  
  //Add the reverse dependency (FIXME: this may be moved in the if above... but we would lose ability to track)
  if (depends_on->rdepends_count >= depends_on->rdepends_size) {
    depends_on->rdepends = realloc(depends_on->rdepends, sizeof(threadqueue_job_t *) * (depends_on->rdepends_size + THREADQUEUE_LIST_REALLOC_SIZE));
    if (!depends_on->rdepends) {
      fprintf(stderr, "Could not realloc rdepends!\n");
      assert(0);
      return 0;
    }
    depends_on->rdepends_size += THREADQUEUE_LIST_REALLOC_SIZE;
  }
  depends_on->rdepends[depends_on->rdepends_count++] = job;
  
  PTHREAD_UNLOCK(&depends_on->lock);
  PTHREAD_UNLOCK(&job->lock);
  
  return 1;
}

int kvz_threadqueue_job_unwait_job(threadqueue_queue_t * const threadqueue, threadqueue_job_t *job) {
  int ndepends = 0;
  
  //NULL job =>  no threads, nothing to do
  if (!job) return 1;
  PTHREAD_LOCK(&job->lock);
  job->ndepends--;
  ndepends = job->ndepends;
  PTHREAD_UNLOCK(&job->lock);
  
  if (ndepends == 0) {
    PTHREAD_LOCK(&threadqueue->lock);
    assert(threadqueue->queue_waiting_dependency > 0);
    --threadqueue->queue_waiting_dependency;
    ++threadqueue->queue_waiting_execution;
    //Hope a thread can do it...
    PTHREAD_COND_SIGNAL(&(threadqueue->cond));
    
    PTHREAD_UNLOCK(&threadqueue->lock);
  }
  
  return 1;
}

#ifdef KVZ_DEBUG
int threadqueue_log(threadqueue_queue_t * threadqueue, const KVZ_CLOCK_T *start, const KVZ_CLOCK_T *stop, const char* debug_description) {
  int i, thread_id = -1;
  FILE* output;
  
  assert(start);
  
  if (threadqueue) {
    //We need to lock to output safely
    PTHREAD_LOCK(&threadqueue->lock);
    
    output = threadqueue->debug_log;
    
    //Find the thread
    for(i = 0; i < threadqueue->threads_count; i++) {
      if(pthread_equal(threadqueue->threads[i], pthread_self()) != 0) {
        thread_id = i;
        break;
      }
    }
  } else {
    thread_id = -1;
    output = stderr;
  }
  
  if (thread_id >= 0) {
    if (stop) {
      fprintf(output, "\t%d\t-\t%lf\t+%lf\t-\t%s\n", thread_id, KVZ_CLOCK_T_AS_DOUBLE(*start), KVZ_CLOCK_T_DIFF(*start, *stop), debug_description);
    } else {
      fprintf(output, "\t%d\t-\t%lf\t-\t-\t%s\n", thread_id, KVZ_CLOCK_T_AS_DOUBLE(*start), debug_description);
    }
  } else {
    if (stop) {
      fprintf(output, "\t\t-\t%lf\t+%lf\t-\t%s\n", KVZ_CLOCK_T_AS_DOUBLE(*start), KVZ_CLOCK_T_DIFF(*start, *stop), debug_description);
    } else {
      fprintf(output, "\t\t-\t%lf\t-\t-\t%s\n", KVZ_CLOCK_T_AS_DOUBLE(*start), debug_description);
    }
  }
  
  if (threadqueue) {
    PTHREAD_UNLOCK(&threadqueue->lock);
  }
  return 1;
}
#endif //KVZ_DEBUG
