#ifndef THREADQUEUE_H_
#define THREADQUEUE_H_
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

/**
 * \ingroup Threading
 * \file
 * Container for worker tasks.
 */

#include <pthread.h>

#include "global.h" // IWYU pragma: keep

typedef enum {
  THREADQUEUE_JOB_STATE_QUEUED = 0,
  THREADQUEUE_JOB_STATE_RUNNING = 1,
  THREADQUEUE_JOB_STATE_DONE = 2
} threadqueue_job_state;

typedef struct threadqueue_job_t {
  pthread_mutex_t lock;
  
  threadqueue_job_state state;
  
  unsigned int ndepends; //Number of active dependencies that this job wait for
  
  struct threadqueue_job_t **rdepends; //array of pointer to jobs that depend on this one. They have to exist when the thread finishes, because they cannot be run before.
  unsigned int rdepends_count; //number of rdepends
  unsigned int rdepends_size; //allocated size of rdepends

  // Reference count
  int refcount;

  //Job function and state to use
  void (*fptr)(void *arg);
  void *arg;
} threadqueue_job_t;


  

typedef struct {
  pthread_mutex_t lock;
  pthread_cond_t cond;
  pthread_cond_t cb_cond;
  
  pthread_t *threads;
  int threads_count;
  int threads_running;

  bool stop; // if true, threads should stop asap

  int fifo;
  
  threadqueue_job_t **queue;
  unsigned int queue_start;
  unsigned int queue_count;
  unsigned int queue_size;
  unsigned int queue_waiting_execution; //Number of jobs without any dependency which could be run
  unsigned int queue_waiting_dependency; //Number of jobs waiting for a dependency to complete
  unsigned int queue_running; //Number of jobs running

} threadqueue_queue_t;

//Init a threadqueue (if fifo, then behave as a FIFO with dependencies, otherwise as a LIFO with dependencies)
int kvz_threadqueue_init(threadqueue_queue_t * threadqueue, int thread_count, int fifo);

//Add a job to the queue, and returs a threadqueue_job handle. If wait == 1, one has to run kvz_threadqueue_job_unwait_job in order to have it run
threadqueue_job_t * kvz_threadqueue_submit(threadqueue_queue_t * threadqueue, void (*fptr)(void *arg), void *arg, int wait);

void kvz_threadqueue_free_job(threadqueue_job_t **job_ptr);

threadqueue_job_t *kvz_threadqueue_copy_ref(threadqueue_job_t *job);

int kvz_threadqueue_job_unwait_job(threadqueue_queue_t * threadqueue, threadqueue_job_t *job);

//Add a dependency between two jobs.
int kvz_threadqueue_job_dep_add(threadqueue_job_t *job, threadqueue_job_t *depends_on);

/**
 * \brief Stop all threads after they finish the current jobs.
 *
 * Blocks until all threads have stopped.
 */
int kvz_threadqueue_stop(threadqueue_queue_t * const threadqueue);

//Blocking call until job is executed. Job handles submitted before job should not be used any more as they are removed from the queue.
int kvz_threadqueue_waitfor(threadqueue_queue_t * threadqueue, threadqueue_job_t * job);

//Free ressources in a threadqueue
int kvz_threadqueue_finalize(threadqueue_queue_t * threadqueue);

/* Constraints: 
 * 
 * - Always first lock threadqueue, than a job inside it
 * - When job A depends on job B, always lock first job B and then job A
 * - Jobs should be submitted in an order which is compatible with serial execution.
 * 
 * */

#endif //THREADQUEUE_H_
