#ifndef THREADQUEUE_H_
#define THREADQUEUE_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include <pthread.h>

typedef enum {
  THREADQUEUE_JOB_STATE_QUEUED = 0,
  THREADQUEUE_JOB_STATE_RUNNING = 1,
  THREADQUEUE_JOB_STATE_DONE = 2
} threadqueue_job_state;

typedef struct threadqueue_job {
  pthread_mutex_t lock;
  
  threadqueue_job_state state;
  
  unsigned int ndepends; //Number of job on which this job depends
  
  struct threadqueue_job **rdepends; //array of pointer to jobs that depend on this one. They have to exist when the thread finishes, because they cannot be run before.
  unsigned int rdepends_count; //number of rdepends
  unsigned int rdepends_size; //allocated size of rdepends
  
  //Job function and state to use
  void (*fptr)(void *arg);
  void *arg;
} threadqueue_job;


  

typedef struct {
  pthread_mutex_t lock;
  pthread_cond_t cond;
  pthread_cond_t cb_cond;
  
  pthread_t *threads;
  int threads_count;
  int threads_running;

  int stop; //=>1: threads should stop asap
  
  threadqueue_job **queue;
  unsigned int queue_count;
  unsigned int queue_size;
  unsigned int queue_waiting;
} threadqueue_queue;

//Init a threadqueue
int threadqueue_init(threadqueue_queue * threadqueue, int thread_count);

//Add a job to the queue, and returs a threadqueue_job handle. If wait == 1, one has to run threadqueue_job_unwait_job in order to have it run
threadqueue_job * threadqueue_submit(threadqueue_queue * threadqueue, void (*fptr)(void *arg), void *arg, int wait);

int threadqueue_job_unwait_job(threadqueue_queue * threadqueue, threadqueue_job *job);

//Add a dependency between two jobs.
int threadqueue_job_dep_add(threadqueue_job *job, threadqueue_job *depends_on);

//Blocking call until the queue is empty. Previously set threadqueue_job handles should not be used anymore
int threadqueue_flush(threadqueue_queue * threadqueue);

//Free ressources in a threadqueue
int threadqueue_finalize(threadqueue_queue * threadqueue);

/* Constraints: 
 * 
 * - Always first lock threadqueue, than a job inside it
 * - When job A depends on job B, always lock first job B and then job A
 * 
 * */

#endif //THREADQUEUE_H_
