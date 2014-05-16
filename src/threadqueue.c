
#include <assert.h>
#include <pthread.h>
#include <stdlib.h>

#ifdef _DEBUG
#include <string.h>
#endif //_DEBUG

#include "global.h"
#include "threadqueue.h"
#include "threads.h"

typedef struct {
  threadqueue_queue * threadqueue;
  int worker_id;
} threadqueue_worker_spec;

#define THREADQUEUE_LIST_REALLOC_SIZE 32

//#define PTHREAD_LOCK(l) fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p) (try)\n", __FUNCTION__, __LINE__, #l, l); if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s=%p) failed!\n", #l, l); assert(0); return 0; } else {fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}
//#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s=%p) failed!\n", #l, l); assert(0); return 0; }  else {fprintf(stderr, "%s:%d pthread_mutex_unlock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}

#define PTHREAD_LOCK(l) if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s) failed!\n", #l); assert(0); return 0; }
#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s) failed!\n", #l); assert(0); return 0; }

static void* threadqueue_worker(void* threadqueue_worker_spec_opaque) {
  threadqueue_worker_spec * const threadqueue_worker_spec = threadqueue_worker_spec_opaque;
  threadqueue_queue * const threadqueue = threadqueue_worker_spec->threadqueue;
  
#ifdef _DEBUG
  GET_TIME(&threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id]);
#endif //_DEBUG

  for(;;) {
    int task_id = -1, i = 0;
    
    PTHREAD_LOCK(&threadqueue->lock);

    while(!threadqueue->stop && threadqueue->queue_waiting == 0) {
      if (pthread_cond_wait(&threadqueue->cond, &threadqueue->lock) != 0) {
        fprintf(stderr, "pthread_cond_wait failed!\n");
        assert(0);
        return 0;
      }
    }

    if(threadqueue->stop) {
        break;
    }
    
    //Find a task (should be fast enough)
    task_id = -1;
    for (i = 0; i < threadqueue->queue_count; ++i) {
      threadqueue_job * const job = threadqueue->queue[i];
      PTHREAD_LOCK(&job->lock);
      
      if (job->state == THREADQUEUE_JOB_STATE_QUEUED && job->ndepends == 0) {
        task_id = i;
        break; //Task remains locked
      }
      
      //Not this task, so unlock it
      PTHREAD_UNLOCK(&job->lock);
    }
    
    //Ok we got a job (and we have a lock on it)
    if (task_id != -1) {
      threadqueue_job * const job = threadqueue->queue[i];

      assert(job->state == THREADQUEUE_JOB_STATE_QUEUED);
      job->state = THREADQUEUE_JOB_STATE_RUNNING;
      
      --threadqueue->queue_waiting;
      
      //We can unlock the job here, since fptr and arg are constant
      PTHREAD_UNLOCK(&job->lock);
      //Unlock the queue
      PTHREAD_UNLOCK(&threadqueue->lock);
      
#ifdef _DEBUG
      job->debug_worker_id = threadqueue_worker_spec->worker_id;
      GET_TIME(&job->debug_clock_start);
#endif //_DEBUG
      
      job->fptr(job->arg);
      
#ifdef _DEBUG
      job->debug_worker_id = threadqueue_worker_spec->worker_id;
      GET_TIME(&job->debug_clock_stop);
#endif //_DEBUG
      
      //Re-lock the job to update its status and treat its dependencies
      PTHREAD_LOCK(&job->lock);
      assert(job->state == THREADQUEUE_JOB_STATE_RUNNING);
      
      job->state = THREADQUEUE_JOB_STATE_DONE;
      
      //Decrease counter of dependencies
      for (i = 0; i < job->rdepends_count; ++i) {
        threadqueue_job * const depjob = job->rdepends[i];
        //Note that we lock the dependency AFTER locking the source. This avoids a deadlock in dep_add
        PTHREAD_LOCK(&depjob->lock);
        
        assert(depjob->state == THREADQUEUE_JOB_STATE_QUEUED);
        assert(depjob->ndepends > 0);
        --depjob->ndepends;
        
        PTHREAD_UNLOCK(&depjob->lock);
      }
      //Unlock the job
      PTHREAD_UNLOCK(&job->lock);
      
      //Signal the queue that we've done a job
      PTHREAD_LOCK(&threadqueue->lock);
      pthread_cond_broadcast(&threadqueue->cb_cond);
      PTHREAD_UNLOCK(&threadqueue->lock);
    } else {
      PTHREAD_UNLOCK(&threadqueue->lock);
    }
  }

  //We got out of the loop because threadqueue->stop == 1. The queue is locked.
  assert(threadqueue->stop);
  --threadqueue->threads_running;
  
#ifdef _DEBUG
  GET_TIME(&threadqueue->debug_clock_thread_end[threadqueue_worker_spec->worker_id]);
  
  fprintf(threadqueue->debug_log, "\t%d\t-\t%lf\t+%lf\t-\tthread\n", threadqueue_worker_spec->worker_id, CLOCK_T_AS_DOUBLE(threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id]), CLOCK_T_DIFF(threadqueue->debug_clock_thread_start[threadqueue_worker_spec->worker_id], threadqueue->debug_clock_thread_end[threadqueue_worker_spec->worker_id]));
#endif //_DEBUG
  
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  pthread_exit(NULL);
  
  return NULL;
}
int threadqueue_init(threadqueue_queue * const threadqueue, int thread_count) {
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
  threadqueue->threads_running = 0;
  threadqueue->threads_count = thread_count;
  
  threadqueue->threads = MALLOC(pthread_t, thread_count);
  if (!threadqueue->threads) {
    fprintf(stderr, "Could not malloc threadqueue->threads!\n");
    return 0;
  }
#ifdef _DEBUG
  threadqueue->debug_clock_thread_start = MALLOC(CLOCK_T, thread_count);
  assert(threadqueue->debug_clock_thread_start);
  threadqueue->debug_clock_thread_end = MALLOC(CLOCK_T, thread_count);
  assert(threadqueue->debug_clock_thread_end);
  threadqueue->debug_log = fopen("threadqueue.log", "w");
#endif //_DEBUG
  
  threadqueue->queue = NULL;
  threadqueue->queue_size = 0;
  threadqueue->queue_count = 0;
  threadqueue->queue_waiting = 0;
  
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
      return 0;
    }
  }
  
  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
}

static void threadqueue_free_jobs(threadqueue_queue * const threadqueue) {
  int i;
  for (i=0; i < threadqueue->queue_count; ++i) {
#ifdef _DEBUG
    int j;
    GET_TIME(&threadqueue->queue[i]->debug_clock_dequeue);
    fprintf(threadqueue->debug_log, "%p\t%d\t%lf\t+%lf\t+%lf\t+%lf\t%s\n", threadqueue->queue[i], threadqueue->queue[i]->debug_worker_id, CLOCK_T_AS_DOUBLE(threadqueue->queue[i]->debug_clock_enqueue), CLOCK_T_DIFF(threadqueue->queue[i]->debug_clock_enqueue, threadqueue->queue[i]->debug_clock_start), CLOCK_T_DIFF(threadqueue->queue[i]->debug_clock_start, threadqueue->queue[i]->debug_clock_stop), CLOCK_T_DIFF(threadqueue->queue[i]->debug_clock_stop, threadqueue->queue[i]->debug_clock_dequeue), threadqueue->queue[i]->debug_description);
    
    for (j = 0; j < threadqueue->queue[i]->rdepends_count; ++j) {
      fprintf(threadqueue->debug_log, "%p->%p\n", threadqueue->queue[i], threadqueue->queue[i]->rdepends[j]);
    }
    
    FREE_POINTER(threadqueue->queue[i]->debug_description);
#endif
    FREE_POINTER(threadqueue->queue[i]);
  }
  threadqueue->queue_count = 0;
#ifdef _DEBUG
  {
    CLOCK_T time;
    GET_TIME(&time);
   
    fprintf(threadqueue->debug_log, "\t\t-\t-\t%lf\t-\tFLUSH\n", CLOCK_T_AS_DOUBLE(time));
  }
#endif
}

int threadqueue_finalize(threadqueue_queue * const threadqueue) {
  int i;
  
  //Flush the queue
  if (!threadqueue_flush(threadqueue)) {
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
  for(i = 0; i < threadqueue->threads_count; i++) {
    if(pthread_join(threadqueue->threads[i], NULL) != 0) {
      fprintf(stderr, "pthread_join failed!\n");
      return 0;
    }
  }
  
#ifdef _DEBUG
  fclose(threadqueue->debug_log);
#endif
  
  //Free allocated stuff
  FREE_POINTER(threadqueue->queue);
  threadqueue->queue_count = 0;
  threadqueue->queue_size = 0;
  
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

int threadqueue_flush(threadqueue_queue * const threadqueue) {
  int notdone = 1;
  int i;
  
  //Lock the queue
  PTHREAD_LOCK(&threadqueue->lock);
  
  do {
    if (threadqueue->queue_waiting > 0) {
      notdone = threadqueue->queue_waiting;
    } else {
      notdone = 0;
      for (i = 0; i < threadqueue->queue_count; ++i) {
        PTHREAD_LOCK(&threadqueue->queue[i]->lock);
        if (threadqueue->queue[i]->state != THREADQUEUE_JOB_STATE_DONE) {
          notdone++;
        }
        PTHREAD_UNLOCK(&threadqueue->queue[i]->lock);
      }
    }

    if (notdone > 0) {
      //Give threads a change to unlock if needed
      if (pthread_cond_broadcast(&(threadqueue->cond)) != 0) {
        fprintf(stderr, "pthread_cond_broadcast failed!\n");
        assert(0);
        return 0;
      }
      SLEEP();
      if (pthread_cond_wait(&threadqueue->cb_cond, &threadqueue->lock) != 0) {
        fprintf(stderr, "pthread_cond_wait failed!\n");
        assert(0); //FIXME
        return 0;
      }
    }
  } while (notdone > 0);
  
  threadqueue_free_jobs(threadqueue);

  assert(threadqueue->queue_waiting == 0);

  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
}

threadqueue_job * threadqueue_submit(threadqueue_queue * const threadqueue, void (*fptr)(void *arg), void *arg, int wait, const char* const debug_description) {
  threadqueue_job *job;
  //No lock here... this should be constant
  if (threadqueue->threads_count == 0) {
    //FIXME: This should be improved in order to handle dependencies
    PERFORMANCE_MEASURE_START();
    fptr(arg);
    PERFORMANCE_MEASURE_END(threadqueue, "%s", debug_description);
    return NULL;
  }
  
  assert(wait == 0 || wait == 1);
  
  job = MALLOC(threadqueue_job, 1);
  
#ifdef _DEBUG
  if (debug_description) {
    int desc_len = MIN(255, strlen(debug_description));
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
  GET_TIME(&job->debug_clock_enqueue);
#endif //_DEBUG
  
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
    threadqueue->queue = realloc(threadqueue->queue, sizeof(threadqueue_job *) * (threadqueue->queue_size + THREADQUEUE_LIST_REALLOC_SIZE));
    if (!threadqueue->queue) {
      fprintf(stderr, "Could not realloc queue!\n");
      assert(0);
      return NULL;
    }
    threadqueue->queue_size += THREADQUEUE_LIST_REALLOC_SIZE;
  }
  threadqueue->queue[threadqueue->queue_count++] = job;
  
  ++threadqueue->queue_waiting;
  
  //Hope a thread can do it...
  if(pthread_cond_signal(&(threadqueue->cond)) != 0) {
      fprintf(stderr, "pthread_cond_signal failed!\n");
      assert(0);
      return NULL;
  }
  
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  return job;
}

int threadqueue_job_dep_add(threadqueue_job *job, threadqueue_job *depends_on) {
  //If we are not using threads, job are NULL pointers, so we can skip that
  if (!job && !depends_on) return 1;
  //Lock first the job, and then the dependency
  PTHREAD_LOCK(&job->lock);
  PTHREAD_LOCK(&depends_on->lock);
  
  if (depends_on->state != THREADQUEUE_JOB_STATE_DONE) {
    job->ndepends++;
  }
  
  //Add the reverse dependency (FIXME: this may be moved in the if above... but we would lose ability to track)
  if (depends_on->rdepends_count >= depends_on->rdepends_size) {
    depends_on->rdepends = realloc(depends_on->rdepends, sizeof(threadqueue_job *) * (depends_on->rdepends_size + THREADQUEUE_LIST_REALLOC_SIZE));
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

int threadqueue_job_unwait_job(threadqueue_queue * const threadqueue, threadqueue_job *job) {
  //NULL job =>  no threads, nothing to do
  if (!job) return 1;
  PTHREAD_LOCK(&job->lock);
  job->ndepends--;
  PTHREAD_UNLOCK(&job->lock);
  
  PTHREAD_LOCK(&threadqueue->lock);
  //Hope a thread can do it...
  if(pthread_cond_signal(&(threadqueue->cond)) != 0) {
      fprintf(stderr, "pthread_cond_signal failed!\n");
      assert(0);
      return 0;
  }
  PTHREAD_UNLOCK(&threadqueue->lock);
  
  return 1;
}

#ifdef _DEBUG
int threadqueue_log(threadqueue_queue * threadqueue, const CLOCK_T *start, const CLOCK_T *stop, const char* debug_description) {
  int i, thread_id = -1;
  
  assert(start);
  
  //We need to lock to output safely
  PTHREAD_LOCK(&threadqueue->lock);
  
  //Find the thread
  for(i = 0; i < threadqueue->threads_count; i++) {
    if(pthread_equal(threadqueue->threads[i], pthread_self()) != 0) {
      thread_id = i;
      break;
    }
  }
  
  if (thread_id >= 0) {
    if (stop) {
      fprintf(threadqueue->debug_log, "\t%d\t-\t%lf\t+%lf\t-\t%s\n", thread_id, CLOCK_T_AS_DOUBLE(*start), CLOCK_T_DIFF(*start, *stop), debug_description);
    } else {
      fprintf(threadqueue->debug_log, "\t%d\t-\t%lf\t-\t-\t%s\n", thread_id, CLOCK_T_AS_DOUBLE(*start), debug_description);
    }
  } else {
    if (stop) {
      fprintf(threadqueue->debug_log, "\t\t-\t%lf\t+%lf\t-\t%s\n", CLOCK_T_AS_DOUBLE(*start), CLOCK_T_DIFF(*start, *stop), debug_description);
    } else {
      fprintf(threadqueue->debug_log, "\t\t-\t%lf\t-\t-\t%s\n", CLOCK_T_AS_DOUBLE(*start), debug_description);
    }
  }
  PTHREAD_UNLOCK(&threadqueue->lock);
  return 1;
}
#endif //_DEBUG