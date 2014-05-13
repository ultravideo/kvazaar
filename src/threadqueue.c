
#include <assert.h>
#include <pthread.h>
#include <stdlib.h>

#include "global.h"
#include "threadqueue.h"
#include "threads.h"

#define THREADQUEUE_LIST_REALLOC_SIZE 32

//#define PTHREAD_LOCK(l) fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p) (try)\n", __FUNCTION__, __LINE__, #l, l); if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s=%p) failed!\n", #l, l); assert(0); return 0; } else {fprintf(stderr, "%s:%d pthread_mutex_lock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}
//#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s=%p) failed!\n", #l, l); assert(0); return 0; }  else {fprintf(stderr, "%s:%d pthread_mutex_unlock(%s=%p)\n", __FUNCTION__, __LINE__, #l, l);}

#define PTHREAD_LOCK(l) if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s) failed!\n", #l); assert(0); return 0; }
#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s) failed!\n", #l); assert(0); return 0; }

static void* threadqueue_worker(void* threadqueue_opaque) {
  threadqueue_queue * const threadqueue = threadqueue_opaque;
  
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
      
      job->fptr(job->arg);
      
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
  
  threadqueue->queue = NULL;
  threadqueue->queue_size = 0;
  threadqueue->queue_count = 0;
  threadqueue->queue_waiting = 0;
  
  //Lock the queue before creating threads, to ensure they all have correct information
  PTHREAD_LOCK(&threadqueue->lock);
  
  for(i = 0; i < thread_count; i++) {
      if(pthread_create(&(threadqueue->threads[i]), NULL, threadqueue_worker, (void*)threadqueue) != 0) {
          fprintf(stderr, "pthread_create failed!\n");
          assert(0);
          return 0;
      }
      threadqueue->threads_running++;
  }
  
  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
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
  for (i=0; i < threadqueue->queue_count; ++i) {
    FREE_POINTER(threadqueue->queue[i]);
  }
  threadqueue->queue_count = 0;

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
  
#if 1
  //technically not needed, but spares memory. On the other hand, it makes debugging harder.
  for (i=0; i < threadqueue->queue_count; ++i) {
    FREE_POINTER(threadqueue->queue[i]);
  }
  threadqueue->queue_count = 0;
#endif
  assert(threadqueue->queue_waiting == 0);

  PTHREAD_UNLOCK(&threadqueue->lock);

  return 1;
}

threadqueue_job * threadqueue_submit(threadqueue_queue * const threadqueue, void (*fptr)(void *arg), void *arg) {
  threadqueue_job *job;
  //No lock here... this should be constant
  if (threadqueue->threads_count == 0) {
    fptr(arg);
    return NULL;
  }
  
  job = MALLOC(threadqueue_job, 1);
  
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
  job->ndepends = 0;
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