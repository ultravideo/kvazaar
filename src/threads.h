#ifndef THREADS_H_
#define THREADS_H_
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
 * Abstractions for operating system specific stuff.
 */

#include "global.h"

#include <pthread.h>

#if defined(__GNUC__) && !defined(__MINGW32__) 
#include <unistd.h>
#include <time.h>

#define KVZ_CLOCK_T struct timespec

#ifdef __MACH__
// Workaround Mac OS not having clock_gettime.
#include <mach/clock.h>
#include <mach/mach.h>
#define KVZ_GET_TIME(clock_t) { \
  clock_serv_t cclock; \
  mach_timespec_t mts; \
  host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock); \
  clock_get_time(cclock, &mts); \
  mach_port_deallocate(mach_task_self(), cclock); \
  (clock_t)->tv_sec = mts.tv_sec; \
  (clock_t)->tv_nsec = mts.tv_nsec; \
}
#else
#define KVZ_GET_TIME(clock_t) { clock_gettime(CLOCK_MONOTONIC, (clock_t)); }
#endif

#define KVZ_CLOCK_T_AS_DOUBLE(ts) ((double)((ts).tv_sec) + (double)((ts).tv_nsec) / (double)1000000000L)
#define KVZ_CLOCK_T_DIFF(start, stop) ((double)((stop).tv_sec - (start).tv_sec) + (double)((stop).tv_nsec - (start).tv_nsec) / (double)1000000000L)

#define KVZ_ATOMIC_INC(ptr)                     __sync_add_and_fetch((volatile int32_t*)ptr, 1)
#define KVZ_ATOMIC_DEC(ptr)                     __sync_add_and_fetch((volatile int32_t*)ptr, -1)
#define KVZ_SLEEP()                             usleep(0)

#else //__GNUC__
//TODO: we assume !GCC => Windows... this may be bad
#include <windows.h>

#define KVZ_CLOCK_T struct _FILETIME
#define KVZ_GET_TIME(clock_t) { GetSystemTimeAsFileTime(clock_t); }
// _FILETIME has 32bit low and high part of 64bit 100ns resolution timestamp (since 12:00 AM January 1, 1601)
#define KVZ_CLOCK_T_AS_DOUBLE(ts) ((double)(((uint64_t)(ts).dwHighDateTime)<<32 | (uint64_t)(ts).dwLowDateTime) / (double)10000000L)
#define KVZ_CLOCK_T_DIFF(start, stop) ((double)((((uint64_t)(stop).dwHighDateTime)<<32 | (uint64_t)(stop).dwLowDateTime) - \
                                  (((uint64_t)(start).dwHighDateTime)<<32 | (uint64_t)(start).dwLowDateTime)) / (double)10000000L)


#define KVZ_ATOMIC_INC(ptr)                     InterlockedIncrement((volatile LONG*)ptr)
#define KVZ_ATOMIC_DEC(ptr)                     InterlockedDecrement((volatile LONG*)ptr)
// Sleep(0) results in bad performance on Windows for some reason,
// As a work around sleep for 10ms.
#define KVZ_SLEEP()                             Sleep(10)

#endif //__GNUC__

#endif //THREADS_H_
