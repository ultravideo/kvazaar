#ifndef THREADS_H_
#define THREADS_H_
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

#ifdef __GNUC__

#define ATOMIC_INC(ptr)                     __sync_add_and_fetch((volatile int32_t*)ptr, 1)
#define ATOMIC_DEC(ptr)                     __sync_add_and_fetch((volatile int32_t*)ptr, -1)

#else //__GNUC__

//TODO: we assume !GCC => Windows... this may be bad
#define ATOMIC_INC(ptr)                     InterlockedIncrement((volatile LONG*)ptr)
#define ATOMIC_DEC(ptr)                     InterlockedDecrement((volatile LONG*)ptr)

#endif //__GNUC__

#endif //THREADS_H_
