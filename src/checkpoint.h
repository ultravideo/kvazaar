#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_
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
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \file
 * Printing of debug information.
 */

#ifdef CHECKPOINTS
#ifdef NDEBUG
#error "CHECKPOINTS require assertions to be enabled!"
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "global.h" // IWYU pragma: keep


extern FILE* g_ckpt_file;
extern int g_ckpt_enabled; //Do we check?
extern int g_ckpt_record; //Do we record?

#define CHECKPOINTS_INIT() do { \
  if (getenv("CHECKPOINTS")) {\
    if (strcmp(getenv("CHECKPOINTS"),"record") == 0) { \
      g_ckpt_file = fopen("__debug_ckpt.log", "w"); assert(g_ckpt_file); \
      g_ckpt_record = 1; \
    } else if (strcmp(getenv("CHECKPOINTS"),"check") == 0) { \
      g_ckpt_file = fopen("__debug_ckpt.log", "r"); assert(g_ckpt_file); \
      g_ckpt_record = 0; \
      g_ckpt_enabled = 0; \
    } else { \
      g_ckpt_file = NULL; \
    } \
  } else {\
    g_ckpt_file = NULL; \
  } \
} while (0)
    
#define CHECKPOINTS_FINALIZE() do {if (g_ckpt_file) fclose(g_ckpt_file); g_ckpt_file = NULL;} while (0)

#define CHECKPOINT_MARK(str, ...) do { \
  if (g_ckpt_file) { \
    if (g_ckpt_record) { \
      fprintf(g_ckpt_file, "MARK: " str "\n", __VA_ARGS__); \
    } else { \
      char buffer_ckpt[4096]; \
      long pos; \
      snprintf(buffer_ckpt, 4095, "MARK: " str "\n", __VA_ARGS__); \
      pos = ftell(g_ckpt_file); \
      g_ckpt_enabled = 0; \
      while (!feof(g_ckpt_file)) { \
        char buffer_file[4096]; \
        assert(fgets(buffer_file, 4095, g_ckpt_file) != NULL); \
        if (strncmp(buffer_file, buffer_ckpt, 4096)==0) { \
          g_ckpt_enabled = 1; \
          break; \
        } \
      } \
      if (!g_ckpt_enabled) fseek(g_ckpt_file, pos, SEEK_SET); \
    } \
  } \
} while (0)
#define CHECKPOINT(str, ...) do { \
  if (g_ckpt_file) { \
    if (g_ckpt_record) { \
      fprintf(g_ckpt_file, str "\n", __VA_ARGS__); \
    } else if (g_ckpt_enabled) { \
      char buffer_file[4096], buffer_ckpt[4096]; \
      assert(fgets(buffer_file, 4095, g_ckpt_file) != NULL); \
      snprintf(buffer_ckpt, 4095, str "\n", __VA_ARGS__); \
      if (strncmp(buffer_file, buffer_ckpt, 4096)!=0) { \
        fprintf(stderr, "Checkpoint failed (at %ld):\nFile: %sExec: %s", ftell(g_ckpt_file), buffer_file, buffer_ckpt); \
        assert(0); \
      } \
    } \
  } \
} while (0)
#endif

#if !defined(CHECKPOINTS)
#define CHECKPOINTS_INIT() 
#define CHECKPOINTS_FINALIZE() 
#define CHECKPOINT_MARK(str, ...) 
#define CHECKPOINT(str, ...) 
#endif



#endif //CHECKPOINT_H_
