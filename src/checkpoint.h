#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

/*****************************************************************************
 * Copyright (C) 2007-2014 Laurent Fasnacht <l@libres.ch>                    *
 *                                                                           *
 * You can redistribute it and/or modify it under the terms of the GNU       *
 * General Public License version 2 as published by the Free Software        *
 * Foundation.                                                               *
 ****************************************************************************/

#ifdef CHECKPOINTS
#ifdef NDEBUG
#error "CHECKPOINTS require assertions to be enabled!"
#endif
#include <stdio.h>
#include <stdlib.h>

extern FILE* g_ckpt_file;
extern int g_ckpt_enabled; //Do we check?
extern int g_ckpt_record; //Do we record?
#endif



#if defined(CHECKPOINTS)
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
    
#define CHECKPOINTS_FINALIZE() do {fclose(g_ckpt_file); g_ckpt_file = NULL;} while (0)

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
        fgets(buffer_file, 4095, g_ckpt_file); \
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
      fgets(buffer_file, 4095, g_ckpt_file); \
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
#define CHECKPOINTS_INIT() do {} while (0)
#define CHECKPOINTS_FINALIZE() do {} while (0)
#define CHECKPOINT_MARK(str, ...) do {} while (0)
#define CHECKPOINT(str, ...) do {} while (0)
#endif



#endif //CHECKPOINT_H_