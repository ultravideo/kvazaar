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

/*
 * \file
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * \brief Allocate memory for config object
 * \return pointer to allocated memory
 */
config *config_alloc()
{
  config *cfg = (config *)malloc(sizeof(config));
  if (!cfg) {
    fprintf(stderr, "Failed to allocate a config object!\n");
    return cfg;
  }

  memset(cfg, 0, sizeof(config));
  return cfg;
}

/**
 * \brief Initialize config structure
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_init(config *cfg)
{  
  cfg->input  = NULL;
  cfg->output = NULL;
  cfg->debug  = NULL;
  cfg->frames = 0;
  cfg->width  = 320;
  cfg->height = 240;
  cfg->qp     = 32;
  cfg->intra_period = 0;

  return 1;
}

/**
 * \brief Free memory allocated to the config
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_destroy(config *cfg)
{
  FREE_POINTER(cfg->input);
  FREE_POINTER(cfg->output);
  free(cfg);

  return 1;
}

/**
 * \brief Allocates memory space for a string, and copies it
 * \param char * string to copy
 * \return a pointer to the copied string on success, null on failure
 */
char *copy_string(char *string)
{
  // Allocate +1 for \0
  char *allocated_string = (char *)malloc(strlen(string) + 1);
  if (!allocated_string) {
    fprintf(stderr, "Failed to allocate a string!\n");
    return allocated_string;
  }

  // Copy the string to the new buffer
  memcpy(allocated_string, string, strlen(string) + 1);

  return allocated_string;
}

/**
 * \brief Read configuration options from argv to the config struct
 * \param cfg config object
 * \param argc argument count
 * \param argv argument list
 * \return 1 on success, 0 on failure
 */
int config_read(config *cfg,int argc, char *argv[])
{
  uint32_t pos = 0;
  int arg = 1;
  char option = 0;
  
  while(arg < argc && pos < strlen(argv[arg])) {
    // Check for an option
    if(argv[arg][0] == '-' && strlen(argv[arg]) > 1) {
      // Second letter of the argument is the option we want to use
      option = argv[arg][1];
      
      // Point to the next argument
      arg++;
      switch(option) {
        case 'i': // Input
          cfg->input = copy_string(argv[arg]);
          break;
        case 'o': // Output
          cfg->output = copy_string(argv[arg]);
          break;
        case 'd': // Debug
          cfg->debug = copy_string(argv[arg]);
          break;
        case 'w': // width
          cfg->width = atoi(argv[arg]);
          break;
        case 'h': // height
          cfg->height = atoi(argv[arg]);
          break;
        case 'n': // number of frames to encode
          cfg->frames = atoi(argv[arg]);
          break;
        case 'q': // QP
          cfg->qp = atoi(argv[arg]);
          break;
        case 'p': // Intra period
          cfg->intra_period = atoi(argv[arg]);
          break;
        default:
          // Unknown command, print error message and ignore
          fprintf(stderr, "%c is not a known option\r\n", option);
          break;
      }
    }
    // Next argument
    arg++;
  }
  
  // Check that the required files were defined
  if(cfg->input == NULL || cfg->output == NULL) return 0;  
  
  return 1;
}
