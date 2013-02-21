/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file config.c
    \brief Configuration related functions
    \author Marko Viitanen
    \date 2012-05
    
    This file has all configuration related functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"


/*!
    \brief Allocate memory for config object
    \return pointer to allocated memory
*/
config* config_alloc()
{
  config* cfg = (config *)malloc(sizeof(config));
  memset(cfg, 0, sizeof(config));
  return cfg;
}

/*!
    \brief Initialize config structure
    \param cfg config object
    \return 1 on success, 0 on failure
*/
int config_init(config* cfg)
{  
  cfg->input  = NULL;
  cfg->output = NULL;
  cfg->debug  = NULL;
  cfg->frames = 0;
  cfg->width  = 320;
  cfg->height = 240;

  return 1;
}

/*!
    \brief Free memory allocated to the config
    \param cfg config object
    \return 1 on success, 0 on failure
*/
int config_destroy(config* cfg)
{
  free_pointer(cfg->input);
  free_pointer(cfg->output);
  free(cfg);

  return 1;
}

/*!
    \brief Read configuration options from argv to the config struct
    \param cfg config object
    \param argc argument count
    \param argv argument list
    \return 1 on success, 0 on failure
*/
int config_read(config* cfg,int argc, char* argv[])
{
  uint32_t pos = 0;
  int arg = 1;
  char option = 0;
  
  while(arg < argc && pos < strlen(argv[arg]))
  {
    /* Check for an option */
    if(argv[arg][0] == '-' && strlen(argv[arg]) > 1)
    {
      /* Second letter of the argument is the option we want to use */
      option = argv[arg][1];
      
      /* Point to the next argument */
      arg++;
      switch(option)
      {
        case 'i': /* Input */
          /* Allocate +1 for \0 */
          cfg->input = (char *)malloc(strlen(argv[arg])+1);
          memcpy(cfg->input, argv[arg], strlen(argv[arg])+1);
          break;
        case 'o': /* Output */
          cfg->output = (char *)malloc(strlen(argv[arg])+1);
          memcpy(cfg->output, argv[arg], strlen(argv[arg])+1);
          break;
        case 'd': /* Debug */
          cfg->debug = (char *)malloc(strlen(argv[arg])+1);
          memcpy(cfg->debug, argv[arg], strlen(argv[arg])+1);
          break;
        case 'w': /* width */
          /* Get picture width */
          cfg->width = atoi(argv[arg]);
          break;
        case 'h': /* height */
          /* Get picture height */
          cfg->height = atoi(argv[arg]);
          break;
        case 'n': /* Framecount */
          /* Get frame count to encode */
          cfg->frames = atoi(argv[arg]);
          break;
        default:
          /* Unknown command, print error message and ignore */
          fprintf(stderr, "%c is not a known option\r\n", option);
          break;
      }
    }
    /* Next argument */
    arg++;
  }
  
  /* Check that the required files were defined */
  if(cfg->input == NULL || cfg->output == NULL)
  {
    return 0;
  }
  
  return 1;
}