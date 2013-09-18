/**
 * \file
 * \brief Handles parsing and storing of configuration of the encoder.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "global.h"


/*!
    \brief Struct which contains all configuration data
*/
typedef struct
{
  char *input;      /*!< \brief Pointer to input filename  */
  char *output;     /*!< \brief Pointer to output filename */
  char *debug;      /*!< \brief Pointer to debug output    */
  int32_t frames;  /*!< \brief Number of frames to decode */
  int32_t width;   /*!< \brief frame width */
  int32_t height;  /*!< \brief frame height */
} config;

/* Function definitions */
config* config_alloc();
int config_init(config* cfg);
int config_destroy(config* cfg);
int config_read(config* cfg,int argc, char* argv[]);

#endif