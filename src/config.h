/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file config.h
    \brief Configuration header
    \author Marko Viitanen
    \date 2012-05
    
    Contains all configuration system related functions and structs
*/

#ifndef _CONFIG_H_
#define _CONFIG_H_

/*!
    \brief Struct which contains all configuration data
*/
typedef struct
{
  char *input;      /*!< \brief Pointer to input filename  */
  char *output;     /*!< \brief Pointer to output filename */
  char *debug;      /*!< \brief Pointer to debug output    */
  uint32_t frames;  /*!< \brief Number of frames to decode */
  uint32_t width;   /*!< \brief frame width */
  uint32_t height;  /*!< \brief frame height */
} config;

/* Function definitions */
config* config_alloc();
int config_init(config* cfg);
int config_destroy(config* cfg);
int config_read(config* cfg,int argc, char* argv[]);

#endif