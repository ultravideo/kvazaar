#ifndef CONFIG_H_
#define CONFIG_H_
/**
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
 * \brief Handles parsing and storing of configuration of the encoder.
 */

#include "global.h"


/*!
    \brief Struct which contains all configuration data
*/
typedef struct
{
  char *input;      /*!< \brief Pointer to input filename  */
  char *output;     /*!< \brief Pointer to output filename */
  char *debug;      /*!< \brief Pointer to debug output    */
  int8_t qp;        /*!< \brief Quantization parameter */
  int16_t intra_period; /*!< \brief the period of intra frames in stream */
  int32_t frames;  /*!< \brief Number of frames to decode */
  int32_t width;   /*!< \brief frame width */
  int32_t height;  /*!< \brief frame height */
  int8_t deblock_enable; /*!< \brief Flag to enable deblocking filter */
  int8_t deblock_beta;   /*!< \brief (deblocking) beta offset (div 2), range -6...6 */
  int8_t deblock_tc;     /*!< \brief (deblocking) tc offset (div 2), range -6...6 */
} config;

/* Function definitions */
config *config_alloc();
int config_init(config *cfg);
int config_destroy(config *cfg);
int config_read(config *cfg,int argc, char *argv[]);

#endif
