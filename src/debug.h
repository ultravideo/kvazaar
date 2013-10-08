#ifndef DEBUG_H_
#define DEBUG_H_
/**
 * \file
 * \brief Tools for visualizing and debugging the encoder.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include <stdio.h>

#include "encoder.h"


FILE * open_cu_file(char *filename);
void close_cu_file(FILE *fp);
unsigned render_cu_file(encoder_control *encoder, picture *pic, unsigned depth, uint16_t x_cu, uint16_t y_cu, FILE *fp);

#endif
