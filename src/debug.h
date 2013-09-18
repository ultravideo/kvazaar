/**
 * HEVC Encoder
 * - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing 2013.
 * - Ari Koivula (ari at koivu.la ), Tampere University of Technology, Department of Pervasive Computing 2013.
 */

#ifndef __DEBUG_H
#define __DEBUG_H

#include "global.h"

#include <stdio.h>
#include "encoder.h"


FILE * open_cu_file(char *filename);
void close_cu_file(FILE *fp);
unsigned render_cu_file(encoder_control *encoder, unsigned depth, uint16_t xCtb, uint16_t yCtb, FILE *fp);

#endif
