/**
 * \file
 * \brief Searching of parameters for intra and inter frames.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#ifndef __SEARCH_H
#define __SEARCH_H

#include "global.h"

#include "encoder.h"


void search_slice_data(encoder_control* encoder);
void search_tree(encoder_control* encoder,uint16_t x_cu,uint16_t y_cu, uint8_t depth);
uint32_t search_best_mode(encoder_control* encoder,uint16_t x_cu,uint16_t y_cu, uint8_t depth);

#endif
