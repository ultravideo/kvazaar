#ifndef INTER_H_
#define INTER_H_
/**
 * \file
 * \brief Handling Coding Units (CU's) for inter frames.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "picture.h"
#include "encoder.h"


void inter_set_block(picture* pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, cu_info *cur_cu);
void inter_recon(picture *ref,int32_t xpos, int32_t ypos,int32_t width, const int16_t mv[2], picture* dst);

void inter_get_mv_cand(encoder_control *encoder, int32_t x_cu, int32_t y_cu, int8_t depth, int16_t mv_cand[2][2]);

#endif
