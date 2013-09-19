#ifndef TRANSFORM_H_
#define TRANSFORM_H_
/**
 * \file
 * \brief Transformations, such as quantization and DST.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "encoder.h"


extern int32_t* g_quant_coeff[4][6][6];
extern const int32_t g_quant_intra_default_8x8[64];

void quant(encoder_control* encoder, int16_t* p_src, int16_t* p_des, int32_t width,
           int32_t height, uint32_t *ac_sum, int8_t type, int8_t scan_idx );
void dequant(encoder_control* encoder, int16_t* q_coef, int16_t* coef, int32_t width, int32_t height,int8_t type);

void transform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);
void itransform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);

void scalinglist_init();
void scalinglist_process_enc( int32_t *coeff, int32_t *quant_coeff, int32_t quant_scales, 
                             uint32_t height,uint32_t width, uint32_t ratio, int32_t size_num, uint32_t dc, uint8_t flat);
void scalinglist_process();
void scalinglist_set(int32_t *coeff, uint32_t list_id, uint32_t size_id, uint32_t qp);
void scalinglist_destroy();

#endif
