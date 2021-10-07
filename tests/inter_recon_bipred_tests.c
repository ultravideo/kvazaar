/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "greatest/greatest.h"

#include "test_strategies.h"
#include "strategies/generic/picture-generic.h"
#include <string.h>
#include <stdlib.h>


static lcu_t expected_test_result;
static lcu_t result;

static lcu_t lcu1;

int temp1, temp2, temp3, temp4;

int16_t mv_param[2][2] = { { 3,3 },{ 3,3 } };
int width = 16;
int height = 16;
int xpos = 0;
int ypos = 0;


kvz_pixel temp_lcu_y[LCU_WIDTH*LCU_WIDTH];
kvz_pixel temp_lcu_u[LCU_WIDTH_C*LCU_WIDTH_C];
kvz_pixel temp_lcu_v[LCU_WIDTH_C*LCU_WIDTH_C];

int hi_prec_luma_rec0;
int hi_prec_luma_rec1;
int hi_prec_chroma_rec0;
int hi_prec_chroma_rec1;

hi_prec_buf_t* high_precision_rec0 = 0;
hi_prec_buf_t* high_precision_rec1 = 0;

int temp_x, temp_y;



static void setup()
{

	memset(lcu1.rec.y, 0, sizeof(kvz_pixel) * 64 * 64);
	memset(lcu1.rec.u, 0, sizeof(kvz_pixel) * 32 * 32);
	memset(lcu1.rec.v, 0, sizeof(kvz_pixel) * 32 * 32);


	memset(expected_test_result.rec.y, 0, sizeof(kvz_pixel) * 64 * 64);
	memset(expected_test_result.rec.u, 0, sizeof(kvz_pixel) * 32 * 32);
	memset(expected_test_result.rec.v, 0, sizeof(kvz_pixel) * 32 * 32);

	memcpy(expected_test_result.rec.y, lcu1.rec.y, sizeof(kvz_pixel) * 64 * 64);
	memcpy(expected_test_result.rec.u, lcu1.rec.u, sizeof(kvz_pixel) * 32 * 32);
	memcpy(expected_test_result.rec.v, lcu1.rec.v, sizeof(kvz_pixel) * 32 * 32);

	// Setup is not optimized working function from picture-generic.c.

	
	int shift = 15 - KVZ_BIT_DEPTH;
	int offset = 1 << (shift - 1);

 hi_prec_luma_rec0 = mv_param[0][0] & 3 || mv_param[0][1] & 3;
 hi_prec_luma_rec1 = mv_param[1][0] & 3 || mv_param[1][1] & 3;

 hi_prec_chroma_rec0 = mv_param[0][0] & 7 || mv_param[0][1] & 7;
 hi_prec_chroma_rec1 = mv_param[1][0] & 7 || mv_param[1][1] & 7;

	if (hi_prec_chroma_rec0) high_precision_rec0 = kvz_hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);
	if (hi_prec_chroma_rec1) high_precision_rec1 = kvz_hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);

	


	for (temp_y = 0; temp_y < height; ++temp_y) {
		int y_in_lcu = ((ypos + temp_y) & ((LCU_WIDTH)-1));
		for (temp_x = 0; temp_x < width; ++temp_x) {
			int x_in_lcu = ((xpos + temp_x) & ((LCU_WIDTH)-1));
			int16_t sample0_y = (hi_prec_luma_rec0 ? high_precision_rec0->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (temp_lcu_y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			int16_t sample1_y = (hi_prec_luma_rec1 ? high_precision_rec1->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (expected_test_result.rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			expected_test_result.rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_y + sample1_y + offset) >> shift);
		}

	}
	for (temp_y = 0; temp_y < height >> 1; ++temp_y) {
		int y_in_lcu = (((ypos >> 1) + temp_y) & (LCU_WIDTH_C - 1));
		for (temp_x = 0; temp_x < width >> 1; ++temp_x) {
			int x_in_lcu = (((xpos >> 1) + temp_x) & (LCU_WIDTH_C - 1));
			int16_t sample0_u = (hi_prec_chroma_rec0 ? high_precision_rec0->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			int16_t sample1_u = (hi_prec_chroma_rec1 ? high_precision_rec1->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (expected_test_result.rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			expected_test_result.rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_u + sample1_u + offset) >> shift);

			int16_t sample0_v = (hi_prec_chroma_rec0 ? high_precision_rec0->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			int16_t sample1_v = (hi_prec_chroma_rec1 ? high_precision_rec1->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (expected_test_result.rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
			expected_test_result.rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_v + sample1_v + offset) >> shift);



		}
	}
}


TEST test_inter_recon_bipred()
{


	memcpy(result.rec.y, lcu1.rec.y, sizeof(kvz_pixel) * 64 * 64);
	memcpy(result.rec.u, lcu1.rec.u, sizeof(kvz_pixel) * 32 * 32);
	memcpy(result.rec.v, lcu1.rec.v, sizeof(kvz_pixel) * 32 * 32);

	
	kvz_inter_recon_bipred_blend(hi_prec_luma_rec0, hi_prec_luma_rec1, hi_prec_chroma_rec0, hi_prec_chroma_rec1, width, height, xpos, ypos, high_precision_rec0, high_precision_rec1, &result, temp_lcu_y, temp_lcu_u, temp_lcu_v); 
 
 for (temp_y = 0; temp_y < height; ++temp_y) {
  int y_in_lcu = ((ypos + temp_y) & ((LCU_WIDTH)-1));
  for (temp_x = 0; temp_x < width; temp_x += 1) {
   int x_in_lcu = ((xpos + temp_x) & ((LCU_WIDTH)-1));
   printf("%d ", result.rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu]);
  }
 }
 printf("\n");
 
 /*
 for (temp_y = 0; temp_y < height >> 1; ++temp_y) {
  int y_in_lcu = (((ypos >> 1) + temp_y) & (LCU_WIDTH_C - 1));
  for (temp_x = 0; temp_x < width >> 1; ++temp_x) {
   int x_in_lcu = (((xpos >> 1) + temp_x) & (LCU_WIDTH_C - 1));
   printf("%d ", result.rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu]);
  }
 }
 printf("\n");
 */

	for (temp_y = 0; temp_y < height; ++temp_y) {
		int y_in_lcu = ((ypos + temp_y) & ((LCU_WIDTH)-1));
		for (temp_x = 0; temp_x < width; temp_x+=1) {
			int x_in_lcu = ((xpos + temp_x) & ((LCU_WIDTH)-1));
			ASSERT_EQ_FMT(expected_test_result.rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu], result.rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu], "%d");
		}
	}

	for (temp_y = 0; temp_y < height >> 1; ++temp_y) {
		int y_in_lcu = (((ypos >> 1) + temp_y) & (LCU_WIDTH_C - 1));
		for (temp_x = 0; temp_x < width >> 1; ++temp_x) {
			int x_in_lcu = (((xpos >> 1) + temp_x) & (LCU_WIDTH_C - 1));
			ASSERT_EQ_FMT(expected_test_result.rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu], result.rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu], "%d");
			ASSERT_EQ_FMT(expected_test_result.rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu], result.rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu], "%d");
		}
	}
	
	PASS();
}

SUITE(inter_recon_bipred_tests)
{
	setup();

	for (volatile int i = 0; i < strategies.count; ++i) {
		if (strcmp(strategies.strategies[i].type, "inter_recon_bipred") != 0) {
			continue;
		}

		kvz_inter_recon_bipred_blend = strategies.strategies[i].fptr;
		RUN_TEST(test_inter_recon_bipred);
	}
}
