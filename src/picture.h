#ifndef PICTURE_H_
#define PICTURE_H_
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
 * \brief Coding Unit (CU) and picture data related functions.
 */

#include "global.h"

//#include "sao.h"
struct sao_info_struct;


//////////////////////////////////////////////////////////////////////////
// CONSTANTS

typedef enum { CU_NOTSET = 0, CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER } cu_type_t;
enum { SLICE_B = 0, SLICE_P = 1, SLICE_I = 2 };
enum { REF_PIC_LIST_0 = 0, REF_PIC_LIST_1 = 1, REF_PIC_LIST_X = 100 };
typedef enum { COLOR_Y = 0, COLOR_U, COLOR_V, NUM_COLORS } color_index;


//////////////////////////////////////////////////////////////////////////
// TYPES

typedef struct {
  int x;
  int y;
} vector2d;

/**
 * \brief Struct for CU intra info
 */
typedef struct
{
  uint32_t cost;
  uint32_t bitcost;
  int8_t mode;
  int8_t mode_chroma;
  int8_t tr_skip;    //!< \brief transform skip flag
} cu_info_intra;

/**
 * \brief Struct for CU inter info
 */
typedef struct
{
  uint32_t cost;
  uint32_t bitcost;
  int16_t mv[2];
  int16_t mvd[2];
  uint8_t mv_cand; // \brief selected MV candidate
  uint8_t mv_ref; // \brief Index of the encoder_control.ref array.
  uint8_t mv_dir; // \brief Probably describes if mv_ref is forward, backward or both. Might not be needed?
  int8_t mode;
} cu_info_inter;

typedef struct
{
  uint8_t y;
  uint8_t u;
  uint8_t v;
} cu_cbf_t;

/**
 * \brief Struct for CU info
 */
typedef struct
{
  int8_t type;       //!< \brief block type, CU_INTER / CU_INTRA
  int8_t depth;      //!< \brief depth / size of this block
  int8_t part_size;  //!< \brief Currently only 2Nx2N, TODO: AMP/SMP/NxN parts
  int8_t tr_depth;   //!< \brief transform depth
  int8_t coded;      //!< \brief flag to indicate this block is coded and reconstructed
  int8_t skipped;    //!< \brief flag to indicate this block is skipped
  int8_t merged;     //!< \brief flag to indicate this block is merged
  int8_t merge_idx;  //!< \brief merge index

  cu_cbf_t cbf;
  cu_info_intra intra[4];
  cu_info_inter inter;
} cu_info;

/**
 * \brief Struct which contains all picture data
 */
typedef struct picture
{
  pixel* y_data;        //!< \brief Pointer to luma pixel array.
  pixel* u_data;        //!< \brief Pointer to chroma U pixel array.
  pixel* v_data;        //!< \brief Pointer to chroma V pixel array.
  pixel *data[NUM_COLORS]; // \brief  Alternate access method to same data.

  pixel* y_recdata;     //!< \brief Pointer to reconstructed Y-data.
  pixel* u_recdata;     //!< \brief Pointer to reconstructed U-data.
  pixel* v_recdata;     //!< \brief Pointer to reconstructed V-data.
  pixel *recdata[NUM_COLORS]; // \brief Alternate access method to same data.

  coefficient* coeff_y;   //!< \brief coefficient pointer Y
  coefficient* coeff_u;   //!< \brief coefficient pointer U
  coefficient* coeff_v;   //!< \brief coefficient pointer V

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.
  uint8_t referenced;     //!< \brief Whether this picture is referenced.
  int32_t refcount;     //!< \brief Number of references in reflist to the picture
  cu_info* cu_array;     //!< \brief Info for each CU at each depth.
  struct sao_info_struct *sao_luma;   //!< \brief Array of sao parameters for every LCU.
  struct sao_info_struct *sao_chroma;   //!< \brief Array of sao parameters for every LCU.
  int32_t poc;           //!< \brief Picture order count
} picture;


#define SUB_SCU_BIT_MASK (64 - 1)
#define SUB_SCU(xy) (xy & SUB_SCU_BIT_MASK)
#define LCU_CU_WIDTH 8
#define LCU_T_CU_WIDTH 9
#define LCU_CU_OFFSET 10

// Width from top left of the LCU, so +1 for ref buffer size.
#define LCU_REF_PX_WIDTH (LCU_WIDTH + LCU_WIDTH / 2)

/**
 * Top and left intra reference pixels for LCU.
 * - Intra needs maximum of 32 to the right and down from LCU border.
 * - First pixel is the top-left pixel.
 */
typedef struct {
  pixel y[LCU_REF_PX_WIDTH + 1];
  pixel u[LCU_REF_PX_WIDTH / 2 + 1];
  pixel v[LCU_REF_PX_WIDTH / 2 + 1];
} lcu_ref_px_t;

typedef struct {
  coefficient y[LCU_LUMA_SIZE];
  coefficient u[LCU_CHROMA_SIZE];
  coefficient v[LCU_CHROMA_SIZE];
} lcu_coeff_t;

typedef struct {
  pixel y[LCU_LUMA_SIZE];
  pixel u[LCU_CHROMA_SIZE];
  pixel v[LCU_CHROMA_SIZE];
} lcu_yuv_t;

typedef struct {
  int size;
  pixel *y;
  pixel *u;
  pixel *v;
} yuv_t;

typedef struct {
  lcu_ref_px_t top_ref;  //!< Reference pixels from adjacent LCUs.
  lcu_ref_px_t left_ref; //!< Reference pixels from adjacent LCUs.
  lcu_yuv_t ref; //!< LCU reference pixels
  lcu_yuv_t rec; //!< LCU reconstructed pixels
  /**
   * We get the coefficients as a byproduct of doing reconstruction during the
   * search. It might be more efficient to recalculate the final coefficients
   * once we know the final modes rather than copying them.
   */
  lcu_coeff_t coeff; //!< LCU coefficients

  /**
   * A 9x9 CU array for the LCU, +1 CU.
   * - Top reference CUs on row 0.
   * - Left reference CUs on column 0.
   * - All of LCUs CUs on 1:9, 1:9.
   * - Top right reference CU on the last slot.
   */
  cu_info cu[9*9+1];
} lcu_t;

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS

/**
 * Check if CBF in a given level >= depth is true.
 */
static INLINE int cbf_is_set(uint8_t cbf_flags, int depth)
{
  // Transform data for 4x4 blocks is stored at depths 4-8 for luma, so masks
  // for those levels don't include the other ones.
  static const uint8_t masks[8] = { 0xff, 0x7f, 0x3f, 0x1f, 0x8, 0x4, 0x2, 0x1 };

  return (cbf_flags & masks[depth]) != 0;
}

/**
 * Set CBF in a level to true.
 */
static INLINE void cbf_set(uint8_t *cbf_flags, int depth)
{
  // Return value of the bit corresponding to the level.
  *cbf_flags |= 1 << (7 - depth);
}

/**
 * Set CBF in a levels <= depth to false.
 */
static INLINE void cbf_clear(uint8_t *cbf_flags, int depth)
{
  static const uint8_t masks[8] = { 0xff, 0x7f, 0x3f, 0x1f, 0x8, 0x4, 0x2, 0x1 };

  *cbf_flags &= ~masks[depth];
}

yuv_t * yuv_t_alloc(int luma_size);
void yuv_t_free(yuv_t * yuv);

picture * picture_alloc(int32_t width, int32_t height,
                       int32_t width_in_lcu, int32_t height_in_lcu);
int picture_free(picture *pic);

void picture_blit_pixels(const pixel* orig, pixel *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);
void picture_blit_coeffs(const coefficient *orig, coefficient *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);


typedef unsigned (*cost_16bit_nxn_func)(const pixel *block1, const pixel *block2);


cost_16bit_nxn_func get_satd_16bit_nxn_func(unsigned n);
cost_16bit_nxn_func get_sad_16bit_nxn_func(unsigned n);

unsigned satd_16bit_nxn(pixel *block1, pixel *block2, unsigned n);
unsigned sad_16bit_nxn(pixel *block1, pixel *block2, unsigned n);

unsigned calc_sad(const picture *pic, const picture *ref,
                  int pic_x, int pic_y, int ref_x, int ref_y,
                  int block_width, int block_height);

unsigned calc_ssd(const pixel *const ref, const pixel *const rec,
                  const int ref_stride, const int rec_stride,
                  const int width);
unsigned calc_abs_coeff(const coefficient *const buf, const int buf_stride,
                        const int width);

double image_psnr(pixel *frame1, pixel *frame2, int32_t x, int32_t y);


//////////////////////////////////////////////////////////////////////////
// MACROS

#define GET_SPLITDATA(CU,curDepth) ((CU)->depth > curDepth)
#define SET_SPLITDATA(CU,flag) { (CU)->split=(flag); }


#endif
