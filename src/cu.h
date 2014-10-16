#ifndef CU_H_
#define CU_H_
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
 * \brief CU and coefficients related functions
 */

#include "global.h"
#include "image.h"

//Cu stuff
//////////////////////////////////////////////////////////////////////////
// CONSTANTS

typedef enum { CU_NOTSET = 0, CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER } cu_type_t;

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
  int8_t mode;
  int8_t mode_chroma;
  int8_t tr_skip;    //!< \brief transform skip flag
} cu_info_intra;

/**
 * \brief Struct for CU inter info
 */
typedef struct
{
  double cost;
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

#define CHECKPOINT_CU(prefix_str, cu) CHECKPOINT(prefix_str " type=%d depth=%d part_size=%d tr_depth=%d coded=%d " \
  "skipped=%d merged=%d merge_idx=%d cbf.y=%d cbf.u=%d cbf.v=%d " \
  "intra[0].cost=%u intra[0].bitcost=%u intra[0].mode=%d intra[0].mode_chroma=%d intra[0].tr_skip=%d " \
  "intra[1].cost=%u intra[1].bitcost=%u intra[1].mode=%d intra[1].mode_chroma=%d intra[1].tr_skip=%d " \
  "intra[2].cost=%u intra[2].bitcost=%u intra[2].mode=%d intra[2].mode_chroma=%d intra[2].tr_skip=%d " \
  "intra[3].cost=%u intra[3].bitcost=%u intra[3].mode=%d intra[3].mode_chroma=%d intra[3].tr_skip=%d " \
  "inter.cost=%u inter.bitcost=%u inter.mv[0]=%d inter.mv[1]=%d inter.mvd[0]=%d inter.mvd[1]=%d " \
  "inter.mv_cand=%d inter.mv_ref=%d inter.mv_dir=%d inter.mode=%d" \
  , (cu).type, (cu).depth, (cu).part_size, (cu).tr_depth, (cu).coded, \
  (cu).skipped, (cu).merged, (cu).merge_idx, (cu).cbf.y, (cu).cbf.u, (cu).cbf.v, \
  (cu).intra[0].cost, (cu).intra[0].bitcost, (cu).intra[0].mode, (cu).intra[0].mode_chroma, (cu).intra[0].tr_skip, \
  (cu).intra[1].cost, (cu).intra[1].bitcost, (cu).intra[1].mode, (cu).intra[1].mode_chroma, (cu).intra[1].tr_skip, \
  (cu).intra[2].cost, (cu).intra[2].bitcost, (cu).intra[2].mode, (cu).intra[2].mode_chroma, (cu).intra[2].tr_skip, \
  (cu).intra[3].cost, (cu).intra[3].bitcost, (cu).intra[3].mode, (cu).intra[3].mode_chroma, (cu).intra[3].tr_skip, \
  (cu).inter.cost, (cu).inter.bitcost, (cu).inter.mv[0], (cu).inter.mv[1], (cu).inter.mvd[0], (cu).inter.mvd[1], \
  (cu).inter.mv_cand, (cu).inter.mv_ref, (cu).inter.mv_dir, (cu).inter.mode)

typedef struct {
  cu_info *data;           //!< \brief cu_info data
  int32_t refcount;        //!< \brief number of references in reflists to this cu_array
} cu_array;

cu_array * cu_array_alloc(int width_in_scu, int height_in_scu);
int cu_array_free(cu_array *cua);
  

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

#define CHECKPOINT_LCU(prefix_str, lcu) do { \
  CHECKPOINT_CU(prefix_str " cu[0]", (lcu).cu[0]); \
  CHECKPOINT_CU(prefix_str " cu[1]", (lcu).cu[1]); \
  CHECKPOINT_CU(prefix_str " cu[2]", (lcu).cu[2]); \
  CHECKPOINT_CU(prefix_str " cu[3]", (lcu).cu[3]); \
  CHECKPOINT_CU(prefix_str " cu[4]", (lcu).cu[4]); \
  CHECKPOINT_CU(prefix_str " cu[5]", (lcu).cu[5]); \
  CHECKPOINT_CU(prefix_str " cu[6]", (lcu).cu[6]); \
  CHECKPOINT_CU(prefix_str " cu[7]", (lcu).cu[7]); \
  CHECKPOINT_CU(prefix_str " cu[8]", (lcu).cu[8]); \
  CHECKPOINT_CU(prefix_str " cu[9]", (lcu).cu[9]); \
  CHECKPOINT_CU(prefix_str " cu[10]", (lcu).cu[10]); \
  CHECKPOINT_CU(prefix_str " cu[11]", (lcu).cu[11]); \
  CHECKPOINT_CU(prefix_str " cu[12]", (lcu).cu[12]); \
  CHECKPOINT_CU(prefix_str " cu[13]", (lcu).cu[13]); \
  CHECKPOINT_CU(prefix_str " cu[14]", (lcu).cu[14]); \
  CHECKPOINT_CU(prefix_str " cu[15]", (lcu).cu[15]); \
  CHECKPOINT_CU(prefix_str " cu[16]", (lcu).cu[16]); \
  CHECKPOINT_CU(prefix_str " cu[17]", (lcu).cu[17]); \
  CHECKPOINT_CU(prefix_str " cu[18]", (lcu).cu[18]); \
  CHECKPOINT_CU(prefix_str " cu[19]", (lcu).cu[19]); \
  CHECKPOINT_CU(prefix_str " cu[20]", (lcu).cu[20]); \
  CHECKPOINT_CU(prefix_str " cu[21]", (lcu).cu[21]); \
  CHECKPOINT_CU(prefix_str " cu[22]", (lcu).cu[22]); \
  CHECKPOINT_CU(prefix_str " cu[23]", (lcu).cu[23]); \
  CHECKPOINT_CU(prefix_str " cu[24]", (lcu).cu[24]); \
  CHECKPOINT_CU(prefix_str " cu[25]", (lcu).cu[25]); \
  CHECKPOINT_CU(prefix_str " cu[26]", (lcu).cu[26]); \
  CHECKPOINT_CU(prefix_str " cu[27]", (lcu).cu[27]); \
  CHECKPOINT_CU(prefix_str " cu[28]", (lcu).cu[28]); \
  CHECKPOINT_CU(prefix_str " cu[29]", (lcu).cu[29]); \
  CHECKPOINT_CU(prefix_str " cu[30]", (lcu).cu[30]); \
  CHECKPOINT_CU(prefix_str " cu[31]", (lcu).cu[31]); \
  CHECKPOINT_CU(prefix_str " cu[32]", (lcu).cu[32]); \
  CHECKPOINT_CU(prefix_str " cu[33]", (lcu).cu[33]); \
  CHECKPOINT_CU(prefix_str " cu[34]", (lcu).cu[34]); \
  CHECKPOINT_CU(prefix_str " cu[35]", (lcu).cu[35]); \
  CHECKPOINT_CU(prefix_str " cu[36]", (lcu).cu[36]); \
  CHECKPOINT_CU(prefix_str " cu[37]", (lcu).cu[37]); \
  CHECKPOINT_CU(prefix_str " cu[38]", (lcu).cu[38]); \
  CHECKPOINT_CU(prefix_str " cu[39]", (lcu).cu[39]); \
  CHECKPOINT_CU(prefix_str " cu[40]", (lcu).cu[40]); \
  CHECKPOINT_CU(prefix_str " cu[41]", (lcu).cu[41]); \
  CHECKPOINT_CU(prefix_str " cu[42]", (lcu).cu[42]); \
  CHECKPOINT_CU(prefix_str " cu[43]", (lcu).cu[43]); \
  CHECKPOINT_CU(prefix_str " cu[44]", (lcu).cu[44]); \
  CHECKPOINT_CU(prefix_str " cu[45]", (lcu).cu[45]); \
  CHECKPOINT_CU(prefix_str " cu[46]", (lcu).cu[46]); \
  CHECKPOINT_CU(prefix_str " cu[47]", (lcu).cu[47]); \
  CHECKPOINT_CU(prefix_str " cu[48]", (lcu).cu[48]); \
  CHECKPOINT_CU(prefix_str " cu[49]", (lcu).cu[49]); \
  CHECKPOINT_CU(prefix_str " cu[50]", (lcu).cu[50]); \
  CHECKPOINT_CU(prefix_str " cu[51]", (lcu).cu[51]); \
  CHECKPOINT_CU(prefix_str " cu[52]", (lcu).cu[52]); \
  CHECKPOINT_CU(prefix_str " cu[53]", (lcu).cu[53]); \
  CHECKPOINT_CU(prefix_str " cu[54]", (lcu).cu[54]); \
  CHECKPOINT_CU(prefix_str " cu[55]", (lcu).cu[55]); \
  CHECKPOINT_CU(prefix_str " cu[56]", (lcu).cu[56]); \
  CHECKPOINT_CU(prefix_str " cu[57]", (lcu).cu[57]); \
  CHECKPOINT_CU(prefix_str " cu[58]", (lcu).cu[58]); \
  CHECKPOINT_CU(prefix_str " cu[59]", (lcu).cu[59]); \
  CHECKPOINT_CU(prefix_str " cu[60]", (lcu).cu[60]); \
  CHECKPOINT_CU(prefix_str " cu[61]", (lcu).cu[61]); \
  CHECKPOINT_CU(prefix_str " cu[62]", (lcu).cu[62]); \
  CHECKPOINT_CU(prefix_str " cu[63]", (lcu).cu[63]); \
  CHECKPOINT_CU(prefix_str " cu[64]", (lcu).cu[64]); \
  CHECKPOINT_CU(prefix_str " cu[65]", (lcu).cu[65]); \
  CHECKPOINT_CU(prefix_str " cu[66]", (lcu).cu[66]); \
  CHECKPOINT_CU(prefix_str " cu[67]", (lcu).cu[67]); \
  CHECKPOINT_CU(prefix_str " cu[68]", (lcu).cu[68]); \
  CHECKPOINT_CU(prefix_str " cu[69]", (lcu).cu[69]); \
  CHECKPOINT_CU(prefix_str " cu[70]", (lcu).cu[70]); \
  CHECKPOINT_CU(prefix_str " cu[71]", (lcu).cu[71]); \
  CHECKPOINT_CU(prefix_str " cu[72]", (lcu).cu[72]); \
  CHECKPOINT_CU(prefix_str " cu[73]", (lcu).cu[73]); \
  CHECKPOINT_CU(prefix_str " cu[74]", (lcu).cu[74]); \
  CHECKPOINT_CU(prefix_str " cu[75]", (lcu).cu[75]); \
  CHECKPOINT_CU(prefix_str " cu[76]", (lcu).cu[76]); \
  CHECKPOINT_CU(prefix_str " cu[77]", (lcu).cu[77]); \
  CHECKPOINT_CU(prefix_str " cu[78]", (lcu).cu[78]); \
  CHECKPOINT_CU(prefix_str " cu[79]", (lcu).cu[79]); \
  CHECKPOINT_CU(prefix_str " cu[80]", (lcu).cu[80]); \
  CHECKPOINT_CU(prefix_str " cu[81]", (lcu).cu[81]); \
} while(0)


void coefficients_blit(const coefficient *orig, coefficient *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);

unsigned coefficients_calc_abs(const coefficient *const buf, const int buf_stride,
                        const int width);



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

#define GET_SPLITDATA(CU,curDepth) ((CU)->depth > curDepth)
#define SET_SPLITDATA(CU,flag) { (CU)->split=(flag); }

#endif
