#ifndef CU_H_
#define CU_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/**
 * \ingroup DataStructures
 * \file
 * Coding Unit data structure and related functions.
 */

#include "global.h" // IWYU pragma: keep
#include "image.h"
#include "kvazaar.h"


//Cu stuff
//////////////////////////////////////////////////////////////////////////
// CONSTANTS

typedef enum { CU_NOTSET = 0, CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER } cu_type_t;

typedef enum {
  SIZE_2Nx2N = 0,
  SIZE_2NxN  = 1,
  SIZE_Nx2N  = 2,
  SIZE_NxN   = 3,
  SIZE_2NxnU = 4,
  SIZE_2NxnD = 5,
  SIZE_nLx2N = 6,
  SIZE_nRx2N = 7,
} part_mode_t;

extern const uint8_t kvz_part_mode_num_parts[];
extern const uint8_t kvz_part_mode_offsets[][4][2];
extern const uint8_t kvz_part_mode_sizes[][4][2];

/**
 * \brief Get the x coordinate of a PU.
 *
 * \param part_mode   partition mode of the containing CU
 * \param cu_width    width of the containing CU
 * \param cu_x        x coordinate of the containing CU
 * \param i           number of the PU
 * \return            location of the left edge of the PU
 */
#define PU_GET_X(part_mode, cu_width, cu_x, i) \
  ((cu_x) + kvz_part_mode_offsets[(part_mode)][(i)][0] * (cu_width) / 4)

/**
 * \brief Get the y coordinate of a PU.
 *
 * \param part_mode   partition mode of the containing CU
 * \param cu_width    width of the containing CU
 * \param cu_y        y coordinate of the containing CU
 * \param i           number of the PU
 * \return            location of the top edge of the PU
 */
#define PU_GET_Y(part_mode, cu_width, cu_y, i) \
  ((cu_y) + kvz_part_mode_offsets[(part_mode)][(i)][1] * (cu_width) / 4)

/**
 * \brief Get the width of a PU.
 *
 * \param part_mode   partition mode of the containing CU
 * \param cu_width    width of the containing CU
 * \param i           number of the PU
 * \return            width of the PU
 */
#define PU_GET_W(part_mode, cu_width, i) \
  (kvz_part_mode_sizes[(part_mode)][(i)][0] * (cu_width) / 4)

/**
 * \brief Get the height of a PU.
 *
 * \param part_mode   partition mode of the containing CU
 * \param cu_width    width of the containing CU
 * \param i           number of the PU
 * \return            height of the PU
 */
#define PU_GET_H(part_mode, cu_width, i) \
  (kvz_part_mode_sizes[(part_mode)][(i)][1] * (cu_width) / 4)

//////////////////////////////////////////////////////////////////////////
// TYPES

typedef struct {
  int x;
  int y;
} vector2d_t;

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
  unsigned type      : 3; //!< \brief block type, CU_INTER / CU_INTRA
  unsigned depth     : 3; //!< \brief depth / size of this block
  unsigned part_size : 3; //!< \brief Currently only 2Nx2N, TODO: AMP/SMP/NxN parts
  unsigned tr_depth  : 3; //!< \brief transform depth
  unsigned coded     : 1; //!< \brief flag to indicate this block is coded and reconstructed
  unsigned skipped   : 1; //!< \brief flag to indicate this block is skipped
  unsigned merged    : 1; //!< \brief flag to indicate this block is merged
  unsigned merge_idx : 3; //!< \brief merge index

  cu_cbf_t cbf;
  struct {
    int8_t mode;
    int8_t mode_chroma;
    int8_t tr_skip;    //!< \brief transform skip flag
  } intra[4];
  struct {
    int16_t mv[2][2];  // \brief Motion vectors for L0 and L1
    int16_t mvd[2][2]; // \brief Motion vector differences for L0 and L1
    uint8_t mv_cand[2]; // \brief selected MV candidate
    uint8_t mv_ref[2]; // \brief Index of the encoder_control.ref array.
    uint8_t mv_ref_coded[2]; // \brief Coded and corrected index of ref picture
    uint8_t mv_dir; // \brief Probably describes if mv_ref is L0, L1 or both (bi-pred)
  } inter;
} cu_info_t;

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
  cu_info_t *data;           //!< \brief cu_info data
  int32_t refcount;        //!< \brief number of references in reflists to this cu_array
} cu_array_t;

cu_array_t * kvz_cu_array_alloc(int width_in_scu, int height_in_scu);
int kvz_cu_array_free(cu_array_t *cua);

/**
 * \brief Return the 7 lowest-order bits of the pixel coordinate.
 *
 * The 7 lower-order bits correspond to the distance from the left or top edge
 * of the containing LCU.
 */
#define SUB_SCU(xy) ((xy) & (LCU_WIDTH - 1))

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
  kvz_pixel y[LCU_REF_PX_WIDTH + 1];
  kvz_pixel u[LCU_REF_PX_WIDTH / 2 + 1];
  kvz_pixel v[LCU_REF_PX_WIDTH / 2 + 1];
} lcu_ref_px_t;

typedef struct {
  coeff_t y[LCU_LUMA_SIZE];
  coeff_t u[LCU_CHROMA_SIZE];
  coeff_t v[LCU_CHROMA_SIZE];
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
   *
   \verbatim

      .-- left reference CUs
      v
       0 |  1  2  3  4  5  6  7  8 | 81 <-- top reference CUs
     ----+-------------------------+----
       9 | 10 11 12 13 14 15 16 17 |
      18 | 19 20 21 22 23 24 25 26 <-- this LCU
      27 | 28 29 30 31 32 33 34 35 |
      36 | 37 38 39 40 41 42 43 44 |
      45 | 46 47 48 49 50 51 52 53 |
      54 | 55 56 57 58 59 60 61 62 |
      63 | 64 65 66 67 68 69 70 71 |
      72 | 73 74 75 76 77 78 79 80 |
     ----+-------------------------+----

   \endverbatim
   */
  cu_info_t cu[9*9+1];
} lcu_t;

/**
 * \brief Return pointer to a given CU.
 *
 * \param lcu   pointer to the containing LCU
 * \param x_cu  x-index of the CU
 * \param y_cu  y-index of the CU
 * \return      pointer to the CU
 */
#define LCU_GET_CU(lcu, x_cu, y_cu) \
  (&(lcu)->cu[LCU_CU_OFFSET + (x_cu) + (y_cu) * LCU_T_CU_WIDTH])

/**
 * \brief Return pointer to the top right reference CU.
 */
#define LCU_GET_TOP_RIGHT_CU(lcu) \
  (&(lcu)->cu[LCU_T_CU_WIDTH * LCU_T_CU_WIDTH])

/**
 * \brief Return pointer to the CU containing a given pixel.
 *
 * \param lcu   pointer to the containing LCU
 * \param x_px  x-coordinate relative to the upper left corner of the LCU
 * \param y_px  y-coordinate relative to the upper left corner of the LCU
 * \return      pointer to the CU at coordinates (x_px, y_px)
 */
#define LCU_GET_CU_AT_PX(lcu, x_px, y_px) LCU_GET_CU(lcu, (x_px) >> 3, (y_px) >> 3)

/**
 * \brief Return pointer to a CU relative to the given CU.
 *
 * \param cu      pointer to a CU in the array at some location (x, y)
 * \param x_offs  x-offset
 * \param y_offs  y-offset
 * \return        pointer to the CU at (x + x_offs, y + y_offs)
 */
#define CU_GET_CU(cu_array, x_offs, y_offs) \
  (&cu_array[(x_offs) + (y_offs) * LCU_T_CU_WIDTH])

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


void kvz_coefficients_blit(const coeff_t *orig, coeff_t *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);

unsigned kvz_coefficients_calc_abs(const coeff_t *const buf, const int buf_stride,
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
