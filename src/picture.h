#ifndef PICTURE_H_
#define PICTURE_H_
/**
 * \file
 * \brief Coding Unit (CU) and picture data related functions.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"


//////////////////////////////////////////////////////////////////////////
// CONSTANTS

enum { CU_NOTSET = 0, CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER };
enum { SLICE_B = 0, SLICE_P = 1, SLICE_I = 2 };
enum { REF_PIC_LIST_0 = 0, REF_PIC_LIST_1 = 1, REF_PIC_LIST_X = 100 };


//////////////////////////////////////////////////////////////////////////
// TYPES

/**
 * \brief Struct for CU intra info
 */
typedef struct
{
  int8_t mode;
  int8_t mode_chroma;
  uint32_t cost;
} cu_info_intra;

/**
 * \brief Struct for CU inter info
 */
typedef struct
{
  int8_t mode;
  uint32_t cost;

  int16_t mv[2];
  int16_t mvd[2];
  uint8_t mv_ref; // \brief Index of the encoder_control.ref array.
  uint8_t mv_dir; // \brief Probably describes if mv_ref is forward, backward or both. Might not be needed?
} cu_info_inter;

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
  int8_t coeff_y;    //!< \brief is there coded coeffs Y
  int8_t coeff_u;    //!< \brief is there coded coeffs U
  int8_t coeff_v;    //!< \brief is there coded coeffs V

  int8_t coeff_top_y[MAX_DEPTH+1];  //!< \brief is there coded coeffs Y in top level
  int8_t coeff_top_u[MAX_DEPTH+1];  //!< \brief is there coded coeffs U in top level
  int8_t coeff_top_v[MAX_DEPTH+1];  //!< \brief is there coded coeffs V in top level
  cu_info_intra intra;
  cu_info_inter inter;
} cu_info;

/**
 * \brief Struct which contains all picture data
 */
typedef struct
{
  pixel* y_data;        //!< \brief Pointer to luma pixel array.
  pixel* u_data;        //!< \brief Pointer to chroma U pixel array.
  pixel* v_data;        //!< \brief Pointer to chroma V pixel array.

  pixel* y_recdata;     //!< \brief Pointer to reconstructed Y-data.
  pixel* u_recdata;     //!< \brief Pointer to reconstructed U-data.
  pixel* v_recdata;     //!< \brief Pointer to reconstructed V-data.

  pixel* pred_y;        //!< \brief Pointer to predicted Y
  pixel* pred_u;        //!< \brief Pointer to predicted U
  pixel* pred_v;        //!< \brief Pointer to predicted V

  coefficient* coeff_y;   //!< \brief coefficient pointer Y
  coefficient* coeff_u;   //!< \brief coefficient pointer U
  coefficient* coeff_v;   //!< \brief coefficient pointer V

  int32_t width;          //!< \brief Luma pixel array width.
  int32_t height;         //!< \brief Luma pixel array height.
  int32_t height_in_lcu;  //!< \brief Picture width in number of LCU's.
  int32_t width_in_lcu;   //!< \brief Picture height in number of LCU's.
  uint8_t referenced;     //!< \brief Whether this picture is referenced.
  cu_info** cu_array;     //!< \brief Info for each CU at each depth.
  uint8_t type;
  uint8_t slicetype;
} picture;

/**
 * \brief Struct which contains array of picture structs
 */
typedef struct
{
  picture** pics;          //!< \brief Pointer to array of picture pointers.
  unsigned int size;       //!< \brief Array size.
  unsigned int used_size;
} picture_list;


//////////////////////////////////////////////////////////////////////////
// FUNCTIONS

picture * picture_init(int32_t width, int32_t height, 
                       int32_t width_in_lcu, int32_t height_in_lcu);
int picture_destroy(picture *pic);
void picture_set_block_coded(picture *pic, uint32_t x_scu, uint32_t y_scu,
                             uint8_t depth, int8_t coded);
void picture_set_block_residual(picture *pic, uint32_t x_scu, uint32_t y_scu,
                                uint8_t depth, int8_t residual);
void picture_set_block_split(picture *pic, uint32_t x_scu, uint32_t y_scu,
                             uint8_t depth, int8_t split);
void picture_set_block_skipped(picture *pic, uint32_t x_scu, uint32_t y_scu,
                                uint8_t depth, int8_t skipped);
picture_list * picture_list_init(int size);
int picture_list_resize(picture_list *list, int size);
int picture_list_destroy(picture_list *list);
int picture_list_add(picture_list *list, picture *pic);
int picture_list_rem(picture_list *list, int n, int8_t destroy);

typedef unsigned (*cost_16bit_nxn_func)(pixel *block1, pixel *block2);


cost_16bit_nxn_func get_satd_16bit_nxn_func(unsigned n);
cost_16bit_nxn_func get_sad_16bit_nxn_func(unsigned n);

unsigned satd_16bit_nxn(pixel *block1, pixel *block2, unsigned n);
unsigned sad_16bit_nxn(pixel *block1, pixel *block2, unsigned n);

unsigned calc_sad(const picture *pic, const picture *ref, 
                  int pic_x, int pic_y, int ref_x, int ref_y, 
                  int block_width, int block_height);

double image_psnr(pixel *frame1, pixel *frame2, int32_t x, int32_t y);


//////////////////////////////////////////////////////////////////////////
// MACROS

#define GET_SPLITDATA(CU,curDepth) ((CU)->depth > curDepth)
#define SET_SPLITDATA(CU,flag) { (CU)->split=(flag); }


#endif
