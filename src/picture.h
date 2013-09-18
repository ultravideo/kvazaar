/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file picture.h
    \brief Picture header
    \author Marko Viitanen
    \date 2013-04
    
    Contains all picture related functions and structs
*/

#ifndef _PICTURE_H_
#define _PICTURE_H_

#include "global.h"


/* Functions */
uint32_t SAD64x64(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD32x32(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD16x16(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD8x8(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD4x4(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD(uint8_t *data1, uint8_t *data2, unsigned width, unsigned height, unsigned stride);

double imagePSNR(uint8_t *frame1, uint8_t *frame2, int32_t x, int32_t y);

/** \defgroup picture_group Picture handler group
 *  This group contains all picture related stuff
 *  @{
 */

enum { CU_NOTSET = 0,CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER };

#define GET_SPLITDATA(CU,curDepth) ((CU)->depth > curDepth)
#define SET_SPLITDATA(CU,flag) { (CU)->split=(flag); }

/*!
    \brief Struct for CU intra info
*/
typedef struct
{
  int8_t mode;
  uint32_t cost;
} CU_info_intra;

/*!
    \brief Struct for CU inter info
*/
typedef struct
{
  int8_t mode;
  uint32_t cost;
  int16_t mv[2];
  uint8_t mv_ref; // \brief Index of the encoder_control.ref array.
  uint8_t mv_dir; // \brief Probably describes if mv_ref is forward, backward or both. Might not be needed?
} CU_info_inter;


/*!
    \brief Struct for CU info
*/
typedef struct
{  
  int8_t type;
  int8_t depth;
  int8_t coded;
  CU_info_intra intra;
  CU_info_inter inter;
  int8_t split;
} CU_info;

/*!
    \brief Struct which contains all picture data
*/
typedef struct
{
  uint8_t* yData;     /*!< \brief Pointer to Y-data  */
  uint8_t* uData;     /*!< \brief Pointer to U-data  */
  uint8_t* vData;     /*!< \brief Pointer to V-data  */

  uint8_t* yRecData;     /*!< \brief Pointer to reconstructed Y-data  */
  uint8_t* uRecData;     /*!< \brief Pointer to reconstructed U-data  */
  uint8_t* vRecData;     /*!< \brief Pointer to reconstructed V-data  */

  int32_t width;          /*!< \brief Picture width */
  int32_t height;         /*!< \brief Picture height  */
  int32_t height_in_lcu;  /*!< \brief input picture width in number of LCU's */
  int32_t width_in_lcu;   /*!< \brief input picture height in number of LCU's */
  uint8_t referenced; /*!< \brief Is this picture referenced */
  CU_info** CU;     /*!< \brief info for each CU at each depth */
  uint8_t type;
  uint8_t slicetype;
} picture;

/*!
    \brief Struct which contains array of picture structs
*/
typedef struct
{
  picture** pics; /*!< \brief Pointer to array of picture pointers  */
  unsigned int size;       /*!< \brief Array size */
  unsigned int used_size;

} picture_list;


picture_list *picture_list_init(int size);
int picture_list_resize(picture_list *list, int size);
int picture_list_destroy(picture_list *list);
int picture_list_add(picture_list *list,picture* pic);
int picture_list_rem(picture_list *list,int n, int8_t destroy);

picture *picture_init(int32_t width, int32_t height, int32_t width_in_LCU, int32_t height_in_LCU);
int picture_destroy(picture *pic);


void picture_setBlockCoded(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t coded);
void picture_setBlockSplit(picture* pic,uint32_t xCtb, uint32_t yCtb, uint8_t depth, int8_t split);

enum { SLICE_B = 0, SLICE_P = 1, SLICE_I = 2 };

enum { REF_PIC_LIST_0 = 0,REF_PIC_LIST_1 = 1,REF_PIC_LIST_X = 100 };


  /** @} */ // end of group1

#endif