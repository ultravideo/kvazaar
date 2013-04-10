/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file picture.h
    \brief Picture header
    \author Marko Viitanen
    \date 2013-04
    
    Contains all picture related functions and structs
*/

#ifndef _PICTURE_H_
#define _PICTURE_H_


/* Functions */
uint32_t SAD64x64(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD32x32(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD16x16(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD8x8(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
uint32_t SAD4x4(int16_t *block,uint32_t stride1,int16_t* block2, uint32_t stride2);
double imagePSNR(uint8_t *frame1, uint8_t *frame2, int32_t x, int32_t y);

/** \defgroup picture_group Picture handler group
 *  This group contains all picture related stuff
 *  @{
 */

enum { CU_NOTSET = 0,CU_PCM, CU_SKIP, CU_SPLIT, CU_INTRA, CU_INTER };

#define GET_SPLITDATA(CU) ((CU)->split)
#define SET_SPLITDATA(CU,flag) { (CU)->split=(flag); }

/*!
    \brief Struct for CU intra info
*/
typedef struct
{
  uint8_t mode;
  uint32_t cost;
} CU_info_intra;

/*!
    \brief Struct for CU inter info
*/
typedef struct
{
  uint8_t mode;
  uint32_t cost;
  int16_t mv[2];
} CU_info_inter;


/*!
    \brief Struct for CU info
*/
typedef struct
{
  uint8_t type;
  int8_t coded;
  CU_info_intra intra;
  CU_info_inter inter;
  uint8_t split;
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

  int width;          /*!< \brief Picture width */
  int height;         /*!< \brief Picture height  */
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

int picture_destroy(picture *pic);


enum { SLICE_P = 0, SLICE_B = 1, SLICE_I = 2 };


  /** @} */ // end of group1

#endif