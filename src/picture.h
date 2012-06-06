/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file picture.h
    \brief Picture header
    \author Marko Viitanen
    \date 2012-06
    
    Contains all picture related functions and structs
*/

#ifndef _PICTURE_H_
#define _PICTURE_H_


/** \defgroup picture_group Picture handler group
 *  This group contains all picture related stuff
 *  @{
 */

/*!
    \brief Struct which contains all picture data
*/
typedef struct
{
  uint8_t* yData;     /*!< \brief Pointer to Y-data  */
  uint8_t* uData;     /*!< \brief Pointer to U-data  */
  uint8_t* vData;     /*!< \brief Pointer to V-data  */
  int width;          /*!< \brief Picture width */
  int height;         /*!< \brief Picture height  */
  uint8_t referenced; /*!< \brief Is this picture referenced */
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