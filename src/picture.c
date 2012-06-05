/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file picture.c
    \brief Functions to handle pictures
    \author Marko Viitanen
    \date 2012-06
    
  This file contains all the needed functions to handle pictures

*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "picture.h"

/** \defgroup picture_group Picture handler group
 *  This group contains all picture related stuff
 *  @{
 */


/*!
    \brief Allocate memory for picture_list
    \param size initial array size
    \return picture_list pointer, NULL on failure
*/
  picture_list *picture_list_init(int size)
  {
    picture_list *list = (picture_list *)malloc(sizeof(picture_list));
    list->size = size;
    if(size > 0)
    {
      list->pics = (picture**)malloc(sizeof(picture*)*size);
    }

    list->used_size = 0;

    return list;
  }

  /*!
    \brief Resize picture_list array
    \param list picture_list pointer
    \param size new array size
    \return 1 on success, 0 on failure
  */
  int picture_list_resize(picture_list *list, int size)
  {
    unsigned int i;
    picture** old_pics = NULL;

    //No need to do anything when resizing to same size
    if(size == list->size)
    {
      return 1;
    }

    //Save the old list
    if(list->used_size > 0)
    {
      old_pics = list->pics;
    }

    //allocate space for the new list
    list->pics = (picture**)malloc(sizeof(picture*)*size);

    //Copy everthing from the old list to the new if needed.
    if(old_pics != NULL)
    {
      for(i = 0; i < list->used_size; i++)
      {
        list->pics[i] = old_pics[i];
      }

      free(old_pics);
    }

    return 1;
  }

  /*!
    \brief Free memory allocated to the picture_list
    \param list picture_list pointer
    \return 1 on success, 0 on failure
  */
  int picture_list_destroy(picture_list *list)
  {
    unsigned int i;
    if(list->used_size > 0)
    {
      for(i = 0; i < list->used_size; i++)
      {
        picture_destroy(list->pics[i]);
      }
    }
    
    if(list->size > 0)
    {
      free(list->pics);
    }
    free(list);
    return 1;
  }

  /*!
    \brief Free memory allocated to picture
    \param pic picture pointer
    \return 1 on success, 0 on failure
  */
  int picture_destroy(picture *pic)
  {
    free(pic->uData);
    free(pic->vData);
    free(pic->yData);
    return 1;
  }


  /** @} */ // end of group1