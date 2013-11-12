#ifndef SAO_H_
#define SAO_H_
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
#include "picture.h"


typedef enum { COLOR_Y = 0, COLOR_U = 1, COLOR_V = 2, NUM_COLORS } color_index;
typedef enum { SAO_TYPE_NONE = 0, SAO_TYPE_BAND, SAO_TYPE_EDGE } sao_type;
typedef enum { SAO_EO0 = 0, SAO_EO1, SAO_EO2, SAO_EO3, SAO_NUM_EO } sao_eo_class;
typedef enum { SAO_EO_CAT0 = 0, SAO_EO_CAT1, SAO_EO_CAT2, SAO_EO_CAT3, SAO_EO_CAT4, NUM_SAO_EDGE_CATEGORIES } sao_eo_cat;

// Offsets of a and b in relation to c.
// dir_offset[dir][a or b]
// |       |   a   | a     |     a |
// | a c b |   c   |   c   |   c   |
// |       |   b   |     b | b     |
static const vector2d g_sao_edge_offsets[SAO_NUM_EO][2] = { 
  { { -1, 0 }, { 1, 0 } },
  { { 0, -1 }, { 0, 1 } },
  { { -1, -1 }, { 1, 1 } },
  { { 1, -1 }, { -1, 1 } }
};


typedef struct sao_info_struct {
  sao_type type;
  sao_eo_class eo_class;
  int ddistortion;
  int merge_left_flag;
  int merge_up_flag;
  int offsets[NUM_SAO_EDGE_CATEGORIES];
} sao_info;


void init_sao_info(sao_info *sao);
void sao_search_chroma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao);
void sao_search_luma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao);
void sao_reconstruct(picture *pic, const pixel *old_rec, 
                     unsigned x_ctb, unsigned y_ctb, 
                     const sao_info *sao, color_index color_i);

#endif