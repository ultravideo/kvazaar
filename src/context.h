#ifndef CONTEXT_H_
#define CONTEXT_H_
/**
 * \file
 * \brief Context derivation for CABAC.

 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "encoder.h"


// Types
typedef struct
{
  uint8_t  uc_state;
  uint32_t bins_coded;
} cabac_ctx;


// Functions
void ctx_init(cabac_ctx* ctx, uint32_t qp, uint32_t init_value);
void init_contexts(encoder_control *encoder, int8_t slice);
int32_t context_calc_pattern_sig_ctx( const uint32_t *sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width);

uint32_t context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,uint32_t pos_x, uint32_t pos_y,int32_t width);


int32_t context_get_sig_ctx_inc(int32_t pattern_sig_ctx,uint32_t scan_idx,int32_t pos_x,
                                int32_t pos_y,int32_t block_type,int32_t width, int8_t texture_type);

// CONTEXTS
extern cabac_ctx g_sao_merge_flag_model;
extern cabac_ctx g_sao_type_idx_model;
extern cabac_ctx g_split_flag_model[3];
extern cabac_ctx g_intra_mode_model;
extern cabac_ctx g_chroma_pred_model[2];
extern cabac_ctx g_trans_subdiv_model[3];
extern cabac_ctx g_qt_cbf_model_luma[3];
extern cabac_ctx g_qt_cbf_model_chroma[3];
extern cabac_ctx g_part_size_model[4];
extern cabac_ctx g_cu_sig_coeff_group_model[4];
extern cabac_ctx g_cu_sig_model_luma[27];
extern cabac_ctx g_cu_sig_model_chroma[15];
extern cabac_ctx g_cu_ctx_last_y_luma[15];
extern cabac_ctx g_cu_ctx_last_y_chroma[15];
extern cabac_ctx g_cu_ctx_last_x_luma[15];
extern cabac_ctx g_cu_ctx_last_x_chroma[15];
extern cabac_ctx g_cu_one_model_luma[16];
extern cabac_ctx g_cu_one_model_chroma[8];
extern cabac_ctx g_cu_abs_model_luma[4];
extern cabac_ctx g_cu_abs_model_chroma[2];
extern cabac_ctx g_cu_pred_mode_model;
extern cabac_ctx g_cu_skip_flag_model[3];
extern cabac_ctx g_cu_merge_idx_ext_model;
extern cabac_ctx g_cu_merge_flag_ext_model;
extern cabac_ctx g_cu_mvd_model[2];
extern cabac_ctx g_cu_ref_pic_model[2];
extern cabac_ctx g_mvp_idx_model[2];
extern cabac_ctx g_cu_qt_root_cbf_model;
#define CNU 154

static const uint8_t INIT_SAO_MERGE_FLAG[3] = { 153, 153, 153 };
static const uint8_t INIT_SAO_TYPE_IDX[3] = { 160, 185, 200 };

static const uint8_t 
INIT_QT_ROOT_CBF[3][1] = 
{
  {  79, }, 
  {  79, }, 
  { CNU, }, 
};

static const uint8_t 
INIT_MVP_IDX[3][2] =  
{
  { 168,  CNU, }, 
  { 168,  CNU, }, 
  { CNU,  CNU, }, 
};

static const uint8_t 
INIT_REF_PIC[3][2] =  
{
  { 153,  153 }, 
  { 153,  153 }, 
  { CNU,  CNU }, 
};

static const uint8_t 
INIT_MVD[3][2] =  
{
  { 169,  198, }, 
  { 140,  198, }, 
  { CNU,  CNU, }, 
};

static const uint8_t
INIT_MERGE_FLAG_EXT[3][1] = 
{
  { 154, }, 
  { 110, }, 
  { CNU, }
};

static const uint8_t 
INIT_MERGE_IDX_EXT[3][1] =  
{
  { 137, }, 
  { 122, }, 
  { CNU, }
};

static const uint8_t 
INIT_SKIP_FLAG[3][3] =  
{
  { 197,  185,  201, }, 
  { 197,  185,  201, }, 
  { CNU,  CNU,  CNU, }
};

static const uint8_t 
INIT_PRED_MODE[3][1] = 
{
  { 134, }, 
  { 149, }, 
  { CNU, }
};


static const uint8_t 
INIT_PART_SIZE[3][4] =  
{
  { 154,  139,  CNU,  CNU, }, 
  { 154,  139,  CNU,  CNU, }, 
  { 184,  CNU,  CNU,  CNU, }, 
};

static const uint8_t  INIT_SPLIT_FLAG[3][3] =  
                       { { 107,  139,  126 },
                         { 107,  139,  126 },
                         { 139,  141,  157 } };

static const uint8_t INIT_INTRA_PRED_MODE[3] = { 183,154,184 };

static const uint8_t INIT_CHROMA_PRED_MODE[3][2] = { { 152,  139 }, { 152,  139 }, {  63,  139 } };


static const uint8_t INIT_TRANS_SUBDIV_FLAG[3][3] = 
{
  { 224,  167,  122 }, 
  { 124,  138,   94 }, 
  { 153,  138,  138 }
};

static const uint8_t INIT_QT_CBF[3][6] =  
{
  { 153,  111,  CNU,  149,   92,  167 },
  { 153,  111,  CNU,  149,  107,  167 },
  { 111,  141,  CNU,   94,  138,  182 }
};

static const uint8_t INIT_SIG_CG_FLAG[3][4] =  
 {  { 121,  140,  61,  154  },  { 121,  140,  61,  154 }, {  91,  171,  134,  141  } };

static const uint8_t INIT_SIG_FLAG[3][42] = 
  {
   {170,154,139,153,139,123,123, 63,124,166,183,140,136,153,154,166,183,140,136,153,154,166,183,140,136,153,154,170,153,138,138,122,121,122,121,167,151,183,140,151,183,140,},
   {155,154,139,153,139,123,123,63,153,166,183,140,136,153,154,166,183,140,136,153,154,166,183,140,136,153,154,170,153,123,123,107,121,107,121,167,151,183,140,151,183,140,},
   {111,111,125,110,110,94,124,108,124,107,125,141,179,153,125,107,125,141,179,153,125,107,125,141,179,153,125,140,139,182,182,152,136,152,136,153,136,139,111,136,139,111,}
  };

static const uint8_t INIT_LAST[3][30] =  
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79,
    108,  123,   93,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94,
    108,  123,  108,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79,
    108,  123,   63,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }
};

static const uint8_t INIT_ONE_FLAG[3][24] = 
{
  {154,196,167,167,154,152,167,182,182,134,149,136,153,121,136,122,169,208,166,167,154,152,167,182},
  {154,196,196,167,154,152,167,182,182,134,149,136,153,121,136,137,169,194,166,167,154,167,137,182},
  {140, 92,137,138,140,152,138,139,153, 74,149, 92,139,107,122,152,140,179,166,182,140,227,122,197}
};

static const uint8_t INIT_ABS_FLAG[3][6] =  
{
  { 107,167, 91,107,107,167}, 
  { 107,167, 91,122,107,167}, 
  { 138,153,136,167,152,152}, 
};


#endif
