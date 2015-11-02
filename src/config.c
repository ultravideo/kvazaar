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

/*
 * \file
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

kvz_config *kvz_config_alloc(void)
{
  kvz_config *cfg = (kvz_config *)malloc(sizeof(kvz_config));
  if (!cfg) {
    fprintf(stderr, "Failed to allocate a config object!\n");
    return cfg;
  }

  FILL(*cfg, 0);

  return cfg;
}

int kvz_config_init(kvz_config *cfg)
{
  cfg->width           = 0;
  cfg->height          = 0;
  cfg->framerate       = 25;
  cfg->qp              = 32;
  cfg->intra_period    = 0;
  cfg->vps_period      = 0;
  cfg->deblock_enable  = 1;
  cfg->deblock_beta    = 0;
  cfg->deblock_tc      = 0;
  cfg->sao_enable      = 1;
  cfg->rdoq_enable     = 1;
  cfg->signhide_enable = true;
  cfg->rdo             = 1;
  cfg->full_intra_search = 0;
  cfg->trskip_enable   = 1;
  cfg->tr_depth_intra  = 0;
  cfg->ime_algorithm   = 0; /* hexbs */
  cfg->fme_level       = 1;
  cfg->source_scan_type = 0; /* progressive */
  cfg->vui.sar_width   = 0;
  cfg->vui.sar_height  = 0;
  cfg->vui.overscan    = 0; /* undef */
  cfg->vui.videoformat = 5; /* undef */
  cfg->vui.fullrange   = 0; /* limited range */
  cfg->vui.colorprim   = 2; /* undef */
  cfg->vui.transfer    = 2; /* undef */
  cfg->vui.colormatrix = 2; /* undef */
  cfg->vui.chroma_loc  = 0; /* left center */
  cfg->aud_enable      = 0;
  cfg->cqmfile         = NULL;
  cfg->ref_frames      = DEFAULT_REF_PIC_COUNT;
  cfg->gop_len         = 0;
  cfg->bipred          = 0;
  cfg->target_bitrate  = 0;

  cfg->tiles_width_count         = 0;
  cfg->tiles_height_count         = 0;
  cfg->tiles_width_split          = NULL;
  cfg->tiles_height_split          = NULL;
  
  cfg->wpp = 0;
  cfg->owf = -1;
  cfg->slice_count = 1;
  cfg->slice_addresses_in_ts = MALLOC(int32_t, 1);
  cfg->slice_addresses_in_ts[0] = 0;
  
  cfg->threads = 0;
  cfg->cpuid = 1;

  // Defaults for what sizes of PUs are tried.
  cfg->pu_depth_inter.min = 0; // 0-3
  cfg->pu_depth_inter.max = 3; // 0-3
  cfg->pu_depth_intra.min = 1; // 0-4
  cfg->pu_depth_intra.max = 4; // 0-4

  cfg->add_encoder_info = true;

  return 1;
}

int kvz_config_destroy(kvz_config *cfg)
{
  if (cfg) {
    FREE_POINTER(cfg->cqmfile);
    FREE_POINTER(cfg->tiles_width_split);
    FREE_POINTER(cfg->tiles_height_split);
    FREE_POINTER(cfg->slice_addresses_in_ts);
  }
  free(cfg);

  return 1;
}

static int atobool(const char *str)
{
  if (!strcmp(str, "1")    ||
      !strcmp(str, "true") ||
      !strcmp(str, "yes"))
    return 1;
  if (!strcmp(str, "0")     ||
      !strcmp(str, "false") ||
      !strcmp(str, "no"))
    return 0;
  return 0;
}

static int parse_enum(const char *arg, const char * const *names, int8_t *dst)
{
  int8_t i;
  for (i = 0; names[i]; i++) {
    if (!strcmp(arg, names[i])) {
      *dst = i;
      return 1;
    }
  }

  return 0;
}

static int parse_tiles_specification(const char* const arg, int32_t * const ntiles, int32_t** const array) {
  const char* current_arg = NULL;
  int32_t current_value;
  int32_t values[MAX_TILES_PER_DIM];
  
  int i;
  
  //Free pointer in any case
  if (*array) {
    FREE_POINTER(*array);
  }
  
  //If the arg starts with u, we want an uniform split
  if (arg[0]=='u') {
    *ntiles = atoi(arg+1)-1;
    if (MAX_TILES_PER_DIM <= *ntiles || 0 > *ntiles) {
      fprintf(stderr, "Invalid number of tiles (0 < %d <= %d = MAX_TILES_PER_DIM)!\n", *ntiles + 1, MAX_TILES_PER_DIM);
      return 0;
    }
    //Done with parsing
    return 1;
  }
  
  //We have a comma-separated list of int for the split...
  current_arg = arg;
  *ntiles = 0;
  do {
    int ret = sscanf(current_arg, "%d", &current_value);
    if (ret != 1) {
      fprintf(stderr, "Could not parse integer \"%s\"!\n", current_arg);
      return 0;
    }
    current_arg = strchr(current_arg, ',');
    //Skip the , if we found one
    if (current_arg) ++current_arg;
    values[*ntiles] = current_value;
    ++(*ntiles);
    if (MAX_TILES_PER_DIM <= *ntiles) break;
  } while (current_arg);
  
  if (MAX_TILES_PER_DIM <= *ntiles || 0 >= *ntiles) {
    fprintf(stderr, "Invalid number of tiles (0 < %d <= %d = MAX_TILES_PER_DIM)!\n", *ntiles + 1, MAX_TILES_PER_DIM);
    return 0;
  }
  
  *array = MALLOC(int32_t, *ntiles);
  if (!*array) {
    fprintf(stderr, "Could not allocate array for tiles\n");
    return 0;
  }
  
  //TODO: memcpy?
  for (i = 0; i < *ntiles; ++i) {
    (*array)[i] = values[i];
  }
  
  return 1;
}

static int parse_slice_specification(const char* const arg, int32_t * const nslices, int32_t** const array) {
  const char* current_arg = NULL;
  int32_t current_value;
  int32_t values[MAX_SLICES];
  
  int i;
  
  //Free pointer in any case
  if (*array) {
    FREE_POINTER(*array);
  }
  
  //If the arg starts with u, we want an uniform split
  if (arg[0]=='u') {
    *nslices = atoi(arg+1);
    if (MAX_SLICES <= *nslices || 0 >= *nslices) {
      fprintf(stderr, "Invalid number of tiles (0 < %d <= %d = MAX_SLICES)!\n", *nslices + 1, MAX_SLICES);
      return 0;
    }
    //Done with parsing
    return 1;
  }
  
  //We have a comma-separated list of int for the split...
  current_arg = arg;
  //We always have a slice starting at 0
  values[0] = 0;
  *nslices = 1;
  do {
    int ret = sscanf(current_arg, "%d", &current_value);
    if (ret != 1) {
      fprintf(stderr, "Could not parse integer \"%s\"!\n", current_arg);
      return 0;
    }
    current_arg = strchr(current_arg, ',');
    //Skip the , if we found one
    if (current_arg) ++current_arg;
    values[*nslices] = current_value;
    ++(*nslices);
    if (MAX_SLICES <= *nslices) break;
  } while (current_arg);
  
  if (MAX_SLICES <= *nslices || 0 >= *nslices) {
    fprintf(stderr, "Invalid number of slices (0 < %d <= %d = MAX_SLICES)!\n", *nslices, MAX_SLICES);
    return 0;
  }
  
  *array = MALLOC(int32_t, *nslices);
  if (!*array) {
    fprintf(stderr, "Could not allocate array for slices\n");
    return 0;
  }
  
  //TODO: memcpy?
  for (i = 0; i < *nslices; ++i) {
    (*array)[i] = values[i];
  }
  
  return 1;
}

int kvz_config_parse(kvz_config *cfg, const char *name, const char *value)
{
  static const char * const me_names[]          = { "hexbs", "tz", NULL };
  static const char * const source_scan_type_names[] = { "progressive", "tff", "bff", NULL };

  static const char * const overscan_names[]    = { "undef", "show", "crop", NULL };
  static const char * const videoformat_names[] = { "component", "pal", "ntsc", "secam", "mac", "undef", NULL };
  static const char * const range_names[]       = { "tv", "pc", NULL };
  static const char * const colorprim_names[]   = { "", "bt709", "undef", "", "bt470m", "bt470bg", "smpte170m",
                                                    "smpte240m", "film", "bt2020", NULL };
  static const char * const transfer_names[]    = { "", "bt709", "undef", "", "bt470m", "bt470bg", "smpte170m",
                                                    "smpte240m", "linear", "log100", "log316", "iec61966-2-4",
                                                    "bt1361e", "iec61966-2-1", "bt2020-10", "bt2020-12", NULL };
  static const char * const colormatrix_names[] = { "GBR", "bt709", "undef", "", "fcc", "bt470bg", "smpte170m",
                                                    "smpte240m", "YCgCo", "bt2020nc", "bt2020c", NULL };

  if (!name)
    return 0;

  if (!value)
    value = "true";

  // Treat "--no-param" as --param 0
  if ((!strncmp(name, "no-", 3))) {
    name += 3;
    value = atobool(value) ? "false" : "true";
  }

#define OPT(STR) (!strcmp(name, STR))
  if OPT("width")
    cfg->width = atoi(value);
  else if OPT("height")
    cfg->height = atoi(value);
  else if OPT("input-res")
    if (!strcmp(value, "auto")) {
      return 1;
    } else {
      return (sscanf(value, "%dx%d", &cfg->width, &cfg->height) == 2);
    }
  else if OPT("input-fps")
    cfg->framerate = atof(value);
  else if OPT("qp")
    cfg->qp = atoi(value);
  else if OPT("period")
    cfg->intra_period = atoi(value);
  else if OPT("vps-period")
    cfg->vps_period = atoi(value);
  else if OPT("ref")
    cfg->ref_frames = atoi(value);
  else if OPT("deblock") {
    int beta, tc;
    if (2 == sscanf(value, "%d:%d", &beta, &tc)) {
      cfg->deblock_enable = 1;
      cfg->deblock_beta   = beta;
      cfg->deblock_tc     = tc;
    } else if (sscanf(value, "%d", &beta)) {
      cfg->deblock_enable = 1;
      cfg->deblock_beta   = beta;
      cfg->deblock_tc     = cfg->deblock_beta;
    } else {
      cfg->deblock_enable = atobool(value);
    }
  }
  else if OPT("sao")
    cfg->sao_enable = atobool(value);
  else if OPT("rdoq")
    cfg->rdoq_enable = atobool(value);
  else if OPT("signhide")
    cfg->signhide_enable = (bool)atobool(value);
  else if OPT("rd")
    cfg->rdo = atoi(value);
  else if OPT("full-intra-search")
    cfg->full_intra_search = atobool(value);
  else if OPT("transform-skip")
    cfg->trskip_enable = atobool(value);
  else if OPT("tr-depth-intra")
    cfg->tr_depth_intra = atoi(value);
  else if OPT("me") {
    int8_t ime_algorithm = 0;
    if (!parse_enum(value, me_names, &ime_algorithm)) return 0;
    cfg->ime_algorithm = ime_algorithm;
  }
  else if OPT("subme")
    cfg->fme_level = atoi(value);
  else if OPT("source-scan-type")
    return parse_enum(value, source_scan_type_names, &cfg->source_scan_type);
  else if OPT("sar")
    return sscanf(value, "%d:%d", &cfg->vui.sar_width, &cfg->vui.sar_height) == 2;
  else if OPT("overscan")
    return parse_enum(value, overscan_names, &cfg->vui.overscan);
  else if OPT("videoformat")
    return parse_enum(value, videoformat_names, &cfg->vui.videoformat);
  else if OPT("range")
    return parse_enum(value, range_names, &cfg->vui.fullrange);
  else if OPT("colorprim")
    return parse_enum(value, colorprim_names, &cfg->vui.colorprim);
  else if OPT("transfer")
    return parse_enum(value, transfer_names, &cfg->vui.transfer);
  else if OPT("colormatrix")
    return parse_enum(value, colormatrix_names, &cfg->vui.colormatrix);
  else if OPT("chromaloc")
    cfg->vui.chroma_loc = atoi(value);
  else if OPT("aud")
    cfg->aud_enable = atobool(value);
  else if OPT("cqmfile")
    cfg->cqmfile = strdup(value);
  else if OPT("tiles-width-split")
    return parse_tiles_specification(value, &cfg->tiles_width_count, &cfg->tiles_width_split);
  else if OPT("tiles-height-split")
    return parse_tiles_specification(value, &cfg->tiles_height_count, &cfg->tiles_height_split);
  else if OPT("wpp")
    cfg->wpp = atobool(value);
  else if OPT("owf") {
    cfg->owf = atoi(value);
    if (cfg->owf == 0 && !strcmp(value, "auto")) {
      // -1 means automatic selection
      cfg->owf = -1;
    }
  }
  else if OPT("slice-addresses")
    return parse_slice_specification(value, &cfg->slice_count, &cfg->slice_addresses_in_ts);
  else if OPT("threads")
    cfg->threads = atoi(value);
  else if OPT("cpuid")
    cfg->cpuid = atoi(value);
  else if OPT("pu-depth-inter")
    return sscanf(value, "%d-%d", &cfg->pu_depth_inter.min, &cfg->pu_depth_inter.max) == 2;
  else if OPT("pu-depth-intra")
    return sscanf(value, "%d-%d", &cfg->pu_depth_intra.min, &cfg->pu_depth_intra.max) == 2;
  else if OPT("info")
    cfg->add_encoder_info = atobool(value);
  else if OPT("gop") {
    // TODO: Defining the whole GOP structure via parameters
    if(atoi(value) == 8) {
      cfg->gop_lowdelay = 0;
      // GOP
      cfg->gop_len = 8;
      cfg->gop[0].poc_offset = 8; cfg->gop[0].qp_offset = 1; cfg->gop[0].layer = 1; cfg->gop[0].qp_factor = 0.442;  cfg->gop[0].is_ref = 1;
      cfg->gop[0].ref_pos_count = 0;
      cfg->gop[0].ref_neg_count = 3; cfg->gop[0].ref_neg[0] = 8; cfg->gop[0].ref_neg[1] = 12; cfg->gop[0].ref_neg[2] = 16;

      cfg->gop[1].poc_offset = 4; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.3536; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 2; cfg->gop[1].ref_neg[0] = 4; cfg->gop[1].ref_neg[1] = 8;
      cfg->gop[1].ref_pos_count = 1; cfg->gop[1].ref_pos[0] = 4;

      cfg->gop[2].poc_offset = 2; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.3536; cfg->gop[2].is_ref = 1;
      cfg->gop[2].ref_neg_count = 2; cfg->gop[2].ref_neg[0] = 2; cfg->gop[2].ref_neg[1] = 6;
      cfg->gop[2].ref_pos_count = 2; cfg->gop[2].ref_pos[0] = 2; cfg->gop[2].ref_pos[1] = 6;

      cfg->gop[3].poc_offset = 1; cfg->gop[3].qp_offset = 4; cfg->gop[3].layer = 4; cfg->gop[3].qp_factor = 0.68;   cfg->gop[3].is_ref = 0;
      cfg->gop[3].ref_neg_count = 1; cfg->gop[3].ref_neg[0] = 1;
      cfg->gop[3].ref_pos_count = 3; cfg->gop[3].ref_pos[0] = 1; cfg->gop[3].ref_pos[1] = 3; cfg->gop[3].ref_pos[2] = 7;

      cfg->gop[4].poc_offset = 3; cfg->gop[4].qp_offset = 4; cfg->gop[4].layer = 4; cfg->gop[4].qp_factor = 0.68;   cfg->gop[4].is_ref = 0;
      cfg->gop[4].ref_neg_count = 2; cfg->gop[4].ref_neg[0] = 1; cfg->gop[4].ref_neg[1] = 3;
      cfg->gop[4].ref_pos_count = 2; cfg->gop[4].ref_pos[0] = 1; cfg->gop[4].ref_pos[1] = 5;

      cfg->gop[5].poc_offset = 6; cfg->gop[5].qp_offset = 3; cfg->gop[5].layer = 3; cfg->gop[5].qp_factor = 0.3536; cfg->gop[5].is_ref = 1;
      cfg->gop[5].ref_neg_count = 2; cfg->gop[5].ref_neg[0] = 2; cfg->gop[5].ref_neg[1] = 6;
      cfg->gop[5].ref_pos_count = 1; cfg->gop[5].ref_pos[0] = 2;

      cfg->gop[6].poc_offset = 5; cfg->gop[6].qp_offset = 4; cfg->gop[6].layer = 4; cfg->gop[6].qp_factor = 0.68;   cfg->gop[6].is_ref = 0;
      cfg->gop[6].ref_neg_count = 2;  cfg->gop[6].ref_neg[0] = 1; cfg->gop[6].ref_neg[1] = 5;
      cfg->gop[6].ref_pos_count = 2; cfg->gop[6].ref_pos[0] = 1; cfg->gop[6].ref_pos[1] = 3;

      cfg->gop[7].poc_offset = 7; cfg->gop[7].qp_offset = 4; cfg->gop[7].layer = 4; cfg->gop[7].qp_factor = 0.68;   cfg->gop[7].is_ref = 0;
      cfg->gop[7].ref_neg_count = 3; cfg->gop[7].ref_neg[0] = 1; cfg->gop[7].ref_neg[1] = 3; cfg->gop[7].ref_neg[2] = 7;
      cfg->gop[7].ref_pos_count = 1; cfg->gop[7].ref_pos[0] = 1;
    } else if(!strcmp(value, "lb")) { // Low-delay B
      cfg->gop_lowdelay = 1;
      cfg->gop_len = 4;
      cfg->gop[0].poc_offset = 1; cfg->gop[0].qp_offset = 3; cfg->gop[0].layer = 3; cfg->gop[0].qp_factor = 0.4624;  cfg->gop[0].is_ref = 1;      
      cfg->gop[0].ref_neg_count = 4; cfg->gop[0].ref_neg[0] = 1; cfg->gop[0].ref_neg[1] = 5; cfg->gop[0].ref_neg[2] = 9; cfg->gop[0].ref_neg[3] = 13;
      cfg->gop[0].ref_pos_count = 0;

      cfg->gop[1].poc_offset = 2; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.4624; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 4; cfg->gop[1].ref_neg[0] = 1; cfg->gop[1].ref_neg[1] = 2; cfg->gop[1].ref_neg[2] = 6; cfg->gop[1].ref_neg[3] = 10;
      cfg->gop[1].ref_pos_count = 0;

      cfg->gop[2].poc_offset = 3; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.4624; cfg->gop[2].is_ref = 1;
      cfg->gop[2].ref_neg_count = 4; cfg->gop[2].ref_neg[0] = 1; cfg->gop[2].ref_neg[1] = 3; cfg->gop[2].ref_neg[2] = 7; cfg->gop[2].ref_neg[3] = 11;
      cfg->gop[2].ref_pos_count = 0;

      cfg->gop[3].poc_offset = 4; cfg->gop[3].qp_offset = 1; cfg->gop[3].layer = 1; cfg->gop[3].qp_factor = 0.578;   cfg->gop[3].is_ref = 1;
      cfg->gop[3].ref_neg_count = 4; cfg->gop[3].ref_neg[0] = 1; cfg->gop[3].ref_neg[1] = 4; cfg->gop[3].ref_neg[2] = 8; cfg->gop[3].ref_neg[3] = 12;
      cfg->gop[3].ref_pos_count = 0;
    } else if (!strcmp(value, "g8r1t1")) { // Low-delay B with 1 ref frame
      cfg->gop_lowdelay = 1;
      cfg->gop_len = 8;
      cfg->gop[0].poc_offset = 1; cfg->gop[0].qp_offset = 3; cfg->gop[0].layer = 3; cfg->gop[0].qp_factor = 0.4624;  cfg->gop[0].is_ref = 1;
      cfg->gop[0].ref_pos_count = 0;
      cfg->gop[0].ref_neg_count = 1; cfg->gop[0].ref_neg[0] = 1;

      cfg->gop[1].poc_offset = 2; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.4624; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 1; cfg->gop[1].ref_neg[0] = 1;
      cfg->gop[1].ref_pos_count = 0;

      cfg->gop[2].poc_offset = 3; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.4624; cfg->gop[2].is_ref = 1;
      cfg->gop[2].ref_neg_count = 1; cfg->gop[2].ref_neg[0] = 1;
      cfg->gop[2].ref_pos_count = 0;

      cfg->gop[3].poc_offset = 4; cfg->gop[3].qp_offset = 2; cfg->gop[3].layer = 2; cfg->gop[3].qp_factor = 0.4624;   cfg->gop[3].is_ref = 1;
      cfg->gop[3].ref_neg_count = 1; cfg->gop[3].ref_neg[0] = 1;
      cfg->gop[3].ref_pos_count = 0;

      cfg->gop[4].poc_offset = 5; cfg->gop[4].qp_offset = 3; cfg->gop[4].layer = 3; cfg->gop[4].qp_factor = 0.4624;   cfg->gop[4].is_ref = 1;
      cfg->gop[4].ref_neg_count = 1; cfg->gop[4].ref_neg[0] = 1;
      cfg->gop[4].ref_pos_count = 0;

      cfg->gop[5].poc_offset = 6; cfg->gop[5].qp_offset = 2; cfg->gop[5].layer = 2; cfg->gop[5].qp_factor = 0.4624; cfg->gop[5].is_ref = 1;
      cfg->gop[5].ref_neg_count = 1; cfg->gop[5].ref_neg[0] = 1;
      cfg->gop[5].ref_pos_count = 0;

      cfg->gop[6].poc_offset = 7; cfg->gop[6].qp_offset = 3; cfg->gop[6].layer = 3; cfg->gop[6].qp_factor = 0.4624;   cfg->gop[6].is_ref = 1;
      cfg->gop[6].ref_neg_count = 1;  cfg->gop[6].ref_neg[0] = 1;
      cfg->gop[6].ref_pos_count = 0;

      cfg->gop[7].poc_offset = 8; cfg->gop[7].qp_offset = 1; cfg->gop[7].layer = 1; cfg->gop[7].qp_factor = 0.578;   cfg->gop[7].is_ref = 1;
      cfg->gop[7].ref_neg_count = 1; cfg->gop[7].ref_neg[0] = 1;
      cfg->gop[7].ref_pos_count = 0;
    } else if (!strcmp(value, "g8r1t2")) { // Low-delay B with 1 ref frame
      cfg->gop_lowdelay = 1;
      cfg->gop_len = 8;
      cfg->gop[0].poc_offset = 1; cfg->gop[0].qp_offset = 3; cfg->gop[0].layer = 3; cfg->gop[0].qp_factor = 0.4624;  cfg->gop[0].is_ref = 0;
      cfg->gop[0].ref_pos_count = 0;
      cfg->gop[0].ref_neg_count = 1; cfg->gop[0].ref_neg[0] = 1;

      cfg->gop[1].poc_offset = 2; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.4624; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 1; cfg->gop[1].ref_neg[0] = 2;
      cfg->gop[1].ref_pos_count = 0;

      cfg->gop[2].poc_offset = 3; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.4624; cfg->gop[2].is_ref = 0;
      cfg->gop[2].ref_neg_count = 1; cfg->gop[2].ref_neg[0] = 1;
      cfg->gop[2].ref_pos_count = 0;

      cfg->gop[3].poc_offset = 4; cfg->gop[3].qp_offset = 2; cfg->gop[3].layer = 2; cfg->gop[3].qp_factor = 0.4624;   cfg->gop[3].is_ref = 1;
      cfg->gop[3].ref_neg_count = 1; cfg->gop[3].ref_neg[0] = 2;
      cfg->gop[3].ref_pos_count = 0;

      cfg->gop[4].poc_offset = 5; cfg->gop[4].qp_offset = 3; cfg->gop[4].layer = 3; cfg->gop[4].qp_factor = 0.4624;   cfg->gop[4].is_ref = 0;
      cfg->gop[4].ref_neg_count = 1; cfg->gop[4].ref_neg[0] = 1;
      cfg->gop[4].ref_pos_count = 0;

      cfg->gop[5].poc_offset = 6; cfg->gop[5].qp_offset = 2; cfg->gop[5].layer = 2; cfg->gop[5].qp_factor = 0.4624; cfg->gop[5].is_ref = 1;
      cfg->gop[5].ref_neg_count = 1; cfg->gop[5].ref_neg[0] = 2;
      cfg->gop[5].ref_pos_count = 0;

      cfg->gop[6].poc_offset = 7; cfg->gop[6].qp_offset = 3; cfg->gop[6].layer = 3; cfg->gop[6].qp_factor = 0.4624;   cfg->gop[6].is_ref = 0;
      cfg->gop[6].ref_neg_count = 1;  cfg->gop[6].ref_neg[0] = 1;
      cfg->gop[6].ref_pos_count = 0;

      cfg->gop[7].poc_offset = 8; cfg->gop[7].qp_offset = 1; cfg->gop[7].layer = 1; cfg->gop[7].qp_factor = 0.578;   cfg->gop[7].is_ref = 1;
      cfg->gop[7].ref_neg_count = 1; cfg->gop[7].ref_neg[0] = 2;
      cfg->gop[7].ref_pos_count = 0;
    }  else if (!strcmp(value, "g8r2t1")) {
      cfg->gop_lowdelay = 1;
      cfg->gop_len = 8;
      cfg->gop[0].poc_offset = 1; cfg->gop[0].qp_offset = 3; cfg->gop[0].layer = 3; cfg->gop[0].qp_factor = 0.4624;  cfg->gop[0].is_ref = 1;
      cfg->gop[0].ref_neg_count = 2; cfg->gop[0].ref_neg[0] = 1; cfg->gop[0].ref_neg[1] = 9;
      cfg->gop[0].ref_pos_count = 0;

      cfg->gop[1].poc_offset = 2; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.4624; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 2; cfg->gop[1].ref_neg[0] = 1; cfg->gop[1].ref_neg[1] = 2;
      cfg->gop[1].ref_pos_count = 0;

      cfg->gop[2].poc_offset = 3; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.4624; cfg->gop[2].is_ref = 1;
      cfg->gop[2].ref_neg_count = 2; cfg->gop[2].ref_neg[0] = 1; cfg->gop[2].ref_neg[1] = 3;
      cfg->gop[2].ref_pos_count = 0;

      cfg->gop[3].poc_offset = 4; cfg->gop[3].qp_offset = 2; cfg->gop[3].layer = 2; cfg->gop[3].qp_factor = 0.4624;   cfg->gop[3].is_ref = 1;
      cfg->gop[3].ref_neg_count = 2; cfg->gop[3].ref_neg[0] = 1; cfg->gop[3].ref_neg[1] = 4;
      cfg->gop[3].ref_pos_count = 0;

      cfg->gop[4].poc_offset = 5; cfg->gop[4].qp_offset = 3; cfg->gop[4].layer = 3; cfg->gop[4].qp_factor = 0.4624;   cfg->gop[4].is_ref = 1;
      cfg->gop[4].ref_neg_count = 2; cfg->gop[4].ref_neg[0] = 1; cfg->gop[4].ref_neg[1] = 5;
      cfg->gop[4].ref_pos_count = 0;

      cfg->gop[5].poc_offset = 6; cfg->gop[5].qp_offset = 2; cfg->gop[5].layer = 2; cfg->gop[5].qp_factor = 0.4624; cfg->gop[5].is_ref = 1;
      cfg->gop[5].ref_neg_count = 2; cfg->gop[5].ref_neg[0] = 1; cfg->gop[5].ref_neg[1] = 6;
      cfg->gop[5].ref_pos_count = 0;

      cfg->gop[6].poc_offset = 7; cfg->gop[6].qp_offset = 3; cfg->gop[6].layer = 3; cfg->gop[6].qp_factor = 0.4624;   cfg->gop[6].is_ref = 1;
      cfg->gop[6].ref_neg_count = 2; cfg->gop[6].ref_neg[0] = 1; cfg->gop[5].ref_neg[1] = 7;
      cfg->gop[6].ref_pos_count = 0;

      cfg->gop[7].poc_offset = 8; cfg->gop[7].qp_offset = 1; cfg->gop[7].layer = 1; cfg->gop[7].qp_factor = 0.578;   cfg->gop[7].is_ref = 1;
      cfg->gop[7].ref_neg_count = 2; cfg->gop[7].ref_neg[0] = 1; cfg->gop[7].ref_neg[1] = 8;
      cfg->gop[7].ref_pos_count = 0;
    } else if (!strcmp(value, "g8r2t2")) {
      cfg->gop_lowdelay = 1;
      cfg->gop_len = 8;
      cfg->gop[0].poc_offset = 1; cfg->gop[0].qp_offset = 3; cfg->gop[0].layer = 3; cfg->gop[0].qp_factor = 0.4624;  cfg->gop[0].is_ref = 0;
      cfg->gop[0].ref_neg_count = 2; cfg->gop[0].ref_neg[0] = 1; cfg->gop[0].ref_neg[1] = 9;
      cfg->gop[0].ref_pos_count = 0;

      cfg->gop[1].poc_offset = 2; cfg->gop[1].qp_offset = 2; cfg->gop[1].layer = 2; cfg->gop[1].qp_factor = 0.4624; cfg->gop[1].is_ref = 1;
      cfg->gop[1].ref_neg_count = 2; cfg->gop[1].ref_neg[0] = 2; cfg->gop[1].ref_neg[1] = 10;
      cfg->gop[1].ref_pos_count = 0;

      cfg->gop[2].poc_offset = 3; cfg->gop[2].qp_offset = 3; cfg->gop[2].layer = 3; cfg->gop[2].qp_factor = 0.4624; cfg->gop[2].is_ref = 0;
      cfg->gop[2].ref_neg_count = 2; cfg->gop[2].ref_neg[0] = 1; cfg->gop[2].ref_neg[1] = 3;
      cfg->gop[2].ref_pos_count = 0;

      cfg->gop[3].poc_offset = 4; cfg->gop[3].qp_offset = 2; cfg->gop[3].layer = 2; cfg->gop[3].qp_factor = 0.4624;   cfg->gop[3].is_ref = 1;
      cfg->gop[3].ref_neg_count = 2; cfg->gop[3].ref_neg[0] = 2; cfg->gop[3].ref_neg[1] = 4;
      cfg->gop[3].ref_pos_count = 0;

      cfg->gop[4].poc_offset = 5; cfg->gop[4].qp_offset = 3; cfg->gop[4].layer = 3; cfg->gop[4].qp_factor = 0.4624;   cfg->gop[4].is_ref = 0;
      cfg->gop[4].ref_neg_count = 2; cfg->gop[4].ref_neg[0] = 1; cfg->gop[4].ref_neg[1] = 5;
      cfg->gop[4].ref_pos_count = 0;

      cfg->gop[5].poc_offset = 6; cfg->gop[5].qp_offset = 2; cfg->gop[5].layer = 2; cfg->gop[5].qp_factor = 0.4624; cfg->gop[5].is_ref = 1;
      cfg->gop[5].ref_neg_count = 2; cfg->gop[5].ref_neg[0] = 2; cfg->gop[5].ref_neg[1] = 6;
      cfg->gop[5].ref_pos_count = 0;

      cfg->gop[6].poc_offset = 7; cfg->gop[6].qp_offset = 3; cfg->gop[6].layer = 3; cfg->gop[6].qp_factor = 0.4624;   cfg->gop[6].is_ref = 0;
      cfg->gop[6].ref_neg_count = 2; cfg->gop[6].ref_neg[0] = 1; cfg->gop[5].ref_neg[1] = 7;
      cfg->gop[6].ref_pos_count = 0;

      cfg->gop[7].poc_offset = 8; cfg->gop[7].qp_offset = 1; cfg->gop[7].layer = 1; cfg->gop[7].qp_factor = 0.578;   cfg->gop[7].is_ref = 1;
      cfg->gop[7].ref_neg_count = 2; cfg->gop[7].ref_neg[0] = 2; cfg->gop[7].ref_neg[1] = 8;
      cfg->gop[7].ref_pos_count = 0;
    } else if (atoi(value)) {
      fprintf(stderr, "Input error: unsupported gop length, must be 0 or 8\n");
      return 0;
    }
  }
  else if OPT("bipred")
    cfg->bipred = atobool(value);
  else if OPT("bitrate")
    cfg->target_bitrate = atoi(value);
  else
    return 0;
#undef OPT

  return 1;
}

/**
 * \brief Check that configuration is sensible.
 *
 * \param cfg   config to check
 * \return      1 if the config is ok, otherwise 1
 */
int kvz_config_validate(const kvz_config *const cfg)
{
  int error = 0;

  if (cfg->width <= 0) {
    fprintf(stderr, "Input error: width must be positive\n");
    error = 1;
  }

  if (cfg->height <= 0) {
    fprintf(stderr, "Input error: height must be positive\n");
    error = 1;
  }

  if (cfg->width % 2 != 0) {
    fprintf(stderr, "Input error: width must be a multiple of two\n");
    error = 1;
  }

  if (cfg->height % 2 != 0) {
    fprintf(stderr, "Input error: height must be a multiple of two\n");
    error = 1;
  }

  if (cfg->framerate <= 0.0) {
    fprintf(stderr, "Input error: --input-fps must be positive\n");
    error = 1;
  }

  if (cfg->gop_len &&
      cfg->intra_period &&
      cfg->intra_period % cfg->gop_len != 0) {
    fprintf(stderr,
            "Input error: intra period (%d) not a multiple of gop length (%d)\n",
            cfg->intra_period,
            cfg->gop_len);
    error = 1;
  }

  if (cfg->ref_frames  < 1 || cfg->ref_frames >= MAX_REF_PIC_COUNT) {
    fprintf(stderr, "Input error: --ref out of range [1..%d]\n", MAX_REF_PIC_COUNT - 1);
    error = 1;
  }

  if (cfg->deblock_beta  < -6 || cfg->deblock_beta  > 6) {
    fprintf(stderr, "Input error: deblock beta parameter out of range [-6..6]\n");
    error = 1;
  }
  if (cfg->deblock_tc < -6 || cfg->deblock_tc > 6) {
    fprintf(stderr, "Input error: deblock tc parameter out of range [-6..6]\n");
    error = 1;
  }

  if (cfg->rdo < 0 || cfg->rdo > 3) {
    fprintf(stderr, "Input error: --rd parameter out of range [0..3]\n");
    error = 1;
  }

  if (cfg->tr_depth_intra < 0 || cfg->tr_depth_intra > 4) {
    // range is 0 .. CtbLog2SizeY - Log2MinTrafoSize
    fprintf(stderr, "Input error: --tr-depth-intra is out of range [0..4]\n");
    error = 1;
  }

  if (cfg->fme_level != 0 && cfg->fme_level != 1) {
    fprintf(stderr, "Input error: invalid --subme parameter (must be 0 or 1)\n");
    error = 1;
  }

  if (cfg->vui.chroma_loc < 0 || cfg->vui.chroma_loc > 5) {
    fprintf(stderr, "Input error: --chromaloc parameter out of range [0..5]\n");
    error = 1;
  }

  if (cfg->owf < -1) {
    fprintf(stderr, "Input error: --owf must be nonnegative or -1\n");
    error = 1;
  }

  if (cfg->target_bitrate < 0) {
      fprintf(stderr, "Input error: --bitrate must be nonnegative\n");
      error = 1;
  }

  if (!WITHIN(cfg->pu_depth_inter.min, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX) ||
      !WITHIN(cfg->pu_depth_inter.max, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX)) 
  {
    fprintf(stderr, "Input error: illegal value for --pu-depth-inter (%d-%d)\n",
            cfg->pu_depth_inter.min, cfg->pu_depth_inter.max);
    error = 1;
  } else if (cfg->pu_depth_inter.min > cfg->pu_depth_inter.max) {
    fprintf(stderr, "Input error: Inter PU depth min (%d) > max (%d)\n",
            cfg->pu_depth_inter.min, cfg->pu_depth_inter.max);
    error = 1;
  }

  if (!WITHIN(cfg->pu_depth_intra.min, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX) ||
      !WITHIN(cfg->pu_depth_intra.max, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX))
  {
    fprintf(stderr, "Input error: illegal value for --pu-depth-intra (%d-%d)\n",
      cfg->pu_depth_intra.min, cfg->pu_depth_intra.max);
    error = 1;
  } else if (cfg->pu_depth_intra.min > cfg->pu_depth_intra.max) {
    fprintf(stderr, "Input error: Intra PU depth min (%d) > max (%d)\n",
            cfg->pu_depth_intra.min, cfg->pu_depth_intra.max);
    error = 1;
  }

  // Tile separation should be at round position in terms of LCU, should be monotonic, and should not start by 0
  if (cfg->tiles_width_split) {
    int i;
    int32_t prev_tile_split = 0;
    for (i=0; i < cfg->tiles_width_count; ++i) {
      if (cfg->tiles_width_split[i] <= prev_tile_split) {
        fprintf(stderr, "Input error: tile separations in width should be strictly monotonic (%d <= %d)\n", cfg->tiles_width_split[i], prev_tile_split);
        error = 1;
        break;
      }
      if ((cfg->tiles_width_split[i] % LCU_WIDTH) != 0) {
        fprintf(stderr, "Input error: tile separation in width %d (at %d) is not at a multiple of LCU_WIDTH (%d)\n", i, cfg->tiles_width_split[i], LCU_WIDTH);
        error = 1;
        break;
      }
      prev_tile_split = cfg->tiles_width_split[i];
    }
    if (cfg->tiles_width_split[cfg->tiles_width_count-1] >= cfg->width) {
      fprintf(stderr, "Input error: last x tile separation in width (%d) should smaller than image width (%d)\n", cfg->tiles_width_split[cfg->tiles_width_count-1], cfg->width);
      error = 1;
    }
  }

  if (cfg->tiles_height_split) {
    int i;
    int32_t prev_tile_split = 0;
    for (i=0; i < cfg->tiles_height_count; ++i) {
      if (cfg->tiles_height_split[i] <= prev_tile_split) {
        fprintf(stderr, "Input error: tile separations in height should be strictly monotonic (%d <= %d)\n", cfg->tiles_height_split[i], prev_tile_split);
        error = 1;
        break;
      }
      if ((cfg->tiles_height_split[i] % LCU_WIDTH) != 0) {
        fprintf(stderr, "Input error: tile separation in height %d (at %d) is not at a multiple of LCU_WIDTH (%d)\n", i, cfg->tiles_height_split[i], LCU_WIDTH);
        error = 1;
        break;
      }
      prev_tile_split = cfg->tiles_height_split[i];
    }

    if (cfg->tiles_height_split[cfg->tiles_height_count-1] >= cfg->height) {
      fprintf(stderr, "Input error: last tile separation in height (%d) should smaller than image height (%d)\n", cfg->tiles_height_split[cfg->tiles_height_count-1], cfg->height);
      error = 1;
    }
  }

  return !error;
}
