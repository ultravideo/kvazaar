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

#include <getopt.h>

/**
 * \brief Allocate memory for config object
 * \return pointer to allocated memory
 */
config_t *config_alloc(void)
{
  config_t *cfg = (config_t *)malloc(sizeof(config_t));
  if (!cfg) {
    fprintf(stderr, "Failed to allocate a config object!\n");
    return cfg;
  }

  FILL(*cfg, 0);

  return cfg;
}

/**
 * \brief Initialize config structure
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_init(config_t *cfg)
{
  cfg->input           = NULL;
  cfg->output          = NULL;
  cfg->debug           = NULL;
  cfg->frames          = 0;
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
  cfg->seek            = 0;
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

/**
 * \brief Free memory allocated to the config
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_destroy(config_t *cfg)
{
  FREE_POINTER(cfg->input);
  FREE_POINTER(cfg->output);
  FREE_POINTER(cfg->cqmfile);
  FREE_POINTER(cfg->tiles_width_split);
  FREE_POINTER(cfg->tiles_height_split);
  FREE_POINTER(cfg->slice_addresses_in_ts);
  free(cfg);

  return 1;
}

/**
 * \brief Allocates memory space for a string, and copies it
 * \param char * string to copy
 * \return a pointer to the copied string on success, null on failure
 */
static char *copy_string(const char *string)
{
  // Allocate +1 for \0
  char *allocated_string = (char *)malloc(strlen(string) + 1);
  if (!allocated_string) {
    fprintf(stderr, "Failed to allocate a string!\n");
    return allocated_string;
  }

  // Copy the string to the new buffer
  memcpy(allocated_string, string, strlen(string) + 1);

  return allocated_string;
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
    if (MAX_TILES_PER_DIM <= *ntiles || 0 >= *ntiles) {
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

int config_parse(config_t *cfg, const char *name, const char *value)
{
  static const char * const me_names[]          = { "hexbs", "tz", NULL };

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
  if OPT("input")
    cfg->input = copy_string(value);
  else if OPT("output")
    cfg->output = copy_string(value);
  else if OPT("debug")
    cfg->debug = copy_string(value);
  else if OPT("width")
    cfg->width = atoi(value);
  else if OPT("height")
    cfg->height = atoi(value);
  else if OPT("input-res")
    return sscanf(value, "%dx%d", &cfg->width, &cfg->height) == 2;
  else if OPT("input-fps")
    cfg->framerate = atof(value);
  else if OPT("frames")
    cfg->frames = atoi(value);
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
  else if OPT("me")
    return parse_enum(value, me_names, &cfg->ime_algorithm);
  else if OPT("subme")
    cfg->fme_level = atoi(value);
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
    cfg->cqmfile = copy_string(value);
  else if OPT("seek")
    cfg->seek = atoi(value);
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
 * \brief Read configuration options from argv to the config struct
 * \param cfg config object
 * \param argc argument count
 * \param argv argument list
 * \return 1 on success, 0 on failure
 */
int config_read(config_t *cfg,int argc, char *argv[])
{
  static char short_options[] = "i:o:d:w:h:n:q:p:r:";
  static struct option long_options[] =
  {
    { "input",              required_argument, NULL, 'i' },
    { "output",             required_argument, NULL, 'o' },
    { "debug",              required_argument, NULL, 'd' },
    { "width",              required_argument, NULL, 'w' },
    { "height",             required_argument, NULL, 'h' }, // deprecated
    { "frames",             required_argument, NULL, 'n' }, // deprecated
    { "qp",                 required_argument, NULL, 'q' },
    { "period",             required_argument, NULL, 'p' },
    { "ref",                required_argument, NULL, 'r' },
    { "vps-period",         required_argument, NULL, 0 },
    { "input-res",          required_argument, NULL, 0 },
    { "input-fps",          required_argument, NULL, 0 },
    { "no-deblock",               no_argument, NULL, 0 },
    { "deblock",            required_argument, NULL, 0 },
    { "no-sao",                   no_argument, NULL, 0 },
    { "no-rdoq",                  no_argument, NULL, 0 },
    { "no-signhide",              no_argument, NULL, 0 },
    { "rd",                 required_argument, NULL, 0 },
    { "full-intra-search",        no_argument, NULL, 0 },
    { "no-transform-skip",        no_argument, NULL, 0 },
    { "tr-depth-intra",     required_argument, NULL, 0 },
    { "me",                 required_argument, NULL, 0 },
    { "subme",              required_argument, NULL, 0 },
    { "sar",                required_argument, NULL, 0 },
    { "overscan",           required_argument, NULL, 0 },
    { "videoformat",        required_argument, NULL, 0 },
    { "range",              required_argument, NULL, 0 },
    { "colorprim",          required_argument, NULL, 0 },
    { "transfer",           required_argument, NULL, 0 },
    { "colormatrix",        required_argument, NULL, 0 },
    { "chromaloc",          required_argument, NULL, 0 },
    { "aud",                      no_argument, NULL, 0 },
    { "cqmfile",            required_argument, NULL, 0 },
    { "seek",               required_argument, NULL, 0 },
    { "tiles-width-split",  required_argument, NULL, 0 },
    { "tiles-height-split", required_argument, NULL, 0 },
    { "wpp",                      no_argument, NULL, 0 },
    { "owf",                required_argument, NULL, 0 },
    { "slice-addresses",    required_argument, NULL, 0 },
    { "threads",            required_argument, NULL, 0 },
    { "cpuid",              required_argument, NULL, 0 },
    { "pu-depth-inter",     required_argument, NULL, 0 },
    { "pu-depth-intra",     required_argument, NULL, 0 },
    { "no-info",                  no_argument, NULL, 0 },
    { "gop",                required_argument, NULL, 0 },
    { "bipred",                   no_argument, NULL, 0 },
    { "bitrate",            required_argument, NULL, 0 },
    {0, 0, 0, 0}
  };

  // Parse command line options
  for (optind = 0;;) {
    int long_options_index = -1;

    int c = getopt_long(argc, argv, short_options, long_options, &long_options_index);
    if (c == -1)
      break;

    if (long_options_index < 0) {
      int i;
      for (i = 0; long_options[i].name; i++)
        if (long_options[i].val == c) {
            long_options_index = i;
            break;
        }
      if (long_options_index < 0) {
        // getopt_long already printed an error message
        return 0;
      }
    }

    if (!config_parse(cfg, long_options[long_options_index].name, optarg)) {
      const char *name = long_options_index > 0 ? long_options[long_options_index].name : argv[optind-2];
      fprintf(stderr, "invalid argument: %s=%s\n", name, optarg );
      return 0;
    }
  }

  // Check that the required files were defined
  if(cfg->input == NULL || cfg->output == NULL) return 0;

  // Add dimensions to the reconstructions file name.
  if (cfg->debug != NULL) {
    char dim_str[50]; // log10(2^64) < 20, so this should suffice. I hate C.
    size_t left_len, right_len;
    sprintf(dim_str, "_%dx%d.yuv", cfg->width, cfg->height);
    left_len = strlen(cfg->debug);
    right_len = strlen(dim_str);
    cfg->debug = realloc(cfg->debug, left_len + right_len + 1);
    if (!cfg->debug) {
      fprintf(stderr, "realloc failed!\n");
      return 0;
    }
    strcpy(cfg->debug + left_len, dim_str);
  }

  return 1;
}

/**
 * \brief Check that configuration is sensible.
 *
 * \param cfg   config to check
 * \return      1 if the config is ok, otherwise 1
 */
int config_validate(config_t const *const cfg)
{
  int error = 0;

  if (cfg->width <= 0) {
    fprintf(stderr, "Input error: width is not positive\n");
    error = 1;
  }

  if (cfg->framerate <= 0.0) {
    fprintf(stderr, "Input error: --input-fps must be positive\n");
    error = 1;
  }

  if (cfg->height <= 0) {
    fprintf(stderr, "Input error: height is not positive\n");
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

  if (cfg->rdo < 0 || cfg->rdo > 2) {
    fprintf(stderr, "Input error: --rd parameter out of range [0..2]\n");
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
