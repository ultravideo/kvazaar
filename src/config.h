#ifndef CONFIG_H_
#define CONFIG_H_
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
 * \brief Handles parsing and storing of configuration of the encoder.
 */

#include "kvazaar.h"
#include "global.h"

/* Function definitions */
config_t *config_alloc(void);
int config_init(config_t *cfg);
int config_destroy(config_t *cfg);
int config_read(config_t *cfg,int argc, char *argv[]);
int config_parse(config_t *cfg, const char *name, const char *value);
int config_validate(config_t const *cfg);
int config_set_owf_auto(config_t *cfg);

#endif
