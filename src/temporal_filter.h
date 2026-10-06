#ifndef TEMPORAL_FILTER_H_
#define TEMPORAL_FILTER_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Control
 * \file
 * \brief Motion compensated temporal pre-filter (MCTF) for key pictures.
 *
 * Similar to the GOP based temporal filter of the HM reference encoder.
 * The pictures that the rest of the GOP is predicted from are replaced by a
 * motion compensated, bilateral weighted average of the picture and its
 * neighbouring pictures. This removes noise that would otherwise be coded
 * into the key pictures and propagated to the whole GOP.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"

/**
 * \brief Filter a picture with its neighbouring pictures.
 *
 * \param orig      picture to filter
 * \param refs      neighbouring pictures
 * \param offsets   display order distance of each neighbour (-2, -1, 1, 2, ...)
 * \param num_refs  number of neighbours
 * \param range     maximum distance used (the strength table depends on it)
 * \param strength  overall filter strength
 * \param qp        QP used for the filter strength
 * \return newly allocated filtered picture, or NULL on failure
 */
kvz_picture *kvz_temporal_filter(const kvz_picture *orig,
                                 const kvz_picture *const *refs,
                                 const int *offsets,
                                 int num_refs,
                                 int range,
                                 double strength,
                                 int qp);

#endif // TEMPORAL_FILTER_H_
