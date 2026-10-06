#ifndef LOOKAHEAD_H_
#define LOOKAHEAD_H_
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
 * \brief Lookahead analysis of future input frames.
 *
 * Input pictures are queued for cfg.lookahead frames before they are given
 * to the encoder. Each queued picture is downscaled to half resolution and
 * analysed with a cheap intra and motion search on 8x8 low resolution blocks
 * (16x16 in full resolution). The results are used for
 *  - scene cut detection (an IDR is inserted at the cut),
 *  - CU-tree: temporal propagation of block importance into per-CTU QP
 *    offsets,
 *  - CRF: complexity based frame QP selection,
 *  - ABR rate control (--rc-algorithm abr): complexity based bits model.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"

struct encoder_control_t;

/**
 * \brief Results of the lookahead analysis for one frame.
 *
 * Created when the frame leaves the lookahead and handed over to the encoder
 * state that encodes the frame (see kvz_lookahead_take_result).
 */
typedef struct kvz_la_result {
  const kvz_picture *pic;   //!< picture the result belongs to (lookup key)
  int64_t index;            //!< display order index in the lookahead
  bool scenecut;            //!< frame starts a new scene (encoded as IDR)
  double cost;              //!< inter (P) cost of the frame per low-res block
  double intra_cost;        //!< intra cost of the frame per low-res block
  double blurred_cost;      //!< cost blurred over neighbouring frames
  double crf_qp;            //!< base QP chosen by the CRF model (if CRF)
  double *ctu_offsets;      //!< CU-tree QP offset per CTU, NULL if not used
  double *qg_offsets;       //!< CU-tree QP offset per 32x32 block, NULL if not used
  kvz_picture *orig;        //!< original picture if pic is temporally filtered
  /**
   * \brief Low resolution motion of each 16x16 block (x, y in full
   *        resolution pixels per picture). mvs: motion to the previous
   *        picture, next_mvs: motion of the next picture to this one.
   *        NULL if not available.
   */
  int16_t *mvs;
  int16_t *next_mvs;
  int mv_bw, mv_bh;         //!< size of the motion fields in blocks
  struct kvz_la_result *next;
} kvz_la_result;

typedef struct kvz_lookahead kvz_lookahead;

kvz_lookahead *kvz_lookahead_alloc(const struct encoder_control_t *encoder);
void kvz_lookahead_free(kvz_lookahead *la);

/** \brief Number of pictures currently queued. */
int kvz_lookahead_count(const kvz_lookahead *la);

/** \brief Queue depth after which a picture is released. */
int kvz_lookahead_depth(const kvz_lookahead *la);

/**
 * \brief Add an input picture to the lookahead and analyse it.
 *
 * Takes a new reference to the picture.
 */
void kvz_lookahead_push(kvz_lookahead *la, kvz_picture *pic);

/** \brief Whether the next picture to be released starts a new scene. */
bool kvz_lookahead_next_is_scenecut(kvz_lookahead *la);

/**
 * \brief Release the oldest picture.
 *
 * Runs the final analysis (scene cut, CU-tree, CRF) for the oldest picture
 * and returns it. The caller owns the returned reference. The analysis
 * result is kept until kvz_lookahead_take_result is called.
 */
kvz_picture *kvz_lookahead_pop(kvz_lookahead *la);

/**
 * \brief Take the analysis result of a released picture.
 *
 * Ownership of the result is moved to the caller, which must free it with
 * kvz_lookahead_free_result. Returns NULL if no result exists.
 */
kvz_la_result *kvz_lookahead_take_result(kvz_lookahead *la, const kvz_picture *pic);
void kvz_lookahead_free_result(kvz_la_result *res);

/** \brief Display index of the newest analysed frame, -1 if none. */
int64_t kvz_lookahead_newest(const kvz_lookahead *la);

/**
 * \brief Complexity of an analysed frame.
 *
 * \param cost   inter (P) cost per low-res block
 * \param intra  intra cost per low-res block
 * \param key    whether the frame is expected to be coded as an intra picture
 * \return false if the frame is not available
 */
bool kvz_lookahead_frame_info(const kvz_lookahead *la, int64_t index,
                              double *cost, double *intra, bool *key);

#endif // LOOKAHEAD_H_
