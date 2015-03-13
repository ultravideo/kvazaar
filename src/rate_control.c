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

#include "rate_control.h"

#include <math.h>

static const int SMOOTHING_WINDOW = 40;

/**
 * \brief Select a QP for encoding the next picture
 * \param state the main encoder state
 * \return the QP for the next picture, in range [0, 51]
 */
int8_t select_picture_QP(const encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  if (encoder->cfg->target_bitrate <= 0) {
    // Rate control disabled.
    return encoder->cfg->qp;
  }

  // At this point, total_bits_coded of the current state contains the
  // number of bits written encoder->owf frames before the current frame.
  const int bits_coded = state->global->total_bits_coded;
  const int pictures_coded = MAX(0, state->global->frame - encoder->owf);

  const double avg_bits_per_picture =
    encoder->cfg->target_bitrate / encoder->cfg->framerate;

  // TODO: use picture weights
  const double target_bits_current_picture =
    (avg_bits_per_picture * (pictures_coded + SMOOTHING_WINDOW) - bits_coded)
    / SMOOTHING_WINDOW;

  // TODO: take the picture headers into account
  const int pixels_per_picture = encoder->in.width * encoder->in.height;
  const double target_bits_per_pixel = target_bits_current_picture / pixels_per_picture;

  // The following magical constants, -5.7420835 and 18.598408755005686 are
  // based on the values given in
  //
  //   K. McCann et al., "High Effiency Video Coding (HEVC) Test Model 16
  //   (HM 16) Improved Encoder Description", JCTVC-S1002, October 2014,
  //   (p. 52 - 54)
  const int QP = (int)(-5.7420835 * log(MAX(target_bits_per_pixel, 0.001)) + 18.598408755005686);
  return CLIP(0, 51, QP);
}
