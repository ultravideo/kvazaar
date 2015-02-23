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

/**
 * \file
 * \brief Entry point for the Visual Studio project.
 *
 * This file is needed for Visual Studio, because it will not link the main
 * function from the .lib if the project has no .c files.
 *
 * \author Marko Viitanen ( fador@iki.fi ),
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ),
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

// This is not actually needed, because the linker will use the main from the
// .lib of the encoder, but I will leave it here in case we encounter some
// problem with that.
/*
int encmain(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  int i = 10;
  while (i) {
    --i;
  }
  return encmain(argc, argv);
}
*/
