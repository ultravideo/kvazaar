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
