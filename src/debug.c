/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "debug.h"

#include <stdio.h>


/**
 * Return a FILE* to a file that can be used with render_cu_file and close_cu_file.
 */
FILE * open_cu_file(char *filename) {
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "<html><head><link rel='stylesheet' type='text/css' href='cu_style.css' /></head><body>");
  return fp;
}

/**
 * Close the FILE* returned by open_cu_file.
 */
void close_cu_file(FILE *fp) {
  fprintf(fp, "</body></html>");
  fclose(fp);
}

/**
 * Print information about the Coding Unit (CU) into the FILE* provided by open_cu_file.
 */
unsigned render_cu_file(encoder_control *encoder, unsigned depth, uint16_t xCtb, uint16_t yCtb, FILE *fp)
{
  CU_info *cu = &encoder->in.cur_pic->CU[depth][xCtb + yCtb * (encoder->in.width_in_lcu<<MAX_DEPTH)];
  unsigned lambda_cost = (4 * g_lambda_cost[encoder->QP]) << 4;
  unsigned sum = 0;
  unsigned best_cost = -1;
  char type = cu->type == CU_INTRA ? 'I' : 'P';

  fprintf(fp, 
    "\n<table class=d%u><tr><td colspan=2>"
    "%u (%u, %u), %d, %c, "
    "c=%u, mv=(%d, %d)</td></tr>\n", 
    depth,
    depth, xCtb, yCtb, cu->split, (cu->type == CU_INTRA ? 'I' : 'P'),
    cu->inter.cost, cu->inter.mv[0], cu->inter.mv[1]);


  if(depth != MAX_INTER_SEARCH_DEPTH)
  {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);

    fprintf(fp, "<tr><td>");
    sum += render_cu_file(encoder, depth + 1, xCtb, yCtb, fp);
    fprintf(fp, "</td><td>");
    sum += render_cu_file(encoder, depth + 1, xCtb + change, yCtb, fp);
    fprintf(fp, "</td></tr>");

    fprintf(fp, "<tr><td>");
    sum += render_cu_file(encoder, depth + 1, xCtb, yCtb + change, fp);
    fprintf(fp, "</td><td>");
    sum += render_cu_file(encoder, depth + 1, xCtb + change, yCtb + change, fp);
    fprintf(fp, "</td></tr>");

    fprintf(fp, "<tr><td colspan=2>sum=%u, sum+lambda=%u</td></tr>",
      sum, sum + lambda_cost);
    if (sum + lambda_cost < cu->inter.cost) {
      best_cost = sum + lambda_cost;
    } else {
      best_cost = cu->inter.cost;
    }
  } else {
    best_cost = cu->inter.cost;
  }

  fprintf(fp, "</table>");
  return best_cost;
}
