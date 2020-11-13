#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define BUFSZ (64 * 64 * sizeof(uint16_t))
#define NUM_COEFF_BUCKETS (4)
#define NUM_OTHER_BUCKETS (0)
#define NUM_TOTAL_BUCKETS ((NUM_COEFF_BUCKETS) + (NUM_OTHER_BUCKETS))
#ifdef ERR_SQUARED
#define STEPSIZE (0.00000001f * 0.000001f)
#else
#define STEPSIZE (0.00000001f)
#endif

#define clz(x) __builtin_clz(x)
#define ilog2(x) (sizeof(x) * 8 - clz(x) - 1)
#define coord(x,y,w) ((x)+((y)*(w)))

void update_result(const uint64_t *buckets, uint64_t ccc, const double *mat, double *res)
{
  for (int y = 0; y < NUM_TOTAL_BUCKETS; y++) {
    double addend = 0.0;
    for (int x = 0; x < NUM_TOTAL_BUCKETS; x++) {
      addend += mat[coord(x, y, NUM_TOTAL_BUCKETS)] * (double)buckets[x];
    }
    addend *= (double)ccc;
    res[y] += addend;
  }
}

void read_matrix(const char *fn, double *mat)
{
  FILE *f = fopen(fn, "r");
  for (int y = 0; y < NUM_TOTAL_BUCKETS; y++) {
    for (int x = 0; x < NUM_TOTAL_BUCKETS; x++) {
      float curr;
      fscanf(f, "%f", &curr);
      mat[x + y * NUM_TOTAL_BUCKETS] = curr;
    }
  }
  fclose(f);
}

void count_coeffs(const int16_t *buf, uint32_t size, uint64_t *buckets, uint64_t *num_signs)
{
  uint32_t i;
  for (i = 0; i < size; i++) {
    int16_t curr = buf[i];
    int16_t is_signed = curr >> 15;
    *num_signs += (is_signed & 1);

    uint16_t abs = (curr ^ is_signed) - is_signed;
    if (abs >= NUM_COEFF_BUCKETS)
      abs = NUM_COEFF_BUCKETS - 1;

    buckets[abs]++;
  }
}

static inline int is_power_of_two(uint32_t u)
{
  return (u & (u - 1)) == 0;
}

int process_rdcosts(FILE *in, FILE *out, const double *mat)
{
  void *buf = malloc(BUFSZ);
  uint32_t *u32buf = (uint32_t *)buf;
  int16_t  *i16buf = (int16_t  *)buf;
  int rv = 0;

  double res[NUM_TOTAL_BUCKETS] = {0.0};

  while (!feof(in)) {
    uint32_t size, ccc, size_sqrt;
    uint64_t cg_buckets[NUM_TOTAL_BUCKETS] = {0};
    uint64_t cg_num_signs = 0;
    size_t   n_read;

    n_read = fread(buf, sizeof(uint32_t), 2, in);
    size = u32buf[0];
    ccc  = u32buf[1];

    // Can't rely on feof() alone when reading from a pipe that might only get
    // closed long after the last data has been poured in
    if (n_read == 0) {
      break;
    }
    if (feof(in) || n_read < sizeof(uint32_t) * 2) {
      fprintf(stderr, "Unexpected EOF when reading header, managed still to read %u bytes\n", n_read);
      rv = 1;
      goto out;
    }
    if (!is_power_of_two(size)) {
      fprintf(stderr, "Errorneous block size %u\n", size);
      rv = 1;
      goto out;
    }

    size_sqrt = 1 << (ilog2(size) >> 1);

    n_read = fread(buf, sizeof(int16_t), size, in);
    if (n_read != size * sizeof(int16_t)) {
      fprintf(stderr, "Unexpected EOF when reading block, managed still to read %u bytes\n", n_read);
      rv = 1;
      goto out;
    }

    count_coeffs(i16buf, size, cg_buckets, &cg_num_signs);
    update_result(cg_buckets, ccc, mat, res);
  }

  for (int y = 0; y < NUM_TOTAL_BUCKETS; y++)
    fprintf(out, "%g\n", (float)(res[y]));

out:
  free(buf);
  return rv;
}

int main(int ar, char **av)
{
  double mat[NUM_TOTAL_BUCKETS * NUM_TOTAL_BUCKETS] = {0.0};
  if (ar != 2) {
    fprintf(stderr, "gib matrix plz\n");
    return 1;
  }
  read_matrix(av[1], mat);
  return process_rdcosts(stdin, stdout, mat);
}

