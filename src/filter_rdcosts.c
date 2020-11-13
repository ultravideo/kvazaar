#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define BUFSZ (64 * 64 * sizeof(uint16_t))
#define NUM_COEFF_BUCKETS (4)
#define NUM_OTHER_BUCKETS (0)
#define NUM_TOTAL_BUCKETS ((NUM_COEFF_BUCKETS) + (NUM_OTHER_BUCKETS))
#define MAX_COEFF_BUCKET  ((NUM_COEFF_BUCKETS) - 1)

#define clz(x) __builtin_clz(x)
#define ilog2(x) (sizeof(x) * 8 - clz(x) - 1)

void print_coeffs(const int16_t *buf, uint32_t size, uint32_t ccc)
{
  uint32_t i;
  printf("Buf size %u, ccc %u\n", size, ccc);
  for (i = 0; i < size; i++)
    printf("%i ", buf[i]);
  printf("\n");
}

void count_coeffs(const int16_t *buf, uint32_t size, uint64_t *buckets, uint64_t *num_signs, uint16_t *excess)
{
  *excess = 0;
  uint32_t i;

  for (i = 0; i < size; i++) {
    int16_t curr = buf[i];
    int16_t is_signed = curr >> 15;
    *num_signs += (is_signed & 1);

    uint16_t abs = (curr ^ is_signed) - is_signed;
    if (abs > MAX_COEFF_BUCKET) {
      *excess += abs - MAX_COEFF_BUCKET;
      abs = MAX_COEFF_BUCKET;
    }

    buckets[abs]++;
  }
}

void print_buckets(const uint64_t *buckets, uint64_t num_signs)
{
  uint32_t i;
  for (i = 0; i < NUM_COEFF_BUCKETS; i++)
    printf("%3u: %lu\n", i, buckets[i]);
  printf("Signs: %lu\n", num_signs);
}

void update_matrix(const uint64_t *buckets, uint64_t *mat)
{
  for (int y = 0; y < NUM_TOTAL_BUCKETS; y++) {
    for (int x = 0; x < NUM_TOTAL_BUCKETS; x++) {
      int curr_pos = y * NUM_TOTAL_BUCKETS + x;
      mat[curr_pos] += buckets[x] * buckets[y];
    }
  }
}

static inline int is_power_of_two(uint32_t u)
{
  return (u & (u - 1)) == 0;
}

int process_rdcosts(FILE *in, FILE *out)
{
  void *buf = malloc(BUFSZ);
  uint32_t *u32buf = (uint32_t *)buf;
  int16_t  *i16buf = (int16_t  *)buf;
  int rv = 0;

  float weights[NUM_TOTAL_BUCKETS] = {0.0f};

  uint64_t mat[NUM_TOTAL_BUCKETS * NUM_TOTAL_BUCKETS] = {0};

  while (!feof(in)) {
    uint32_t size, ccc, size_sqrt;
    uint64_t cg_buckets[NUM_TOTAL_BUCKETS] = {0};
    uint64_t cg_num_signs = 0;
    uint16_t excess = 0;
    size_t   n_read;

    n_read = fread(buf, sizeof(uint32_t), 2, in);
    size = u32buf[0];
    ccc  = u32buf[1];

    // Can't rely on feof() alone when reading from a pipe that might only get
    // closed long after the last data has been poured in
    if (n_read == 0) {
      break;
    }
    if (feof(in) || n_read < 2) {
      fprintf(stderr, "Unexpected EOF when reading header, managed still to read %u u32's\n", n_read);
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
    if (n_read != size) {
      fprintf(stderr, "Unexpected EOF when reading block, managed still to read %u i16's\n", n_read);
      rv = 1;
      goto out;
    }

    count_coeffs(i16buf, size, cg_buckets, &cg_num_signs, &excess);
    update_matrix(cg_buckets, mat);
  }
  for (int y = 0; y < NUM_TOTAL_BUCKETS; y++) {
    for (int x = 0; x < NUM_TOTAL_BUCKETS; x++) {
      int curr_pos = y * NUM_TOTAL_BUCKETS + x;
      printf("%lu ", mat[curr_pos]);
    }
    printf("\n");
  }
  fflush(stdout);

out:
  free(buf);
  return rv;
}

int main(int ar, char **av)
{
  return process_rdcosts(stdin, stdout);
}
