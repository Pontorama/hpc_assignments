#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double drand(double limit) {
  return (double)rand() / RAND_MAX * limit - limit / 2;
}

void row_sums(double *sums, const double **matrix, size_t nrs, size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum = 0;
    for (size_t j = 0; j < ncs; j++) {
      sum += matrix[i][j];
    }
    sums[i] = sum;
  }
}

void col_sums(double *sums, const double **matrix, size_t nrs, size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum = 0;
    for (size_t j = 0; j < ncs; j++) {
      sum += matrix[j][i];
    }
    sums[i] = sum;
  }
}

int main() {
  // Seed random number
  srand(time(NULL));
  // Allocate (contigous) memory
  uint size = 1000;
  double *asentries = (double *)malloc(sizeof(double) * size * size);
  double **as = (double **)malloc(sizeof(double *) * size);

  for (uint i = 0, j = 0; i < size; i++, j += size) {
    as[i] = asentries + j;
  }

  for (uint i = 0; i < size; ++i) {
    for (uint j = 0; j < size; ++j) {
      as[i][j] = drand(2000);
    }
  }

  // Compute row/col sums
  uint iterations = 5000;
  // Row sum
  struct timespec start_time_row;
  struct timespec stop_time_row;

  double **all_sums = (double **)malloc(sizeof(double *) * iterations);
  double *sums = (double *)malloc(sizeof(double) * size * iterations);
  for (uint i = 0, j = 0; i < iterations; ++i, j += size) {
    all_sums[i] = sums + j;
  }

  for (uint i = 0; i < iterations; i++) {
    for (uint j = 0; j < size; j++) {
      all_sums[i][j] = 0;
    }
  }
  timespec_get(&start_time_row, TIME_UTC);
  for (uint iter = 0; iter < iterations; iter++) {
    row_sums(all_sums[iter], (const double **)as, size, size);
  }
  timespec_get(&stop_time_row, TIME_UTC);

  double diff_time_row =
      difftime(stop_time_row.tv_sec, start_time_row.tv_sec) * 1e6 +
      (stop_time_row.tv_nsec - start_time_row.tv_nsec) / 1000.0;
  int rand_index = (int)((double)rand() / (RAND_MAX * iterations));
  int rand_index_row = (int)((double)rand() / (RAND_MAX * size));
  printf("(Row) sum is %f\n", all_sums[rand_index][rand_index_row]);
  printf("Time: %fmus\n", diff_time_row);

  // Free memory
  free(sums);
  free(all_sums);

  // Col sum
  struct timespec start_time_col;
  struct timespec stop_time_col;

  all_sums = (double **)malloc(sizeof(double *) * iterations);
  sums = (double *)malloc(sizeof(double) * size * iterations);
  for (uint i = 0, j = 0; i < iterations; ++i, j += size) {
    all_sums[i] = sums + j;
  }

  for (uint i = 0; i < iterations; i++) {
    for (uint j = 0; j < size; j++) {
      all_sums[i][j] = 0;
    }
  }
  timespec_get(&start_time_col, TIME_UTC);
  for (uint iter = 0; iter < iterations; iter++) {
    col_sums(all_sums[iter], (const double **)as, size, size);
  }
  timespec_get(&stop_time_col, TIME_UTC);

  double diff_time_col =
      difftime(stop_time_col.tv_sec, start_time_col.tv_sec) * 1e6 +
      (stop_time_col.tv_nsec - start_time_col.tv_nsec) / 1000.0;
  rand_index = (int)((double)rand() / (RAND_MAX * iterations));
  rand_index_row = (int)((double)rand() / (RAND_MAX * size));
  printf("(Col) sum is %f\n", all_sums[rand_index][rand_index_row]);
  printf("Time: %fmus\n", diff_time_col);

  // Free memory
  free(sums);
  free(all_sums);
  free(asentries);
  free(as);
}
