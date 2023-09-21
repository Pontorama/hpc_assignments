#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double drand(double limit) {
  return (double)rand() / RAND_MAX * limit - limit / 2;
}

double **create_matrix(size_t n_rows, size_t n_cols) {
  double **a = (double **)malloc(sizeof(double *) * n_rows);
  double *a_entries = (double *)malloc(sizeof(double) * n_rows * n_cols);
  for (size_t i = 0, j = 0; i < n_rows; ++i, j += n_cols) {
    a[i] = a_entries + j;
  }

  // init to 0
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      a[i][j] = 0;
    }
  }
  return a;
}

double difftime_mus(struct timespec start_time, struct timespec stop_time) {
  return difftime(stop_time.tv_sec, start_time.tv_sec) * 1000000 +
         (stop_time.tv_nsec - start_time.tv_nsec) / 1000.0;
}

void row_sums(double *sums, const double **matrix, size_t nrs, size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum = 0;
    for (size_t j = 0; j < ncs; j++) {
      sum += matrix[i][j];
    }
    *sums += sum;
  }
}

void row_sums_unrolled2(double *sums, const double **matrix, size_t nrs,
                        size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum0 = 0;
    double sum1 = 0;
    for (size_t j = 0; j < ncs; j += 2) {
      sum0 += matrix[i][j];
      sum1 += matrix[i][j];
    }
    *sums = sum0 + sum1;
  }
}

int main(int argc, char **argv) {
  /*
  // Get cmd line args
  int iterations = 5000;
  int opt;
  while ((opt = getopt(argc, argv, "i:")) != 1) {
    switch (opt) {
    case 'i':
      iterations = (int)atoi(optarg);
      break;
    default:
      abort();
    }
  }
  */
  int iterations = 5000;
  // Seed the random function
  srand(time(NULL));
  // Create (contigous) matrix
  size_t size = 1000;
  double **a = create_matrix(size, size);
  // assign random values
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      a[i][j] = drand(2000);
    }
  }

  // Create vectors for the sums
  double **regular_sums = create_matrix(iterations, size);
  struct timespec start_time;
  struct timespec stop_time;

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; ++iter) {
    for (size_t i = 0; i < size; ++i) {
      row_sums(regular_sums[iter], (const double **)a, size, size);
    }
  }
  timespec_get(&stop_time, TIME_UTC);
  double difftime = difftime_mus(start_time, stop_time);
  int random_index = (int)((rand() / (1.0 * RAND_MAX)) * size);
  int random_iter_index = (int)((rand() / (1.0 * RAND_MAX)) * iterations);
  printf("%i %i\n", random_index, random_iter_index);
  printf("Regular sum:\n");
  printf("Index %i is %f\n", random_index,
         regular_sums[random_iter_index][random_index]);
  printf("Time per iteration (%i total): %f\n", iterations,
         difftime / (1.0 * iterations));

  // Free memory
  free(a[0]);
  free(a);
  free(regular_sums[0]);
  free(regular_sums);

  return 0;
}
