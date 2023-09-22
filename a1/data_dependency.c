#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double drand(double limit) {
  return (double)rand() / (1.0 * RAND_MAX) * limit - limit / 2.0;
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
    sums[i] = sum;
  }
}

void row_sums_unrolled2(double *sums, const double **matrix, size_t nrs,
                        size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum0 = 0;
    double sum1 = 0;
    for (size_t j = 0; j < ncs; j += 2) {
      sum0 += matrix[i][j];
      sum1 += matrix[i][j + 1];
    }
    sums[i] = sum0 + sum1;
  }
}

void row_sums_unrolled4(double *sums, const double **matrix, size_t nrs,
                        size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    for (size_t j = 0; j < ncs; j += 4) {
      sum0 += matrix[i][j];
      sum1 += matrix[i][j + 1];
      sum2 += matrix[i][j + 2];
      sum3 += matrix[i][j + 3];
    }
    sums[i] = sum0 + sum1 + sum2 + sum3;
  }
}

void row_sums_unrolled4_vector(double *sums, const double **matrix, size_t nrs,
                               size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum[4] = {0, 0, 0, 0};
    for (size_t j = 0; j < ncs; j += 4) {
      sum[0] += matrix[i][j];
      sum[1] += matrix[i][j + 1];
      sum[2] += matrix[i][j + 2];
      sum[3] += matrix[i][j + 3];
    }
    sums[i] = sum[0] + sum[1] + sum[2] + sum[3];
  }
}

void row_sums_unrolled8(double *sums, const double **matrix, size_t nrs,
                        size_t ncs) {
  for (size_t i = 0; i < nrs; i++) {
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    double sum7 = 0;
    for (size_t j = 0; j < ncs; j += 8) {
      sum0 += matrix[i][j];
      sum1 += matrix[i][j + 1];
      sum2 += matrix[i][j + 2];
      sum3 += matrix[i][j + 3];
      sum4 += matrix[i][j + 4];
      sum5 += matrix[i][j + 5];
      sum6 += matrix[i][j + 6];
      sum7 += matrix[i][j + 7];
    }
    sums[i] = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
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
      break;
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

  // Regular sum
  // Create vectors for the sums
  double **regular_sums = create_matrix(iterations, size);
  struct timespec start_time;
  struct timespec stop_time;

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; ++iter) {
    row_sums(regular_sums[iter], (const double **)a, size, size);
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
  free(regular_sums[0]);
  free(regular_sums);

  // Unrolled2 sum
  double **unrolled_2_sum = create_matrix(iterations, size);
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    row_sums_unrolled2(unrolled_2_sum[iter], (const double **)a, size, size);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = difftime_mus(start_time, stop_time);

  random_index = (int)((rand() / (1.0 * RAND_MAX)) * size);
  random_iter_index = (int)((rand() / (1.0 * RAND_MAX)) * iterations);
  printf("Unrolled2 sum:\n");
  printf("Index %i is %f\n", random_index,
         unrolled_2_sum[random_iter_index][random_index]);
  printf("Time per iteration (%i total): %f\n", iterations,
         difftime / (1.0 * iterations));
  free(unrolled_2_sum[0]);
  free(unrolled_2_sum);

  // Unrolled4 sum
  double **unrolled_4_sum = create_matrix(iterations, size);
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    row_sums_unrolled4(unrolled_4_sum[iter], (const double **)a, size, size);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = difftime_mus(start_time, stop_time);

  random_index = (int)((rand() / (1.0 * RAND_MAX)) * size);
  random_iter_index = (int)((rand() / (1.0 * RAND_MAX)) * iterations);
  printf("Unrolled4 sum:\n");
  printf("Index %i is %f\n", random_index,
         unrolled_4_sum[random_iter_index][random_index]);
  printf("Time per iteration (%i total): %f\n", iterations,
         difftime / (1.0 * iterations));
  free(unrolled_4_sum[0]);
  free(unrolled_4_sum);

  // Unrolled8 sum
  double **unrolled_8_sum = create_matrix(iterations, size);
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    row_sums_unrolled8(unrolled_8_sum[iter], (const double **)a, size, size);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = difftime_mus(start_time, stop_time);

  random_index = (int)((rand() / (1.0 * RAND_MAX)) * size);
  random_iter_index = (int)((rand() / (1.0 * RAND_MAX)) * iterations);
  printf("Unrolled8 sum:\n");
  printf("Index %i is %f\n", random_index,
         unrolled_8_sum[random_iter_index][random_index]);
  printf("Time per iteration (%i total): %f\n", iterations,
         difftime / (1.0 * iterations));
  free(unrolled_8_sum[0]);
  free(unrolled_8_sum);

  // Unrolled4_vector sum
  double **unrolled_4_vector_sum = create_matrix(iterations, size);
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    row_sums_unrolled4_vector(unrolled_4_vector_sum[iter], (const double **)a,
                              size, size);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = difftime_mus(start_time, stop_time);

  random_index = (int)((rand() / (1.0 * RAND_MAX)) * size);
  random_iter_index = (int)((rand() / (1.0 * RAND_MAX)) * iterations);
  printf("Unrolled4 (vector implementation) sum:\n");
  printf("Index %i is %f\n", random_index,
         unrolled_4_vector_sum[random_iter_index][random_index]);
  printf("Time per iteration (%i total): %f\n", iterations,
         difftime / (1.0 * iterations));
  free(unrolled_4_vector_sum[0]);
  free(unrolled_4_vector_sum);

  // Free memory
  free(a[0]);
  free(a);

  return 0;
}
