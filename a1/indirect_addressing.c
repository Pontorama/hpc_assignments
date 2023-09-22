#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double drand(double limit) {
  return (double)rand() / RAND_MAX * limit - limit / 2;
}

double difftime_mus(struct timespec start_time, struct timespec stop_time) {
  return difftime(stop_time.tv_sec, start_time.tv_sec) * 1000000 +
         (stop_time.tv_nsec - start_time.tv_nsec) / 1000.0;
}

void vector_add(double *y, double a, double *x, int *index_vector, int size) {
  for (size_t k = 0; k < size; ++k) {
    size_t j = index_vector[k];
    y[j] += a * x[j];
  }
}

int main() {
  // Seed the random function
  srand(time(NULL));
  int iterations = 1000;
  // Initialize vectors
  size_t size = 1000000;
  int *p = (int *)malloc(sizeof(int) * size);
  for (size_t i = 0; i < size; i++)
    p[i] = i;

  double *x = (double *)malloc(sizeof(double) * size);
  double *y = (double *)malloc(sizeof(double) * size);

  for (int i = 0; i < size; i++) {
    x[i] = drand(2000);
    y[i] = drand(2000);
  }

  // Benchmarking
  struct timespec start_time;
  struct timespec stop_time;

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    vector_add(y, 3, x, p, size);
  }
  timespec_get(&stop_time, TIME_UTC);
  double difftime = difftime_mus(start_time, stop_time);

  int random_index = (int)(rand() / (1. * RAND_MAX) * size);

  printf("Linear assign index vector\n");
  printf("%f\n", y[random_index]);
  printf("Took %f mus\n", difftime);

  // Jump Initialization
  size_t size_jump = 1000;
  for (size_t j = 0, k = 0; j < size_jump; ++j) {
    for (size_t i = j; i < size; i += size_jump, ++k) {
      p[i] = k;
    }
  }

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    vector_add(y, 3, x, p, size);
  }
  timespec_get(&stop_time, TIME_UTC);

  difftime = difftime_mus(start_time, stop_time);

  random_index = (int)(rand() / (1. * RAND_MAX) * size);

  printf("Jump assign index vector\n");
  printf("%f\n", y[random_index]);
  printf("Took %f mus\n", difftime);

  // No index vector
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    for (int i = 0; i < size; i++) {
      y[i] += 3 * x[i];
    }
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = difftime_mus(start_time, stop_time);
  random_index = (int)(rand() / (1. * RAND_MAX) * size);
  printf("No index vector\n");
  printf("%f\n", y[random_index]);
  printf("Took %f mus\n", difftime);

  // Free mem
  free(y);
  free(x);
  free(p);
}
