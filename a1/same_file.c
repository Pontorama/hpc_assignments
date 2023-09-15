#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void mul_cpx(double *a_re, double *a_im, double *b_re, double *b_im,
             double *c_re, double *c_im);

void mul_cpx(double *a_re, double *a_im, double *b_re, double *b_im,
             double *c_re, double *c_im) {
  *a_re = *b_re * *c_re - *b_im * *c_im;
  *a_im = *b_re * *c_im + *c_re * *b_im;
}

double drand(double limit) {
  return (double)rand() / RAND_MAX * limit - limit / 2;
}

int main() {
  // Create new seed for random numbers
  srand(time(NULL));

  int size = 30000;
  double *as_re = (double *)malloc(sizeof(double) * size);
  double *as_im = (double *)malloc(sizeof(double) * size);
  double *bs_re = (double *)malloc(sizeof(double) * size);
  double *bs_im = (double *)malloc(sizeof(double) * size);
  double *cs_re = (double *)malloc(sizeof(double) * size);
  double *cs_im = (double *)malloc(sizeof(double) * size);

  // Generate entries
  for (uint i = 0; i < size; i++) {
    as_re[i] = 0;
    as_im[i] = 0;
    bs_re[i] = drand(50);
    bs_im[i] = drand(50);
    cs_re[i] = drand(50);
    cs_im[i] = drand(50);
  }
  // Benchmark
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  timespec_get(&bench_start_time, TIME_UTC);
  uint bench_iterations = 1000000;
  for (uint j = 0; j < bench_iterations; j++) {
    for (uint i = 0; i < size; i++) {
      mul_cpx(&as_re[i], &as_im[i], &bs_re[i], &bs_im[i], &cs_re[i], &cs_im[i]);
    }
  }
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_diff_time =
      difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec) * 1000000 +
      (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) / 1000.0;

  uint rand_index = (int)((double)rand() / (RAND_MAX)*size);
  printf("Index %i is %f + i %f\n", rand_index, as_re[rand_index],
         as_im[rand_index]);

  printf("Benchmark time: %fmus\n", bench_diff_time);

  return 0;
}
