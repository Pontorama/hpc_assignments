#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
  unsigned long sum = 0;
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  double bench_diff_time;

  timespec_get(&bench_start_time, TIME_UTC);
  for (size_t i = 0; i < 1000000000; i++) {
    sum += i;
  }
  timespec_get(&bench_stop_time, TIME_UTC);

  printf("The sum is %lu\n", sum);
  bench_diff_time =
      difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec) * 1000000 +
      (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) / 1000.0;
  printf("Benchmark time: %fmus\n", bench_diff_time);
}
