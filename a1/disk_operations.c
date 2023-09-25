#include "../lib/benchmark_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
  char filename[13] = "testfile.dat\0";
  int iterations = 10;
  // Prepare numbers
  int big_num = 1;
  for (size_t i = 0; i < 20; i++) {
    big_num *= 2;
  }
  // Open file for use
  FILE *file = fopen(filename, "w");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  // Benchmark writing
  struct timespec start_time;
  struct timespec stop_time;

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    for (int i = 0; i < big_num; i++) {
      fwrite((void *)&i, sizeof(int), 1, file);
      fflush(file);
    }
  }
  timespec_get(&stop_time, TIME_UTC);
  double difftime = get_timediff_mus(start_time, stop_time);

  fclose(file);

  printf("Writing to file, individual write\n");
  printf("Took %f mus\n", difftime);

  // Benchmark vector write
  file = fopen(filename, "w");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  int *vector = (int *)malloc(sizeof(int) * big_num);
  for (int i = 0; i < big_num; i++) {
    vector[i] = i;
  }
  timespec_get(&start_time, TIME_UTC);
  for (size_t iter = 0; iter < iterations; iter++) {
    fwrite((void *)vector, sizeof(int), big_num, file);
    fflush(file);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = get_timediff_mus(start_time, stop_time);
  fclose(file);
  free(vector);

  printf("Writing to file, vector\n");
  printf("Took %f mus\n", difftime);

  // Benchmark reading
  file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  int read_big_num;
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    fread(&read_big_num, sizeof(int), 1, file);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = get_timediff_mus(start_time, stop_time);
  printf("Reading from file, individual read\n");
  printf("Took %f mus\n", difftime);

  fclose(file);

  // Benchmark vector read
  file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  vector = (int *)malloc(sizeof(int) * big_num);
  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    fread(vector, sizeof(int), big_num, file);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = get_timediff_mus(start_time, stop_time);
  printf("Reading from file, individual read\n");
  printf("Took %f mus\n", difftime);

  fclose(file);

  return 0;
}
