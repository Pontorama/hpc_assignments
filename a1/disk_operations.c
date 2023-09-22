#include "../lib/benchmark_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void write_to_disk(FILE *file) {
  int big_num = 1;
  for (int i = 0; i < 20; i++) {
    big_num *= 2;
    fwrite((void *)&big_num, sizeof(int), 1, file);
    fflush(file);
  }
}

void read_from_disk(FILE *file) {
  int big_num;
  for (int i = 0; i < 20; i++) {
    fread((void *)&big_num, sizeof(int), 1, file);
  }
}

int main() {
  int iterations = 10;
  // Open file for use
  FILE *file = fopen("testfile.dat", "w");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  // Benchmark writing
  struct timespec start_time;
  struct timespec stop_time;

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    write_to_disk(file);
  }
  timespec_get(&stop_time, TIME_UTC);
  double difftime = get_timediff_mus(start_time, stop_time);

  fclose(file);

  printf("Writing to file, individual write\n");
  printf("Took %f mus\n", difftime);

  // Benchmark reading
  file = fopen("testfile.dat", "r");
  if (file == NULL) {
    printf("Error opening file\n");
    return -1;
  }

  timespec_get(&start_time, TIME_UTC);
  for (int iter = 0; iter < iterations; iter++) {
    read_from_disk(file);
  }
  timespec_get(&stop_time, TIME_UTC);
  difftime = get_timediff_mus(start_time, stop_time);
  printf("Reading from file, individual read\n");
  printf("Took %f mus\n", difftime);

  fclose(file);

  return 0;
}
