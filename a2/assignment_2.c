#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static inline float get_distance(char *point_1, char *point_2) {
  // Assumes points in R^3
  // Convert from text to floats
  float point_1_num[3];
  float point_2_num[3];

  sscanf(point_1, "%f %f %f\n", &point_1_num[0], &point_1_num[1],
         &point_1_num[2]);
  sscanf(point_2, "%f %f %f\n", &point_2_num[0], &point_2_num[1],
         &point_2_num[2]);
  float x = point_2_num[0] - point_1_num[0];
  float y = point_2_num[1] - point_1_num[1];
  float z = point_2_num[2] - point_1_num[2];
  return sqrtf(x * x + y * y + z * z);
}

int main(int argc, char **argv) {
  // Read cmd line inputs here
  unsigned short n_threads = atoi(&argv[1][2]);
  // Open file, get file length
  char filename[10] = "cells";
  FILE *file = fopen(filename, "r");
  fseek(file, 0L, SEEK_END);
  const unsigned short line_size = 24;
  const unsigned int file_size = ftell(file) / line_size;
  rewind(file);

  omp_set_num_threads(n_threads);

  unsigned int chunk_size = file_size / n_threads;
  static const unsigned int MAX_CHUNK_SIZE = 600000;
  if (chunk_size > MAX_CHUNK_SIZE)
    chunk_size = MAX_CHUNK_SIZE;
  if (chunk_size < file_size)
    chunk_size = file_size;

  // TODO: Duplicate histogram into as many pieces as there are threads
  const unsigned short max_distance = 3465;
  unsigned long long int *distances = (unsigned long long int *)malloc(
      sizeof(unsigned long long int) * max_distance);

  char format[10];
  sprintf(format, "%%%ic", line_size * chunk_size);

  // Init chunk buffers
  char *current_chunk_str =
      (char *)malloc(sizeof(char) * line_size * chunk_size);
  char *other_chunk_str = (char *)malloc(sizeof(char) * line_size * chunk_size);

  unsigned int n_chunks = file_size / chunk_size;

  for (unsigned int i = 0; i < n_chunks; i++) {
    // Read first chunk
    fscanf(file, format, current_chunk_str);
// Read next chunk(s) and get distances
#pragma omp parallell for
    for (unsigned int j = i; j < n_chunks; j++) {
      if (i == j) {
        fseek(file, i * chunk_size * line_size, SEEK_SET);
      }
      fscanf(file, format, other_chunk_str);
      // Get distances
      for (unsigned int k = 0; k < chunk_size; k++) {
        for (unsigned int l = k + 1; l < chunk_size; l += 4) {
          short distance_1 =
              (unsigned short)(100 *
                               (get_distance(current_chunk_str + k * line_size,
                                             other_chunk_str + l * line_size) +
                                0.005));
          short distance_2 =
              (unsigned short)(100 *
                               (get_distance(current_chunk_str + k * line_size,
                                             other_chunk_str +
                                                 (l + 1) * line_size) +
                                0.005));
          short distance_3 =
              (unsigned short)(100 *
                               (get_distance(current_chunk_str + k * line_size,
                                             other_chunk_str +
                                                 (l + 2) * line_size) +
                                0.005));
          short distance_4 =
              (unsigned short)(100 *
                               (get_distance(current_chunk_str + k * line_size,
                                             other_chunk_str +
                                                 (l + 3) * line_size) +
                                0.005));
          distances[distance_1] += 1;
          distances[distance_2] += 1;
          distances[distance_3] += 1;
          distances[distance_4] += 1;
        }
      }
    }
  }

  free(current_chunk_str);
  free(other_chunk_str);
  for (size_t i = 0; i < max_distance; i++) {
    if (distances[i] != 0)
      printf("%f: %llu\n", i / 100.0, distances[i]);
  }
  free(distances);
}
