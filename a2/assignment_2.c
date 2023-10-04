#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static inline float get_distance(float *point_1, float *point_2) {
  // Assumes points in R^3
  return sqrtf((point_2[0] - point_1[0]) * (point_2[0] - point_1[0]) +
               (point_2[1] - point_1[1]) * (point_2[1] - point_1[1]) +
               (point_2[2] - point_1[2]) * (point_2[2] - point_1[2]));
}

int main(int argc, char **argv) {
  // Read cmd line inputs here

  // Open file, get file length
  char filename[10] = "cells.txt";
  FILE *file = fopen(filename, "r");
  fseek(file, 0L, SEEK_END);
  const unsigned short line_size = 24;
  const unsigned int file_size = ftell(file) / line_size;
  rewind(file);

  omp_set_num_threads(4);

  const unsigned int chunk_size = 10;

  float *all_cell_entries = (float *)malloc(sizeof(float) * chunk_size * 3);
  float **all_cells = (float **)malloc(sizeof(float *) * chunk_size);
  for (short i = 0; i < chunk_size; i++) {
    all_cells[i] = all_cell_entries + i * 3;
  }

  for (size_t i = 0; i < chunk_size; i++) {
    for (size_t j = 0; j < 3; j++) {
      all_cells[i][j] = 0;
    }
  }

  const unsigned short max_distance = 3465;
  unsigned long long int *distances = (unsigned long long int *)malloc(
      sizeof(unsigned long long int) * max_distance);

  unsigned int current_cell_point = 0;
  float *current_cell = (float *)malloc(sizeof(float) * 3);
  char *chunk_string =
      (char *)malloc(sizeof(char) * (chunk_size * line_size + 1));

#pragma omp parallell shared(all_cells, distances, current_cell, chunk_string, \
                                 current_cell_point)
  for (size_t i = 0; i < file_size; i++) {
    fseek(file, current_cell_point, SEEK_SET);
    fscanf(file, "%f %f %f\n", &current_cell[0], &current_cell[1],
           &current_cell[2]);
    current_cell_point = ftell(file);
    fgets(chunk_string, chunk_size * line_size + 1, file);

    for (size_t j = 0; j < (file_size - current_cell_point) / chunk_size; ++j) {
// Loop over all whole chunks
#pragma omp for
      for (size_t k = 0; k < chunk_size; k++) {
        sscanf(chunk_string + k * line_size, "%f %f %f\n", &all_cells[k][0],
               &all_cells[k][1], &all_cells[k][2]);
        short distance =
            (short)(100 * get_distance(current_cell, all_cells[k]));
#pragma omp atomic
        distances[distance] += 1;
      }
    }
#pragma omp for
    for (size_t j = 0; j < file_size % chunk_size; ++j) {
      // Loop over remaining lines
      sscanf(chunk_string + j * line_size, "%f %f %f\n", &all_cells[j][0],
             &all_cells[j][1], &all_cells[j][2]);
      short distance = (short)(100 * get_distance(current_cell, all_cells[j]));
#pragma omp atomic
      distances[distance] += 1;
    }
  }
  free(current_cell);
  free(all_cell_entries);
  free(all_cells);

  for (size_t i = 0; i < max_distance; i++) {
    if (distances[i] != 0)
      printf("%f: %llu\n", i / 100.0, distances[i]);
  }
  free(distances);
}
