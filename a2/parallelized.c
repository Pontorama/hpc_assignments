/*
 * Compute a histogram of distances between points from a given file
 * Constraints:
 * - Coordinates between -10 & 10
 * - Program must consume less than 5MiB = 5* 1024^2 bytes of heap memory
 * - Maximum number of points = 2^32 = 4294967296
 *
 * */
#include "../lib/matrix_utils.h"
#include "read_text_from_file.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

float get_distance(float *point_1, float *point_2) {
  // Assumes points in R^3
  return sqrtf((point_2[0] - point_1[0]) * (point_2[0] - point_1[0]) +
               (point_2[1] - point_1[1]) * (point_2[1] - point_1[1]) +
               (point_2[2] - point_1[2]) * (point_2[2] - point_1[2]));
}

int main(int argc, char **argv) {
  char filename[10] = "cells.txt";
  short line_size = 24; // 24 characters per line in input file
  short point_char_size =
      8; // 8 characters per point in file (including spaces)

  // Init omp
  omp_set_num_threads(5);

  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    printf("Could not open file %s\n", filename);
    return -1;
  }

  // Split file into regions and then parallelize
  fseek(file, 0L, SEEK_END);
  unsigned int file_size = ftell(file);
  printf("File size: %i\n", file_size);
  // Convert number of characters to lines
  file_size = file_size / line_size;
  // Chunk size needs to be adjusted to fit memory limitations and number of
  // threads
  short chunk_size = 2;
  short n_threads = file_size / chunk_size;
  rewind(file);

  float **points = create_matrixf(file_size, 3);
#pragma omp parallel
  {
    unsigned int starting_line = omp_get_thread_num() * chunk_size;
    // Line number being read will correspond to row index
    for (short i = starting_line; i < starting_line + chunk_size; i++) {
      fseek(file, i * line_size, SEEK_SET);
      for (short j = 0; j < 3; j++) {
        char buffer[point_char_size];
        // fgets reads n-1 chars
        fgets(buffer, point_char_size + 1, file);
        points[i][j] = atof(buffer);
      }
    }
  }
  print_matrix((const float **)points, file_size, 3);
  free_matrixf(points);
}
