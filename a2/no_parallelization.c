/*
 * This is a "baseline" to the problem, mostly to familiarize myself with the
 * problem
 * */
#include "../lib/matrix_utils.h"
#include "read_text_from_file.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float get_distance(float *point_1, float *point_2) {
  // Assumes points in R^3
  return sqrtf((point_2[0] - point_1[0]) * (point_2[0] - point_1[0]) +
               (point_2[1] - point_1[1]) * (point_2[1] - point_1[1]) +
               (point_2[2] - point_1[2]) * (point_2[2] - point_1[2]));
}

int main() {
  // Assignment 2:
  // Calculate distance between points
  char filename[10] = "cells.txt";
  size_t n_points = 10;
  float **points = create_matrixf(n_points, 3);
  int result = read_points_from_txt_file(filename, points, n_points);
  if (result == -1) {
    printf("Failed to load file %s", filename);
    return -1;
  }

  float **distances = create_matrixf(n_points, n_points);
  // Maximum distance is 34.64 (Between -10 & 10)
  unsigned short max_distance = 3465;
  unsigned int *histogram =
      (unsigned int *)malloc(sizeof(unsigned int) * max_distance);
  unsigned short distance = 0;
  for (size_t i = 0; i < n_points; i++) {
    for (size_t j = i; j < n_points; j++) {
      // Get distance and round it to 2 decimals, multiply by 100
      // Rounding errors, should do the add .5 trick to properly round
      distance = (unsigned short)(100 * get_distance(points[i], points[j]));
      histogram[distance] += 1;
    }
  }

  for (size_t i = 0; i < max_distance; i++) {
    if (histogram[i] != 0)
      printf("%.2f: %i\n", (float)i / 100, histogram[i]);
  }
}
