#include "../lib/matrix_utils.h"
#include "read_text_from_file.h"
#include <stdio.h>

int main() {
  size_t n_cols = 10;
  char filename[10] = "cells.txt";
  float **matrix = create_matrixf(n_cols, 3);

  int read_ok = read_points_from_txt_file(filename, matrix, n_cols);
  if (read_ok == -1) {
    printf("Read not possible.\n");
    return -1;
  }

  print_matrix((const float **)matrix, n_cols, 3);

  FILE *file = fopen("cells.txt", "r");
  char buffer[7];
  fgets(buffer, 7, file);
  float num = atof(buffer);
  printf("\n%f\n", num);
  free(matrix[0]);
  free(matrix);
}
