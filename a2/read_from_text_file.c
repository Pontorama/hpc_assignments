#include "../lib/matrix_utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int read_points_from_txt_file(char *file_name, float **output_matrix,
                              size_t n_rows) {
  // Assumes output_matrix has 3 columns
  FILE *file = fopen(file_name, "r");
  if (file == NULL) {
    printf("Unable to open file\n");
    return -1;
  }

  // The format for assignment 2 is 3 digits per row (R^3)
  // ±dd.ddd, space separated

  char c;
  size_t iterations = 0;
  const size_t n_digits = 5; // 5 as specified in assignment
  do {
    c = fgetc(file);
    // First char is sign
    int sign = 1;
    if ('-' == c) {
      sign = -1;
    }
    double number = 0;
    for (size_t i = 0; i < n_digits; i++) {
      c = fgetc(file);
      if (c != '.') {
        number += sign * (10.0 / pow(10, i)) * atoi(&c);
      } else {
        i--;
        continue;
      }
    }
    // Read the separating space char but don't use it
    c = fgetc(file);
    output_matrix[iterations / 3][iterations % 3] = number;
    iterations++;
  } while (c != EOF && iterations / 3 < n_rows);
  return 0;
}

void print_matrix(const float **matrix, size_t n_rows, size_t n_columns) {
  for (size_t i = 0; i < n_rows; i++) {
    printf("------------------------\n");
    printf("| ");
    for (size_t j = 0; j < n_columns; j++) {
      printf(" %.2f |", matrix[i][j]);
    }
    printf("\n");
  }
}
