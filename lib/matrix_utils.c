#include <stdlib.h>

double **create_matrix(size_t n_rows, size_t n_columns) {
  // Creates a (row major) matrix of size n_rows * n_columns with all elements
  // set to 0 max_val controls how large the random elements can be (-max_val,
  // max_val)
  double **as = (double **)malloc(sizeof(double *) * n_rows);
  double *asentries = (double *)malloc(sizeof(double) * n_rows * n_columns);

  for (size_t i = 0; i < n_rows; ++i) {
    as[i] = asentries + i * n_columns;
  }

  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_columns; ++j) {
      as[i][j] = 0;
    }
  }

  return as;
}

void free_matrix(double **matrix) {
  // Free memory allocated for matrix by create_matrix function
  free(matrix[0]);
  free(matrix);
}
