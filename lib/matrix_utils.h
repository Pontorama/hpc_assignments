#ifndef Matrix_utils_h
#define Matrix_utils_h

#include <stdlib.h>

double **create_matrix(size_t n_rows, size_t n_columns);
float **create_matrixf(size_t n_rows, size_t n_columns);
void free_matrix(double **matrix);
void free_matrixf(float **matrix);

#endif
