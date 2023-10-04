#ifndef read_from_txt_file_h
#define read_from_txt_file_h
#include <stdlib.h>

int read_points_from_txt_file(char *file_name, float **output_matrix,
                              size_t n_rows);

void print_matrix(const float **matrix, size_t n_rows, size_t n_columns);
#endif
