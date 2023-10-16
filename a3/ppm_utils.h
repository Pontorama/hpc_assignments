#ifndef ppm_utils_h
#define ppm_utils_h
#include "newton.h"
#include <stdio.h>

const unsigned int MAX_COLOR_VAL = 128;
unsigned int N_LINES;
FILE *file;

int write_header(char *file_to_write, unsigned int n_lines);
int write_row(unsigned int *values);
void close_ppm();

#endif
