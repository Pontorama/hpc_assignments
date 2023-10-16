#include "ppm_utils.h"

int write_header(char *file_to_write, unsigned int n_lines) {
  file = fopen(file_to_write, "w");
  N_LINES = n_lines;
  if (file == NULL) {
    printf("Failed to open file %s\n", file_to_write);
    return -1;
  }

  fprintf(file, "P3\n");
  fprintf(file, "%i %i\n", n_lines, n_lines);
  fprintf(file, "%i\n", MAX_COLOR_VAL);

  return 0;
}

void close_ppm() { fclose(file); }

int write_row(unsigned int *values) {
  char str_to_write[4 * N_LINES + 1] return 0;
}
