#include "ppm_utils.h"

char **row_str;

int write_header(char *file_to_write, unsigned int n_lines) {
  file = fopen(file_to_write, "w");
  N_LINES = n_lines;
  if (file == NULL) {
    printf("Failed to open file %s\n", file_to_write);
    return -1;
  }
  // Set row_str, only need to allocate once
  char *row_str_entries = (char *)malloc(4 * n_lines * sizeof(char));
  for (size_t i = 0; i < n_lines; ++i) {
    row_str[i] = row_str_entries + (4 * i);
  }

  fprintf(file, "P3\n");
  fprintf(file, "%i %i\n", n_lines, n_lines);
  fprintf(file, "%i\n", MAX_COLOR_VAL);

  return 0;
}

void close_ppm() {
  fclose(file);
  free(row_str[0]);
  free(row_str);
}

int write_row(unsigned int *values) {
  char tmp[4];
  for (size_t i = 0; i < N_LINES - 1; ++i) {
    sprintf(tmp, "%i ", values[i]);
    row_str[i] = tmp;
  }
  // Write with \n
  sprintf(tmp, "%i\n", values[N_LINES - 1]);
  row_str[N_LINES - 1] = tmp;
  fwrite(row_str, N_LINES * 4 * sizeof(char), 1, file);
  return 0;
}
