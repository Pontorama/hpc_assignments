#include <stdio.h>
#include <stdlib.h>

void contiguous_memory() {
  int size = 10;
  int *asentries = (int *)malloc(sizeof(int) * size * size);
  int **as = (int **)malloc(sizeof(int *) * size);

  for (size_t i = 0, j = 0; i < size; ++i, j += size) {
    as[i] = asentries + j;
  }

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      as[i][j] = 0;
    }
  }

  printf("%d\n", as[0][0]);
  free(asentries);
  free(as);
}

void non_contiguous_mem() {
  int size = 10;

  int **as = (int **)malloc(sizeof(int *) * size);
  for (size_t i = 0; i < size; ++i) {
    as[i] = (int *)malloc(sizeof(int) * size);
  }

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      as[i][j] = 0;
    }
  }

  printf("%d\n", as[0][0]);

  for (size_t i = 0; i < size; ++i) {
    free(as[i]);
  }
  free(as);
}

int main() {
  contiguous_memory();
  non_contiguous_mem();
}
