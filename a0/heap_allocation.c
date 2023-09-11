#include <stdio.h>
#include <stdlib.h>

int main() {
  int size = 10;
  int *as = (int *)malloc(sizeof(int) * size);

  for (size_t i = 0; i < size; ++i) {
    as[i] = 0;
  }

  printf("%d\n", as[0]);
  // Will not produce segfault for large arrays
  free(as);
}
