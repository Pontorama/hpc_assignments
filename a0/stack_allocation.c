#include <stdio.h>

int main() {
  int size = 8192 / 4 * 1000;
  int as[size];
  for (size_t i = 0; i < size; ++i) {
    as[i] = 0;
  }
  // Segfault if size is greater than int size 2^32
  // or if size > than stack size (~8MB)
  printf("%d\n", as[0]);
}
