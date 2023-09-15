#include <stdlib.h>

int main() {
  int *as = NULL;
  int sum = 0;

  for (int i = 0; i < 10; ++i) {
    as[i] = i;
  }

  for (int i = 0; i < 10; ++i) {
    sum += as[i];
  }
  free(as);
}
