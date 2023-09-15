#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  int *as;
  as = (int *)malloc(10 * sizeof(int));
  int sum = 0;
  /*
  for (int i = 0; i < 10; ++i) {
    as[i] = i;
  }
  */

  for (int i = 0; i < 10; ++i) {
    sum += as[i];
  }

  printf("%d\n", sum);
  free(as);
}
