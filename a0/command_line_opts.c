#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char **argv) {
  // parse command line opts
  int a;
  int b;

  int opt;
  printf("Parsing\n");
  while ((opt = getopt(argc, argv, "a:b:")) != -1) {
    switch (opt) {
    case 'a':
      a = (int)(*optarg - '0');
      break;
    case 'b':
      b = (int)(*optarg - '0');
      break;
    default:
      abort();
    }
  }

  printf("A is %i, B is %i\n", a, b);
}
