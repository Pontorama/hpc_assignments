#include "newton.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

unsigned int n_threads;
unsigned int n_lines;
const float conv_lim = 10e-3;
const unsigned int max_iter = 128;
const float max_num_size =
    10e10; // TODO: this will cause overflow? Number too big prob.

int main(int argc, char **argv) {
  unsigned int degree = 1; // 2 as standard for now
  int opt;
  while (optind < argc) {
    if ((opt = getopt(argc, argv, "t:l:")) != -1) {
      // Handle optional arguments
      switch (opt) {
      case 't':
        n_threads = (unsigned int)atoi(optarg);
        break;
      case 'l':
        n_lines = (unsigned int)atoi(optarg);
      default:
        abort();
        break;
      }
    } else {
      // Arguments without flags (mandatory args)
      degree = (unsigned int)atoi(argv[optind]);
    }
  }

  // Set function pointer(s) to correct degree
  complex_float (*f)(complex_float *);
  complex_float (*df)(complex_float *);
  switch (degree) {
  case 2:
    f = &poly_2;
    df = &dpoly_2;
    break;
  case 3:
    f = &poly_3;
    df = &dpoly_3;
    break;
  case 4:
    f = &poly_4;
    df = &dpoly_4;
    break;
  case 5:
    f = &poly_5;
    df = &dpoly_5;
    break;
  case 6:
    f = &poly_6;
    df = dpoly_6;
    break;
  case 7:
    f = &poly_7;
    df = &dpoly_7;
    break;
  case 8:
    f = &poly_8;
    df = &dpoly_8;
    break;
  case 9:
    f = &poly_9;
    df = &dpoly_9;
    break;
  default:
    break;
  }

  // Newton iterations
}
