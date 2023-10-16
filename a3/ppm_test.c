#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

unsigned int n_threads;
unsigned int n_lines;

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
        break;
      }
    } else {
      // Arguments without flags (mandatory args)
      degree = (unsigned int)atoi(argv[argc - 1]);
      break;
    }
  }
}
