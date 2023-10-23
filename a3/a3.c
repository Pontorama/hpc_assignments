#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>
#include <unistd.h>

#define MAX_ITERS 128
#define CONV_LIMIT 0.003
#define MAX_NUM 10000000000
#define line_size 12
#define max_degree 10
#define PI 3.14159
#define x_max 2
#define x_min -2

unsigned int n_threads;
unsigned int n_lines;
unsigned short degree;

void create_color_palette(char *color_palette) {
  memcpy(color_palette, "158   1  66\n", line_size);
  memcpy(color_palette + line_size, "213  62  79\n", line_size);
  memcpy(color_palette + 2 * line_size, "244 109  67\n", line_size);
  memcpy(color_palette + 3 * line_size, "253 174  97\n", line_size);
  memcpy(color_palette + 4 * line_size, "254 224 139\n", line_size);
  memcpy(color_palette + 5 * line_size, "230 245 152\n", line_size);
  memcpy(color_palette + 6 * line_size, "171 221 164\n", line_size);
  memcpy(color_palette + 7 * line_size, "102 194 165\n", line_size);
  memcpy(color_palette + 8 * line_size, " 50 136 189\n", line_size);
  memcpy(color_palette + 9 * line_size, " 94  79 162\n", line_size);
  memcpy(color_palette + 10 * line_size, "  0   0   0\n", line_size);
}
void create_grey_palette(char *grey_palette) {
  for (size_t i = 0; i < MAX_ITERS; i++) {
    sprintf(grey_palette + line_size * i, "%3i %3i %3i\n", i, i, i);
  }
}

static inline int newton_iterations_unrolled4(
    double complex z1, double complex z2, double complex z3, double complex z4,
    double complex *root_list, unsigned int row_number, unsigned int col_number,
    unsigned char **attractors, unsigned char **convergences, mtx_t *mtx) {
  unsigned char root;
  double check_x1;
  double check_x2;
  double check_x3;
  double check_x4;
  double check_y1;
  double check_y2;
  double check_y3;
  double check_y4;
  double complex check1;
  double complex check2;
  double complex check3;
  double complex check4;
  double norm;
  double complex conj_z1;
  double complex conj_z2;
  double complex conj_z3;
  double complex conj_z4;
  double complex z_k1;
  double complex z_k2;
  double complex z_k3;
  double complex z_k4;

  double complex final_z1;
  double complex final_z2;
  double complex final_z3;
  double complex final_z4;

  bool done_z1 = 0;
  bool done_z2 = 0;
  bool done_z3 = 0;
  bool done_z4 = 0;
  size_t i;
  if (degree == 1) {
    z1 = 1;
    z2 = 1;
    z3 = 1;
    z4 = 1;
  } else {
    for (i = 0; i < MAX_ITERS; i++) {
      for (size_t j = 0; j < degree; j++) {
        double complex root = root_list[j];
        check1 = (root - z1);
        check2 = (root - z2);
        check3 = (root - z3);
        check4 = (root - z4);
        check_x1 = creal(check1);
        check_x2 = creal(check2);
        check_x3 = creal(check3);
        check_x4 = creal(check4);
        check_y1 = cimag(check1);
        check_y2 = cimag(check2);
        check_y3 = cimag(check3);
        check_y4 = cimag(check4);
        if (1000 * (check_x1 * check_x1 + check_y1 * check_y1) < 0.001 &&
            !done_z1) {
          attractors[row_number][col_number] = j;
          convergences[row_number][col_number] = i;
          final_z1 = z1;
          done_z1 = 1;
        } else if ((creal(z1) > 10000000000 || cimag(z1) > 10000000000 ||
                    creal(z1) < -10000000000 || cimag(z1) < -10000000000 ||
                    creal(z1 * conj(z1)) <= 0.001) &&
                   !done_z1) {
          attractors[row_number][col_number] = 10;
          convergences[row_number][col_number] = MAX_ITERS - 1;

          final_z1 = z1;
          done_z1 = 1;
        }
        if (1000 * (check_x2 * check_x2 + check_y2 * check_y2) < 0.001 &&
            !done_z2) {
          attractors[row_number][col_number + 1] = j;
          convergences[row_number][col_number + 1] = i;

          final_z2 = z2;
          done_z2 = 1;
        } else if ((creal(z2) > 10000000000 || cimag(z2) > 10000000000 ||
                    creal(z2) < -10000000000 || cimag(z2) < -10000000000 ||
                    creal(z2 * conj(z2)) <= 0.001) &&
                   !done_z2) {
          attractors[row_number][col_number + 1] = 10;
          convergences[row_number][col_number + 1] = MAX_ITERS - 1;

          final_z2 = z2;
          done_z2 = 1;
        }
        if (1000 * (check_x3 * check_x3 + check_y3 * check_y3) < 0.001 &&
            !done_z3) {
          attractors[row_number][col_number + 2] = j;
          convergences[row_number][col_number + 2] = i;
          final_z3 = z3;
          done_z3 = 1;
        } else if ((creal(z3) > 10000000000 || cimag(z3) > 10000000000 ||
                    creal(z3) < -10000000000 || cimag(z3) < -10000000000 ||
                    creal(z3 * conj(z3)) <= 0.001) &&
                   !done_z3) {
          attractors[row_number][col_number + 2] = 10;
          convergences[row_number][col_number + 2] = MAX_ITERS - 1;

          final_z3 = z3;
          done_z3 = 1;
        }
        if (1000 * (check_x4 * check_x4 + check_y4 * check_y4) < 0.001 &&
            !done_z4) {
          attractors[row_number][col_number + 3] = j;
          convergences[row_number][col_number + 3] = i;
          final_z4 = z4;
          done_z4 = 1;
        } else if ((creal(z4) > 10000000000 || cimag(z4) > 10000000000 ||
                    creal(z4) < -10000000000 || cimag(z4) < -10000000000 ||
                    creal(z4 * conj(z4)) <= 0.001) &&
                   !done_z4) {
          attractors[row_number][col_number + 3] = 10;
          convergences[row_number][col_number + 3] = MAX_ITERS - 1;

          final_z4 = z4;
          done_z4 = 1;
        }
      }
      if (done_z1 && done_z2 && done_z3 && done_z4) {
        return 0;
      }
      z_k1 = z1;
      z_k2 = z2;
      z_k3 = z3;
      z_k4 = z4;
      for (unsigned short k = 0; k < degree - 2; k++) {
        z_k1 *= z1;
        z_k2 *= z2;
        z_k3 *= z3;
        z_k4 *= z4;
      }
      conj_z1 = conj(z_k1);
      conj_z2 = conj(z_k2);
      conj_z3 = conj(z_k3);
      conj_z4 = conj(z_k4);
      z1 = z1 - (z_k1 * z1 - 1) * conj_z1 / (degree * creal(z_k1 * conj_z1));
      z2 = z2 - (z_k2 * z2 - 1) * conj_z2 / (degree * creal(z_k2 * conj_z2));
      z3 = z3 - (z_k3 * z3 - 1) * conj_z3 / (degree * creal(z_k3 * conj_z3));
      z4 = z4 - (z_k4 * z4 - 1) * conj_z4 / (degree * creal(z_k4 * conj_z4));

      if (done_z1) {
        z1 = final_z1;
      }
      if (done_z2) {
        z2 = final_z2;
      }
      if (done_z3) {
        z3 = final_z3;
      }
      if (done_z4) {
        z4 = final_z4;
      }
    }
  }

  if (i == MAX_ITERS) {
    if (!done_z1) {
      attractors[row_number][col_number] = 10;
      convergences[row_number][col_number] = MAX_ITERS - 1;
    }
    if (!done_z2) {
      attractors[row_number][col_number + 1] = 10;
      convergences[row_number][col_number + 1] = MAX_ITERS - 1;
    }
    if (!done_z3) {
      attractors[row_number][col_number + 2] = 10;
      convergences[row_number][col_number + 2] = MAX_ITERS - 1;
    }
    if (!done_z4) {
      attractors[row_number][col_number + 3] = 10;
      convergences[row_number][col_number + 3] = MAX_ITERS - 1;
    }
  }

  return 0;
}
static inline int newton_iterations(double complex z, double complex *root_list,
                                    unsigned int row_number,
                                    unsigned int col_number,
                                    unsigned char **attractors,
                                    unsigned char **convergences, mtx_t *mtx) {
  double check_x;
  double check_y;
  double complex check;
  double norm;
  double complex conj_z;
  double complex z_k;
  size_t i;
  if (degree == 1) {
    z = 1;
    mtx_lock(mtx);
    attractors[row_number][col_number] = 0;
    convergences[row_number][col_number] = 1;
    mtx_unlock(mtx);
    return 0;
  } else {
    for (i = 0; i < MAX_ITERS; i++) {
      for (size_t j = 0; j < degree; j++) {
        check = (root_list[j] - z);
        check_x = creal(check);
        check_y = cimag(check);
        if (1000 * (check_x * check_x + check_y * check_y) < CONV_LIMIT) {
          attractors[row_number][col_number] = j;
          convergences[row_number][col_number] = MAX_ITERS - 1;
          return 0;
        }
      }
      if (creal(z) > 10000000000 || cimag(z) > 10000000000 ||
          creal(z) < -10000000000 || cimag(z) < -10000000000 ||
          creal(z * conj(z)) <= 0.001) {
        attractors[row_number][col_number] = 10;
        convergences[row_number][col_number] = 2 * i - 1;
        return 0;
      }

      z_k = z;
      for (unsigned short k = 0; k < degree - 2; k++) {
        z_k *= z;
      }
      conj_z = conj(z_k);
      z = z - (z_k * z - 1) * conj_z / (degree * creal(z_k * conj_z));
    }
  }

  if (i == MAX_ITERS) {
    attractors[row_number][col_number] = 10;
    convergences[row_number][col_number] = MAX_ITERS - 1;
  }

  return 0;
}

typedef struct {
  mtx_t *mtx;
  cnd_t *cnd;
  unsigned int thread_index;
  int *last_processed_row; // Used to keep track of what rows are ready
                           // to be written
  unsigned int n_lines;
  unsigned int n_threads;
  unsigned char **attractors;
  unsigned char **convergence;
  double complex *root_list;
  double dx;
} thrd_info_compute_t;

typedef struct {
  mtx_t *mtx;
  cnd_t *cnd;
  int *last_processed_row;
  unsigned int n_threads;
  unsigned char **attractors;
  unsigned char **convergence;
  char *color_palette;
  char *greyscale_palette;
  FILE *attractors_file;
  FILE *convergence_file;
} thrd_info_write_t;

int main_compute_thread(void *args) {
  // Unpack arguments from struct
  thrd_info_compute_t *comp_args = (thrd_info_compute_t *)args;
  unsigned int n_lines = comp_args->n_lines;
  unsigned int n_threads = comp_args->n_threads;
  unsigned int thread_index = comp_args->thread_index;
  int *last_processed_row = comp_args->last_processed_row;
  unsigned char **attractors = comp_args->attractors;
  unsigned char **convergences = comp_args->convergence;
  double dx = comp_args->dx;
  double complex *root_list = comp_args->root_list;
  mtx_t *mtx = comp_args->mtx;
  cnd_t *cnd = comp_args->cnd;

  double complex z1;
  double complex z2;
  double complex z3;
  double complex z4;
  int r;
  for (unsigned int i = thread_index; i < n_lines; i += n_threads) {
    for (unsigned int j = 0; j < n_lines; j += 4) {
      z1 = x_min + dx * j + I * (x_max - i * dx);
      z2 = x_min + dx * (j + 1) + I * (x_max - i * dx);
      z3 = x_min + dx * (j + 2) + I * (x_max - i * dx);
      z4 = x_min + dx * (j + 3) + I * (x_max - i * dx);
      r = newton_iterations_unrolled4(z1, z2, z3, z4, root_list, i, j,
                                      attractors, convergences, mtx);
    }
    for (unsigned int j = n_lines / 4; j < n_lines; j++) {
      z1 = x_min + dx * j + I * (x_max - i * dx);
      r = newton_iterations(z1, root_list, i, j, attractors, convergences, mtx);
    }
    mtx_lock(mtx);
    // Interweaved processing
    last_processed_row[thread_index] = i + n_threads;
    mtx_unlock(mtx);
    // Signal
    cnd_signal(cnd);
  }
  return 0;
}

int main_write_thread(void *args) {
  // Unpack input args from struct
  thrd_info_write_t *write_args = (thrd_info_write_t *)args;
  mtx_t *mtx = write_args->mtx;
  cnd_t *cnd = write_args->cnd;
  int *last_processed_row = write_args->last_processed_row;
  unsigned int n_threads = write_args->n_threads;
  unsigned char **attractors = write_args->attractors;
  unsigned char **convergences = write_args->convergence;
  char *color_palette = write_args->color_palette;
  char *greyscale_palette = write_args->greyscale_palette;
  FILE *attractors_file = write_args->attractors_file;
  FILE *convergence_file = write_args->convergence_file;

  unsigned int written_lines = 0;
  unsigned int block_size;
  char *buffer_attractors;
  char *buffer_convergences;
  unsigned short root_index;
  unsigned short greyscale_index;
  unsigned int lowest_index = n_lines;
  while (written_lines < n_lines) {
    for (mtx_lock(mtx);;) {
      lowest_index = n_lines;
      // Find lowest line index not yet processed
      for (unsigned int t = 0; t < n_threads; t++) {
        if (last_processed_row[t] < lowest_index) {
          lowest_index = last_processed_row[t];
        }
      }

      if (lowest_index <= written_lines) {
        cnd_wait(cnd, mtx);
      } else {
        mtx_unlock(mtx);
        break;
      }
    }
    block_size = (lowest_index - written_lines);
    buffer_attractors = (char *)malloc(block_size * line_size * n_lines);
    buffer_convergences = (char *)malloc(block_size * line_size * n_lines);

    for (unsigned int i = 0; i < block_size; i++) {
      // Write processed line with lowest index
      for (unsigned int j = 0; j < n_lines; j++) {
        root_index = attractors[i + written_lines][j];
        greyscale_index = convergences[i + written_lines][j];
        memcpy(buffer_attractors + (i * n_lines * line_size) + j * line_size,
               color_palette + root_index * line_size, line_size);
        memcpy(buffer_convergences + (i * n_lines * line_size) + j * line_size,
               greyscale_palette + greyscale_index * line_size, line_size);
      }
    }

    written_lines += block_size;
    fwrite(buffer_attractors, sizeof(char), block_size * line_size * n_lines,
           attractors_file);
    fwrite(buffer_convergences, sizeof(char), block_size * line_size * n_lines,
           convergence_file);
    free(buffer_attractors);
    free(buffer_convergences);
  }
  return 0;
}

int main(int argc, char **argv) {
  degree = 1;
  int opt;
  while (optind < argc - 1) {
    if ((opt = getopt(argc, argv, "t:l:")) != -1) {
      // Handle optional arguments
      switch (opt) {
      case 't':
        n_threads = (unsigned int)atoi(optarg);
        break;
      case 'l':
        n_lines = (unsigned int)atoi(optarg);
        break;
      default:
        abort();
        break;
      }
    }
  }
  degree = (unsigned int)atoi(argv[argc - 1]);
  // Create matricies
  unsigned char **attractors =
      (unsigned char **)malloc(sizeof(unsigned char *) * n_lines);
  unsigned char *attractors_entries =
      (unsigned char *)malloc(sizeof(unsigned char) * n_lines * n_lines);
  unsigned char **convergences =
      (unsigned char **)malloc(sizeof(unsigned char *) * n_lines);
  unsigned char *convergences_entries =
      (unsigned char *)malloc(sizeof(unsigned char) * n_lines * n_lines);
  for (int i = 0; i < n_lines; i++) {
    attractors[i] = attractors_entries + i * n_lines;
    convergences[i] = convergences_entries + i * n_lines;
  }

  // Calculate roots
  double complex *root_list =
      (double complex *)malloc(sizeof(double complex) * max_degree);
  for (size_t i = 0; i < degree; i++) {
    root_list[i] =
        cos(2 * PI * i / (double)degree) + I * sin(2 * PI * i / (double)degree);
  }
  // Create color palette
  char *color_palette = (char *)malloc(line_size * 11 + 1); // 11 colors
  create_color_palette(color_palette);
  char *greyscale_palette =
      (char *)malloc(sizeof(char) * MAX_ITERS * line_size + 1);
  create_grey_palette(greyscale_palette);

  // Create compute threads
  mtx_t mtx;
  cnd_t cnd;
  mtx_init(&mtx, mtx_plain);
  cnd_init(&cnd);

  const double dx = (x_max - x_min) / (1. * n_lines + 1);

  thrd_t threads[n_threads];
  thrd_info_compute_t compute_args[n_threads];
  // last_processed_row will be used to tell the write thread what rows are
  // ready to write
  int *last_processed_row = (int *)malloc(sizeof(int) * n_threads);
  for (unsigned int t = 0; t < n_threads; ++t) {
    last_processed_row[t] = t; // negative num means no rows ready yet
    compute_args[t].last_processed_row = last_processed_row;
    compute_args[t].n_lines = n_lines;
    compute_args[t].n_threads = n_threads;
    compute_args[t].thread_index = t;
    compute_args[t].mtx = &mtx;
    compute_args[t].cnd = &cnd;
    compute_args[t].attractors = attractors;
    compute_args[t].convergence = convergences;
    compute_args[t].dx = dx;
    compute_args[t].root_list = root_list;

    int r = thrd_create(threads + t, main_compute_thread,
                        (void *)(compute_args + t));
    if (r != thrd_success) {
      printf("Failed to create thread %i\n", t);
      exit(1);
    }
  }

  FILE *attractors_file;
  FILE *convergence_file;

  char *a_file_name = (char *)malloc(sizeof(char) * 100);
  char *c_file_name = (char *)malloc(sizeof(char) * 100);
  sprintf(a_file_name, "newton_attractors_x%d.ppm", degree);
  sprintf(c_file_name, "newton_convergence_x%d.ppm", degree);

  attractors_file = fopen(a_file_name, "w");
  convergence_file = fopen(c_file_name, "w");
  // Write header
  fprintf(attractors_file, "P3\n%d %d\n255\n", n_lines, n_lines);
  fprintf(convergence_file, "P3\n%d %d\n255\n", n_lines, n_lines);

  thrd_info_write_t write_args;
  thrd_t write_thread;
  write_args.mtx = &mtx;
  write_args.cnd = &cnd;
  write_args.last_processed_row = last_processed_row;
  write_args.n_threads = n_threads;
  write_args.attractors = attractors;
  write_args.convergence = convergences;
  write_args.attractors_file = attractors_file;
  write_args.convergence_file = convergence_file;
  write_args.color_palette = color_palette;
  write_args.greyscale_palette = greyscale_palette;

  int r = thrd_create(&write_thread, main_write_thread, (void *)(&write_args));
  if (r != thrd_success) {
    printf("Unable to create write thread.\n");
    exit(1);
  }

  // Join threads
  for (unsigned int i = 0; i < n_threads; i++) {
    int r;
    thrd_join(threads[i], &r);
  }
  {
    int r;
    thrd_join(write_thread, &r);
  }
  // Free variables
  free(convergences_entries);
  free(attractors_entries);
  free(convergences);
  free(attractors);
  free(color_palette);
  free(greyscale_palette);
  free(root_list);
  fclose(attractors_file);
  fclose(convergence_file);
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
}
