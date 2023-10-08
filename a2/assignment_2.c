#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
static inline unsigned short get_distance(short *restrict point_1,
                                          short *restrict point_2) {

  // Assumes points in R^3
  short sub_one = point_2[0] - point_1[0];
  short sub_two = point_2[1] - point_1[1];
  short sub_three = point_2[2] - point_1[2];
  return sqrtf(sub_one * sub_one + sub_two * sub_two + sub_three * sub_three);
}

static inline short convert_to_short(char a[2], char b[3]) {
  return (short)(a[0] - '0') * 10000 + (short)(a[1] - '0') * 1000 +
         (short)(b[0] - '0') * 100 + (short)(b[1] - '0') * 10 +
         (short)(b[2] - '0');
}

int main(int argc, char *argv[]) {
  // Open file, get file length
  char filename[10] = "cells";
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    printf("Could not read file.\n");
    return -1;
  }
  fseek(file, 0L, SEEK_END);
  const unsigned short line_size = 24;
  const unsigned int file_size = ftell(file) / line_size;
  rewind(file);
  const unsigned short n_threads = atoi(argv[1] + 2);

  omp_set_num_threads(n_threads);

  const unsigned short max_distance = 3465;
  unsigned long long int *distance_entries = (unsigned long long int *)malloc(
      sizeof(unsigned long long int) * max_distance * n_threads);
  for (size_t i = 0; i < max_distance * n_threads; i++) {
    distance_entries[i] = 0;
  }

  unsigned long long int **distances = (unsigned long long int **)malloc(
      sizeof(unsigned long long int) * max_distance);
  for (size_t i = 0; i < n_threads; i++) {
    distances[i] = distance_entries + i * max_distance;
  }

  unsigned int max_chunk = 600000;
  unsigned int chunk_size;
  // Case 1: Whole file read
  if (file_size < max_chunk) {
    chunk_size = file_size;
  }

  // Case 2: Only two chunks
  else if (file_size * 2 < max_chunk) {
    chunk_size = max_chunk / 2;
  }

  // Case 3: Several chunks
  else {
    chunk_size = max_chunk / 3;
  }

  //	printf("Chunk size: %i\n", chunk_size);
  short *base_entries = (short *)malloc(sizeof(short) * chunk_size * 3);
  short **base = (short **)malloc(sizeof(short *) * chunk_size);
  // short *buffer_entries = (short *)aligned_alloc(sizeof(short)*3,
  // sizeof(short) * chunk_size * 3); short **buffer = (short
  // **)malloc(sizeof(short *) * chunk_size);
  short *moving_entries = (short *)malloc(sizeof(short) * chunk_size * 3);
  short **moving = (short **)malloc(sizeof(short *) * chunk_size);

  for (unsigned int i = 0; i < chunk_size; i++) {
    base[i] = base_entries + i * 3;
    // buffer[i] = buffer_entries + i * 3;
    moving[i] = moving_entries + i * 3;
  }

  unsigned int current_cell_point = 0;
  unsigned int file_point = 0;
  // short *current_cell = (short *)aligned_alloc(sizeof(short) * 3,
  // sizeof(short)*3);

  // float *compare_cell = (float*)malloc(sizeof(float)*3);
  char *chunk_string =
      (char *)malloc(sizeof(char) * (chunk_size * line_size + 1));
  char *moving_string =
      (char *)malloc(sizeof(char) * (chunk_size * line_size + 1));
  unsigned int compare;

  char *format = (char *)malloc(sizeof(char) * 100);
  sprintf(format, "%%%ic", chunk_size * line_size);

  short distance_one;
  short distance_two;
  short distance_three;
  short distance_four;

#pragma omp parallel shared(                                                   \
        base, base_entries, distances, chunk_string, current_cell_point,       \
            file_point, format) private(distance_one, distance_two,            \
                                            distance_three, distance_four)
  {
    for (unsigned int base_iter = 0; base_iter < file_size / chunk_size;
         base_iter++) {
      // chunk_size = file_size;
      fseek(file, base_iter * chunk_size, SEEK_SET);
      fscanf(file, format, chunk_string);
      for (unsigned int moving_iter = base_iter;
           moving_iter < file_size / chunk_size; moving_iter++) {
        if (base_iter == moving_iter) {
          fseek(file, -chunk_size * line_size, SEEK_CUR);
        }
        fscanf(file, format, moving_string);
#pragma omp for
        for (unsigned int i = 0; i < chunk_size; i++) {
          // printf("%i\n", i);
          char buffer[18];
          // Could be optimized by not using sscanf
          sscanf(chunk_string + i * line_size,
                 "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", buffer, buffer + 1,
                 buffer + 3, buffer + 6, buffer + 7, buffer + 9, buffer + 12,
                 buffer + 13, buffer + 15);
          base[i][0] = convert_to_short(buffer + 1, buffer + 3) *
                       (buffer[0] == '+' ? 1 : -1);
          base[i][1] = convert_to_short(buffer + 7, buffer + 9) *
                       (buffer[6] == '+' ? 1 : -1);
          base[i][2] = convert_to_short(buffer + 13, buffer + 15) *
                       (buffer[12] == '+' ? 1 : -1);
          sscanf(moving_string + i * line_size,
                 "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", buffer, buffer + 1,
                 buffer + 3, buffer + 6, buffer + 7, buffer + 9, buffer + 12,
                 buffer + 13, buffer + 15);
          moving[i][0] = convert_to_short(buffer + 1, buffer + 3) *
                         (buffer[0] == '+' ? 1 : -1);
          moving[i][1] = convert_to_short(buffer + 7, buffer + 9) *
                         (buffer[6] == '+' ? 1 : -1);
          moving[i][2] = convert_to_short(buffer + 13, buffer + 15) *
                         (buffer[12] == '+' ? 1 : -1);
        }
        short thread;
        for (unsigned int i = 0; i < chunk_size; i++) {
#pragma omp for
          for (unsigned int j = i + 1; j < chunk_size - 3; j += 4) {
            // printf("base[i]: %i %i %i\n", base[i][0], base[i][1],
            // base[i][2]); printf("base[i+j]: %i %i %i\n", base[i+j][0],
            // base[i+j][1], base[i+j][2]);
            distance_one = get_distance(base[i], moving[j]) / 10;
            distance_two = get_distance(base[i], moving[j + 1]) / 10;
            distance_three = get_distance(base[i], moving[j + 2]) / 10;
            distance_four = get_distance(base[i], moving[j + 3]) / 10;

            thread = omp_get_thread_num();

            distances[thread][distance_one] += 1;
            distances[thread][distance_two] += 1;
            distances[thread][distance_three] += 1;
            distances[thread][distance_four] += 1;

            // #pragma omp atomic
          }
        }
      }

      unsigned int rest_size = file_size % chunk_size;
      sprintf(format, "%%%ic", rest_size * line_size);
      fscanf(file, format, chunk_string);
      unsigned int rest_size_quad = rest_size / 4;
      for (unsigned int rest_iter = 0; rest_iter < rest_size; rest_iter++) {
        char buffer[18];
        // Could be optimized by not using sscanf
        sscanf(chunk_string + rest_iter * line_size,
               "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", buffer, buffer + 1,
               buffer + 3, buffer + 6, buffer + 7, buffer + 9, buffer + 12,
               buffer + 13, buffer + 15);
        base[rest_iter][0] = convert_to_short(buffer + 1, buffer + 3) *
                             (buffer[0] == '+' ? 1 : -1);
        base[rest_iter][1] = convert_to_short(buffer + 7, buffer + 9) *
                             (buffer[6] == '+' ? 1 : -1);
        base[rest_iter][2] = convert_to_short(buffer + 13, buffer + 15) *
                             (buffer[12] == '+' ? 1 : -1);
      }
      unsigned short thread;
      for (unsigned int i = 0; i < rest_size_quad; i += 4) {
#pragma omp for
        for (unsigned int j = i + 1; j < rest_size_quad - 3; j += 4) {
          // printf("base[i]: %i %i %i\n", base[i][0], base[i][1],
          // base[i][2]); printf("base[i+j]: %i %i %i\n", base[i+j][0],
          // base[i+j][1], base[i+j][2]);
          distance_one = get_distance(base[i], base[j]) / 10;
          distance_two = get_distance(base[i], base[j + 1]) / 10;
          distance_three = get_distance(base[i], base[j + 2]) / 10;
          distance_four = get_distance(base[i], base[j + 3]) / 10;

          thread = omp_get_thread_num();

          distances[thread][distance_one] += 1;
          distances[thread][distance_two] += 1;
          distances[thread][distance_three] += 1;
          distances[thread][distance_four] += 1;
        }
      }
      for (unsigned int i = rest_size_quad; i < rest_size; i++) {
        for (unsigned int j = i + 1; j < rest_size; j++) {
          distance_one = get_distance(base[i], base[j]) / 10;
          distances[0][distance_one] += 1;
        }
      }
    }
  }
  free(base_entries);
  free(base);
  free(moving_entries);
  free(moving);
  free(format);
  free(chunk_string);
  free(moving_string);

  unsigned long long int *result =
      (unsigned long long *)malloc(sizeof(unsigned long long) * max_distance);
  for (size_t i = 0; i < max_distance; i++) {
    result[i] = 0;
  }
  for (size_t j = 0; j < n_threads; j++) {
    for (size_t i = 0; i < max_distance; i++) {
      result[i] += distances[j][i];
    }
  }

  for (size_t i = 0; i < max_distance; i++) {
    if (result[i] != 0)
      printf("%05.2f %llu\n", i / 100.0, result[i]);
    // printf("%i : %i (%i)\n", i, i*i, i*i - (i-1)*(i-1));
  }
  free(distance_entries);
  free(distances);
  free(result);
}
