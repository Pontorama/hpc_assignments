#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
static inline 
unsigned short get_distance(short* restrict point_1, short* restrict point_2) {
  
// Assumes points in R^3
	short sub_one = point_2[0] - point_1[0];
	short sub_two = point_2[1] - point_1[1];
	short sub_three = point_2[2] - point_1[2];
  return sqrtf(sub_one*sub_one + sub_two*sub_two + sub_three*sub_three);
//simd simdlen(3) aligned(current_cell, all_cells:3)
/*
	__m128i ar = _mm_loadu_si128((__m128i*)(point_1));
	__m128i br = _mm_loadu_si128((__m128i*)(point_2));

	__m128i cr = _mm_sub_epi64(ar,br);
	_mm_storeu_si128((__m128i*)(tmp), cr);

	return sqrtf(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
*/
}

static inline short convert_to_short(char a[2], char b[3]){
	return (short)(a[0] - '0')*10000+(short)(a[1] - '0')*1000+(short)(b[0] - '0')*100
			+ (short)(b[1] - '0')*10 + (short)(b[2] - '0');
}


int main(int argc, char *argv[]) {
  // Read cmd line inputs here

  // Open file, get file length
  char filename[10] = "cells";
  FILE *file = fopen(filename, "r");
  if (file == NULL){
    printf("Could not read file.\n");
    return -1;
  }
  fseek(file, 0L, SEEK_END);
  const unsigned short line_size = 24;
  const unsigned int file_size = ftell(file) / line_size;
  rewind(file);
  const unsigned short n_threads = atoi(argv[1]+2);
//	printf("n_threads = %i\n", n_threads);
  omp_set_num_threads(n_threads);
/*
  short *all_cell_entries = (short *)aligned_alloc(sizeof(short)*3, sizeof(short) * chunk_size * 3);
  short **all_cells = (short **)malloc(sizeof(short *) * chunk_size);
  for (short i = 0; i < chunk_size; i++) {
    all_cells[i] = all_cell_entries + i * 3;
  }

  for (size_t i = 0; i < chunk_size; i++) {
    for (size_t j = 0; j < 3; j++) {
      all_cells[i][j] = 0;
    }
  }
*/
  const unsigned short max_distance = 3465;
/*
  unsigned long long int *distances = (unsigned long long int *)malloc(
      sizeof(unsigned long long int) * max_distance);
	for (size_t i = 0; i < max_distance; i++){
		distances[i] = 0;
	}
*/
  unsigned long long int *distance_entries = (unsigned long long int *)malloc(
      sizeof(unsigned long long int) * max_distance*n_threads);
	for (size_t i = 0; i < max_distance*n_threads; i++){
		distance_entries[i] = 0;
	}

	unsigned long long int **distances = (unsigned long long int**)malloc(sizeof(unsigned long long int) * max_distance);
	for(size_t i = 0; i < n_threads; i++){
		distances[i] = distance_entries + i*max_distance;
	}

	unsigned int max_chunk = 600000;
  unsigned int chunk_size;
//	printf("File size: %i\n", file_size);
//Case 1: Whole file read
	if (file_size < max_chunk){
		chunk_size = file_size;
	}

//Case 2: Only two chunks
	else if(file_size*2 < max_chunk){
		chunk_size = max_chunk/2;
	}

//Case 3: Several chunks
	else{
		chunk_size = max_chunk/3;
	}

//	printf("Chunk size: %i\n", chunk_size);
  short *base_entries = (short *)malloc(sizeof(short) * chunk_size * 3);
  short **base = (short **)malloc(sizeof(short *) * chunk_size);
  //short *buffer_entries = (short *)aligned_alloc(sizeof(short)*3, sizeof(short) * chunk_size * 3);
  //short **buffer = (short **)malloc(sizeof(short *) * chunk_size);
  short *moving_entries = (short *)malloc(sizeof(short) * chunk_size * 3);
  short **moving = (short **)malloc(sizeof(short *) * chunk_size);

  for (unsigned int i = 0; i < chunk_size; i++) {
    base[i] = base_entries + i * 3;
    //buffer[i] = buffer_entries + i * 3;
    moving[i] = moving_entries + i * 3;
  }

  unsigned int current_cell_point = 0;
  unsigned int file_point = 0;
  //short *current_cell = (short *)aligned_alloc(sizeof(short) * 3, sizeof(short)*3);
  

	//float *compare_cell = (float*)malloc(sizeof(float)*3);
  char *chunk_string =
      (char *)malloc(sizeof(char) * (chunk_size * line_size + 1));
  unsigned int compare;

  char* format = (char*)malloc(sizeof(char)*100);
	sprintf(format, "%%%ic", chunk_size*line_size);

	short* dist_cell;
	short distance_one; short distance_two; short distance_three; short distance_four;
#pragma omp parallel shared(base, base_entries, distances, chunk_string, \
                                 current_cell_point, file_point, format) \
					  private(dist_cell, distance_one, distance_two, distance_three, distance_four)
{
//Case 1: Whole file read
	if (file_size < max_chunk){
		//chunk_size = file_size;
		fscanf(file, format, chunk_string);
#pragma omp for
      for (unsigned int i = 0; i < chunk_size; i++) {
				//printf("%i\n", i);
				char a[2]; char b[3]; char c[2]; char d[3]; char e[2]; char f[3];
				char ca; char cb; char cc;
//Could be optimized by not using sscanf
        sscanf(chunk_string + i * line_size, "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", &ca, a, b, &cb, c, d,
               &cc, e, f);
				base[i][0] = convert_to_short(a,b) * (ca == '+' ? 1:-1);
				base[i][1] = convert_to_short(c,d) * (cb == '+' ? 1:-1);
				base[i][2] = convert_to_short(e,f) * (cc == '+' ? 1:-1);
			}
			dist_cell = (short*)aligned_alloc(sizeof(short)*3, sizeof(short)*3);
			short thread;
			for (unsigned int i = 0; i < chunk_size; i++){
#pragma omp for
				for(unsigned int j = i+1; j < chunk_size-3; j+=4){
					//printf("base[i]: %i %i %i\n", base[i][0], base[i][1], base[i][2]);
					//printf("base[i+j]: %i %i %i\n", base[i+j][0], base[i+j][1], base[i+j][2]);
					distance_one = get_distance(base[i], base[j])/10;
					distance_two = get_distance(base[i], base[j+1])/10;	
					distance_three = get_distance(base[i], base[j+2])/10;
					distance_four = get_distance(base[i], base[j+3])/10;	
					thread = omp_get_thread_num();
	        distances[thread][distance_one] += 1;
	        distances[thread][distance_two] += 1;
	        distances[thread][distance_three] += 1;
	        distances[thread][distance_four] += 1;

//#pragma omp atomic
				}
			}
			free(dist_cell);

	}

//Case 2: Only two chunks
	else if(file_size*2 < max_chunk){
		//chunk_size = max_chunk/2;
	}

//Case 3: Several chunks
	else{
		//chunk_size = max_chunk/3;
	}
}
  //free(current_cell);
  //free(compare_cell);  
  free(base_entries);
  free(base);
  //free(buffer_entries);
  //free(buffer);
  free(moving_entries);
  free(moving);
	free(format);
	free(chunk_string);

	unsigned long long int* result = (unsigned long long*)malloc(sizeof(unsigned long long)*max_distance);
	for (size_t i = 0; i < max_distance; i++){	
		result[i] = 0;
	}
		for(size_t j = 0; j < n_threads; j++){
			for(size_t i = 0; i < max_distance; i++){
					result[i] += distances[j][i];
		}
	}


  for (size_t i = 0; i < max_distance; i++) {
    if (result[i] != 0)
      printf("%05.2f %llu\n", i / 100.0, result[i]);
			//printf("%i : %i (%i)\n", i, i*i, i*i - (i-1)*(i-1));
  }
	free(distance_entries);
  free(distances);
	free(result);

}
/*
  for (size_t i = 0; i < file_size; i++) {
#pragma omp single 
    {
    printf("On cell %i\n", i);
    fseek(file, current_cell_point*line_size, SEEK_SET);
		char a[2]; char b[3]; char c[2]; char d[3]; char e[2]; char f[3];
		char ca; char cb; char cc;

    fscanf(file, "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", &ca, a, b, &cb, c, d,
               &cc, e, f);
		current_cell[0] = convert_to_short(a,b) * (ca == '+' ? 1:-1);
		current_cell[1] = convert_to_short(c,d) * (cb == '+' ? 1:-1);
		current_cell[2] = convert_to_short(e,f) * (cc == '+' ? 1:-1);

	//	for(size_t j = 0; j < 3; j++){
	//		current_cell[j] = -current_cell[j];
	//	}
    current_cell_point = ftell(file)/line_size;
    //printf("current cell point: %i\n", current_cell_point);
    }
    file_point = current_cell_point;
    //fseek(file, current_cell_point*line_size + line_size, SEEK_SET);
    while(file_point < file_size) {
#pragma omp single
    {
//fgets(chunk_string, chunk_size * line_size + 1, file);
    file_point = ftell(file)/line_size;
    //file_point = current_cell_point;
		fscanf(file, format, chunk_string);
    //printf("%i\n", file_point);
    compare = file_size - file_point;
    //printf("compare: %i\n", compare);
    if (chunk_size > file_size - file_point){
      compare = file_size - file_point;
    } else {
      compare = chunk_size;
    }
    }
// Loop over all whole chunks
#pragma omp for //simd simdlen(3) aligned(current_cell, all_cells:3)
      for (size_t k = 0; k < compare; k++) {
				char a[2]; char b[3]; char c[2]; char d[3]; char e[2]; char f[3];
				char ca; char cb; char cc;
				dist_cell = (short*)aligned_alloc(sizeof(short)*3, sizeof(short)*3);
//      printf("%i\n", k);
//          fscanf(file, "%f %f %f\n", &compare_cell[0],
//               &compare_cell[1], &compare_cell[2]);

        sscanf(chunk_string + k * line_size, "%c%2c.%3c %c%2c.%3c %c%2c.%3c\n", &ca, a, b, &cb, c, d,
               &cc, e, f);
//        printf("%f\n", get_distance(current_cell, all_cells[k])); 
				all_cells[k][0] = convert_to_short(a,b) * (ca == '+' ? 1:-1);
				all_cells[k][1] = convert_to_short(c,d) * (cb == '+' ? 1:-1);
				all_cells[k][2] = convert_to_short(e,f) * (cc == '+' ? 1:-1);
				short distance = get_distance(current_cell, all_cells[k], dist_cell)/10;
		//		printf("distance: %i\n", distance);

#pragma omp atomic
        distances[distance] += 1;
				free(dist_cell);
        }
      }  
    }

#pragma omp for
    for (size_t j = 0; j < file_size % chunk_size; ++j) {
      // Loop over remaining lines
      sscanf(chunk_string + j * line_size, "%f %f %f\n", &all_cells[j][0],
             &all_cells[j][1], &all_cells[j][2]);
      short distance = (short)(100 * get_distance(current_cell, all_cells[j]));
#pragma omp atomic
      distances[distance] += 1;
    }
*/

