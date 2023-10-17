#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

#define PI 3.14
unsigned int n_lines;
unsigned short n_threads;
char *color_palette;
unsigned short line_size;
FILE *attractors_file;

void create_color_palette(char *color_palette, const unsigned short line_size) {
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
  memcpy(color_palette + 9 * line_size, "  0   0   0\n", line_size);
}

// In case we'd want to change to something else
static inline double calc_derivative(double x, const unsigned short degree) {
  return degree * pow(x, degree - 1);
}
typedef struct {
	unsigned int n_thread;
	unsigned int *n_written;
  unsigned int row_start;
	unsigned int row_end;
	unsigned int *row_index;
  unsigned char **attractors;
  unsigned char **convergences;
	bool *ready;
  mtx_t *mtx;
  cnd_t *cnd;
} compute_struct_t;

typedef struct {
  unsigned int *row_index;
  unsigned char **attractors;
  unsigned char **convergences;
	unsigned int *n_written;
  bool *ready;
  cnd_t *cnd;
  mtx_t *mtx;
} write_struct_t;

int compute_thread(void *args) { 
	const compute_struct_t *thrd_info = (compute_struct_t*)args;
	const unsigned int row_start = thrd_info->row_start;
	const unsigned int row_end = thrd_info->row_end;
	unsigned int* row_index = thrd_info->row_index;
	unsigned char** attractors = thrd_info->attractors;
	unsigned char** convergences = thrd_info->convergences;
	bool* ready = thrd_info->ready;
	cnd_t* cnd = thrd_info->cnd;
	mtx_t* mtx = thrd_info->mtx;
	unsigned int n_thread = thrd_info->n_thread;
	unsigned int* n_written = thrd_info->n_written;
	//char* attractor = (char*)malloc(sizeof(char)*(row_end-row_start)+1);
	for(unsigned int i = row_start; i < row_end; i++){
		for(unsigned int j = 0; j < n_lines; j++){
			attractors[i][j] = n_thread;
		}
		mtx_lock(mtx);
//		printf("Done with row index %i\n");
//		*row_index = i;
//		*ready = true;
		printf("%u\n",*n_written);
		*n_written += 1;
		mtx_unlock(mtx);
		cnd_signal(cnd);
	}

	return 0; 
}

int write_thread(void *args) {
  const write_struct_t *thrd_write_info = (write_struct_t *)args;
//  unsigned int n_written_rows = 0;
	unsigned int* row_index;
	unsigned char** attractors = thrd_write_info->attractors;
	unsigned char** convergences = thrd_write_info->convergences;
	bool* ready = thrd_write_info->ready;
	mtx_t* mtx = thrd_write_info->mtx;
	cnd_t* cnd = thrd_write_info->cnd;
	unsigned int *n_written = thrd_write_info->n_written;

  mtx_lock(mtx);
  while (*n_written != n_lines) {
    cnd_wait(cnd, mtx);
/*    if (*ready) {
			printf("Line %d ready\n", *row_index);
			fseek(attractors_file, line_size*(*row_index)*n_lines, SEEK_SET);
      for (size_t j = 0; j < n_lines; ++j) {
        fwrite(color_palette +
                   attractors[*row_index][j] *
                       line_size,
               sizeof(char), line_size, attractors_file);
      }
      n_written_rows++;
			*ready = false;
    }
*/
  }
  mtx_unlock(mtx);

	printf("Completed all threads\n");
  fseek(attractors_file, 15, SEEK_SET);
  for (unsigned int i = 0; i < n_lines; i++) {
    for (unsigned int j = 0; j < n_lines; j++) {
        fwrite(color_palette +
                   attractors[i][j] *
                       line_size,
               sizeof(char), line_size, attractors_file);
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  // Read command lines
  const unsigned short degree = atoi(argv[3]);

  char first_char = argv[1][1];

  if (first_char == 't') {
    n_threads = atoi(argv[1] + 2);
    n_lines = atoi(argv[2] + 2);
  } else {
    n_threads = atoi(argv[2] + 2);
    n_lines = atoi(argv[1] + 2);
  }

  // Initialize constants
  const unsigned short max_degree = 10;
  line_size = 12;
  const double x_max = 2;
  const double x_min = -2;
  const double d_x = (x_max - x_min) / n_lines;

  color_palette = (char *)malloc(sizeof(char) * max_degree * line_size + 1);
  create_color_palette(color_palette, line_size);

  double *root_list = (double *)malloc(sizeof(double) * max_degree);
  for (size_t i = 0; i < max_degree; i++) {
    root_list[i] = cos(2 * PI * i / (double)max_degree);
  }

  attractors_file = fopen("newton_attractors_xd.ppm", "w");

  fprintf(attractors_file, "P3\n%d %d\n255\n", n_lines, n_lines);
	for(unsigned int i = 0; i < n_lines*n_lines; i++){
		fwrite("           \n", sizeof(char), line_size, attractors_file);
	}
  //	char* file_header = (char*)malloc(sizeof(char)*100);
  //	sprintf(file_header, "P3\n%d %d\nM\n", n_lines, n_lines);
  //	fwrite(file_header, sizeof(char), 9, file);
  //	free(file_header);

  printf(color_palette);

  // Global result arrays
  unsigned char *attractors_entries =
      (unsigned char *)malloc(sizeof(char) * n_lines * n_lines + 1);
  unsigned char **attractors =
      (unsigned char **)malloc(sizeof(char *) * n_lines);
  unsigned char *convergences_entries =
      (unsigned char *)malloc(sizeof(char) * n_lines * n_lines + 1);
  unsigned char **convergences =
      (unsigned char **)malloc(sizeof(char *) * n_lines);

  for (unsigned int i = 0; i < n_lines; i++) {
    attractors[i] = attractors_entries + i * n_lines;
    convergences[i] = convergences_entries + i * n_lines;
  }

  // Threading: Use threads to calculate one row at a time for the arrays
  // Connect with a write thread that writes row for row
  // Need to keep track of which line to set the file pointer correctly
  // Set up the framwork for computational threads, but just return a set value
  // to check synchronization and writing

	mtx_t mtx;
	mtx_init(&mtx, mtx_plain);

	cnd_t cnd;
	cnd_init(&cnd);

	bool ready = false;
	unsigned int row_index = 0;
	
	thrd_t threads[n_threads];
	thrd_t w_thread;
	compute_struct_t thrds_info[n_threads];
	write_struct_t thrd_info_write;
	const unsigned int n_lines_loc = n_lines/n_threads;
	unsigned int row_start;
	unsigned int n_written = 0;
	for(int i = 0, row_start = 0;i < n_threads; ++i, row_start += n_lines_loc){
		thrds_info[i].attractors = attractors;
		thrds_info[i].convergences = convergences;
		thrds_info[i].n_written = &n_written;
		thrds_info[i].n_thread = i;
		thrds_info[i].row_start = row_start;
		thrds_info[i].row_end = i != n_threads - 1 ? row_start + n_lines_loc : n_lines;
		thrds_info[i].ready = &ready;
		thrds_info[i].row_index = &row_index;
		thrds_info[i].mtx = &mtx;
		thrds_info[i].cnd = &cnd;

		int r = thrd_create(threads+i, compute_thread, (void*)(thrds_info+i));
		if(r != thrd_success){
			printf("Failed to create thread.\n");
			exit(1);
		}
		printf("Successfully created thread %d\n", i);
	}

	{
		thrd_info_write.attractors = attractors;
		thrd_info_write.convergences = convergences;
		thrd_info_write.n_written = &n_written;
		thrd_info_write.ready = &ready;
		thrd_info_write.row_index = &row_index;
		thrd_info_write.mtx = &mtx;
		thrd_info_write.cnd = &cnd;

		int r = thrd_create(&w_thread, write_thread, (void*)(&thrd_info_write));
		if(r != thrd_success){
			printf("Failed to create thread.\n");
			exit(1);
		}
		printf("Successfully created write thread\n");
	}

	for(int i = 0; i < n_threads; i++){
		int r;
		thrd_join(threads[i], &r);
	}
	{
		int r;
		thrd_join(w_thread, &r);
	}

	mtx_destroy(&mtx);
	cnd_destroy(&cnd);
  // Test writing to file
  // Assume n_lines = 10 and write all them colors for the number of lines

/*
  for (unsigned int i = 0; i < n_lines; i++) {
    for (unsigned int j = 0; j < n_lines; j++) {
      attractors[i][j] = i;
    }
  }

  for (unsigned int i = 0; i < n_lines; i++) {
    for (unsigned int j = 0; j < n_lines; j++) {
    }
  }
*/
  fclose(attractors_file);
  free(color_palette);
  free(root_list);
  free(attractors);
  free(attractors_entries);
  free(convergences);
  free(convergences_entries);
  return 0;
}
