#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>
#include <complex.h>

#define PI 3.14159265359
#define MAX_ITERS 128

unsigned int n_lines;
unsigned short n_threads;
char *color_palette;
char *grey_palette;
unsigned short line_size;
unsigned short degree;

char* write_attractors;
char* write_convergence;

FILE *attractors_file;
FILE *convergences_file;

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
  memcpy(color_palette + 10 * line_size, "  0   0   0\n", line_size);
}
void create_grey_palette(char *grey_palette, const unsigned short line_size) {
	for(size_t i = 0; i < MAX_ITERS; i++){
		sprintf(grey_palette+line_size*i, "%3i %3i %3i\n", 2*i,2*i,2*i);
	}

}


static inline int newton_iterations_4(double complex z[4], const double complex* root_list,
		 const int row_number, const int col_number[4]){

	unsigned short loc_degree = degree;
	double check_x[4];
	double check_y[4];
	double complex check[4];
	double norm[4];
	double complex conj_z[4];
	double complex z_k[4];
	double complex div[4];
	double complex div_r[4];
	//double tmp;
	short bools[4];
	for(int k = 0; k < 4; k++){
		bools[k] = 0;
	}
	short done = 0;
	size_t i = 0;
	if(loc_degree == 1){
		//There is only one possibility for degree = 1
		for(size_t l = 0; l < 4; l++){
		strncpy(write_attractors + line_size*(n_lines*row_number + col_number[l]),
														color_palette, line_size);
		strncpy(write_convergence + line_size*(n_lines*row_number + col_number[l]),
														grey_palette, line_size);
		}
	}
	else{
		for(i = 0; i < MAX_ITERS; i++){
			done = 0;
			for(size_t j = 0; j < 4; j++){
				done += bools[j];
				if(done == 4){
					return 0;
				}
			}

			conj_z[0] = conj(z[0]);
			conj_z[1] = conj(z[1]);
			conj_z[2] = conj(z[2]);
			conj_z[3] = conj(z[3]);
			for(size_t l = 0; l < 4; l++){
					if(bools[l] == 0){
					norm[l] = creal(z[l]*conj_z[l]);
					if(norm[l] - 1 < 0.001){

					for(size_t j = 0; j < loc_degree; j++){
							check[l] = (root_list[j]-z[l])*1000000;
							check_x[l] = creal(check[l]);
							check_y[l] = cimag(check[l]);
							if(check_x[l]*check_x[l] + check_y[l]*check_y[l] < 1){
								strncpy(write_attractors + line_size*(n_lines*row_number + col_number[l]),
											color_palette + line_size*j, line_size);
								strncpy(write_convergence + line_size*(n_lines*row_number + col_number[l]),
											grey_palette + line_size*i, line_size);

								bools[l] = 1;
								break;
								}
							}
					}
					if(bools[l] == 1)
						continue;
					if(norm[l] < 0.000001){
										strncpy(write_attractors + line_size*(n_lines*row_number + col_number[l]),
														color_palette + line_size*10, line_size);
										strncpy(write_convergence + line_size*(n_lines*row_number + col_number[l]),
														grey_palette, line_size);
								bools[l] = 1;
								continue;
								}
						if(norm[l] > 1e20){
										strncpy(write_attractors + line_size*(n_lines*row_number + col_number[l]),
														color_palette + line_size*10, line_size);
										strncpy(write_convergence + line_size*(n_lines*row_number + col_number[l]),
														grey_palette, line_size);
									bools[l] = 1;
									continue;
								}
							}
						}

/*			if (bools[0] == 0)
				norm[0] = creal(z[0]*conj_z[0]);
			if (bools[1] == 0)
				norm[1] = creal(z[1]*conj_z[1]);
			if (bools[2] == 0)
				norm[2] = creal(z[2]*conj_z[2]);
			if (bools[3] == 0)
				norm[3] = creal(z[3]*conj_z[3]);
*/
			if (bools[0] == 0)
				div[0] = conj_z[0]/norm[0];
			if (bools[1] == 0)
				div[1] = conj_z[1]/norm[1];				
			if (bools[2] == 0)
				div[2] = conj_z[2]/norm[2];
			if (bools[3] == 0)
				div[3] = conj_z[3]/norm[3];
			for(size_t k = 0; k < 4; k++){
				div_r[k] = 1;
			}
			for(size_t k = 0; k < loc_degree - 1; ++k){
				if (bools[0] == 0)
						div_r[0] *= div[0];
				if (bools[1] == 0)
						div_r[1] *= div[1];
				if (bools[2] == 0)
						div_r[2] *= div[2];
				if (bools[3] == 0)
						div_r[3] *= div[3];
			}
			if (bools[0] == 0)
				z[0] = ((loc_degree-1)*z[0]+div_r[0])/loc_degree;
			if (bools[1] == 0)
				z[1] = ((loc_degree-1)*z[1]+div_r[1])/loc_degree;
			if (bools[2] == 0)
				z[2] = ((loc_degree-1)*z[2]+div_r[2])/loc_degree;
			if (bools[3] == 0)
				z[3] = ((loc_degree-1)*z[3]+div_r[3])/loc_degree;
		}
	}
	
	if(i == MAX_ITERS){
			for(size_t k = 0; k < 4; k++){
				if(bools[k] == 0){
					strncpy(write_attractors + line_size*(n_lines*row_number + col_number[k]),
										color_palette + line_size*10, line_size);
					strncpy(write_convergence + line_size*(n_lines*row_number + col_number[k]),
										grey_palette, line_size);
			}
		}
	}

	return 0;
	
}
static inline int newton_iterations(double complex z, const double complex* root_list,
		 const int row_number, const int col_number){
	unsigned short loc_degree = degree;
	double check_x;
	double check_y;
	double complex check;
	double norm;
	double complex conj_z;
	double complex z_k;
	double complex div;
	double complex div_r;
	size_t i;
	if(degree == 1){
	 	z=1;
	}
	else{
		for(i = 0; i < MAX_ITERS; i++){
			conj_z = conj(z);
			norm = creal(z*conj_z);
			if(norm - 1 < 0.001){
				for(size_t j = 0; j < degree; j++){
						check = (root_list[j]-z)*1000;
						check_x = creal(check);
						check_y = cimag(check);
						if(check_x*check_x + check_y*check_y < 0.001){
							strncpy(write_attractors + line_size*(n_lines*row_number + col_number),
											color_palette + line_size*j, line_size);
							strncpy(write_convergence + line_size*(n_lines*row_number + col_number),
											grey_palette + line_size*i, line_size);
						return 0;
						}
				}
			}
			if(norm <= 0.001){
						strncpy(write_attractors + line_size*(n_lines*row_number + col_number),
										color_palette + line_size*10, line_size);
						strncpy(write_convergence + line_size*(n_lines*row_number + col_number),
										grey_palette, line_size);

					return 0;
				}
			if(norm > 1e20){
						strncpy(write_attractors + line_size*(n_lines*row_number + col_number),
										color_palette + line_size*10, line_size);
						strncpy(write_convergence + line_size*(n_lines*row_number + col_number),
										grey_palette, line_size);

					return 0;
				}
			div = conj_z/norm;
			div_r = 1;
			for(size_t k = 0; k < loc_degree - 1; ++k){
				div_r *= div;
			}
			z = ((loc_degree-1)*z+div_r)/loc_degree;
		}
	}
	
	if(i == MAX_ITERS){

				strncpy(write_attractors + line_size*(n_lines*row_number + col_number),
									color_palette + line_size*10, line_size);
				strncpy(write_convergence + line_size*(n_lines*row_number + col_number),
									grey_palette, line_size);

	}
	return 0;
	
}
	
typedef struct {
	unsigned int *n_written;
  unsigned int row_start;
	unsigned int row_end;
	double complex *root_list;
  mtx_t *mtx;
  cnd_t *cnd;
	double dx;
	double x_min;
} compute_struct_t;

typedef struct {
	unsigned int *n_written;
  cnd_t *cnd;
  mtx_t *mtx;
} write_struct_t;

int compute_thread(void *args) { 
	const compute_struct_t *thrd_info = (compute_struct_t*)args;
	const unsigned int row_start = thrd_info->row_start;
	const unsigned int row_end = thrd_info->row_end;
	double complex* root_list = thrd_info->root_list;
	cnd_t* cnd = thrd_info->cnd;
	mtx_t* mtx = thrd_info->mtx;
	unsigned int* n_written = thrd_info->n_written;
	const double x_min = thrd_info->x_min;
	const double dx = thrd_info->dx;
	double x;
	double y[4];
  double y_lone;
	double complex z[4];
	double complex z_lone;
	int remainder = n_lines % 4;
	int main_loop_end = n_lines-remainder;
	for(unsigned int i = row_start; i < row_end; i++){
		x = x_min + i*dx;
	for(unsigned int j = 0; j < main_loop_end; j+=4){
			y[0] = (x_min + j*dx);
			y[1] = (x_min + (j+1)*dx);
			y[2] = (x_min + (j+2)*dx);
			y[3] = (x_min + (j+3)*dx);

			z[0] = x + I*y[0];
			z[1] = x + I*y[1];
			z[2] = x + I*y[2];
			z[3] = x + I*y[3];
			int r1;
			int col_number[] = {j,j+1,j+2,j+3};
			r1 = newton_iterations_4(z, root_list ,i ,col_number);
		}
/*
		for(unsigned int j = 0; j < main_loop_end; j++){
			y_lone = x_min + j*dx;
			z_lone = x + I*y_lone;
			int r;
			r = newton_iterations(z_lone,root_list,i,j);
		}
*/
		mtx_lock(mtx);
		*n_written += 1;
		mtx_unlock(mtx);
		cnd_signal(cnd);
	}

	return 0; 
}

int write_thread(void *args) {
  const write_struct_t *thrd_write_info = (write_struct_t *)args;
	mtx_t* mtx = thrd_write_info->mtx;
	cnd_t* cnd = thrd_write_info->cnd;
	unsigned int *n_written = thrd_write_info->n_written;

  mtx_lock(mtx);
  while (*n_written != n_lines) {
    cnd_wait(cnd, mtx);
  }
  mtx_unlock(mtx);

	printf("Completed all threads\n");
		write_attractors[n_lines*n_lines*line_size] = '\0';
		write_convergence[n_lines*n_lines*line_size] = '\0';
		fwrite(write_attractors, sizeof(char), n_lines*n_lines*line_size,attractors_file);
		fwrite(write_convergence, sizeof(char), n_lines*n_lines*line_size,convergences_file);
  return 0;
}

int main(int argc, char *argv[]) {
  // Read command lines
  degree = atoi(argv[3]);
	if(degree > 9){
		printf("Invalid degree\n");
		exit(1);
	}
  const char first_char = argv[1][1];

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
  const double dx = (x_max - x_min) / n_lines;

  color_palette = (char *)malloc(sizeof(char) * max_degree * line_size + 1);
  create_color_palette(color_palette, line_size);

  grey_palette = (char *)malloc(sizeof(char) * MAX_ITERS * line_size + 1);
  create_grey_palette(grey_palette, line_size);

  double complex* root_list = (double complex*)malloc(sizeof(double complex) * max_degree);
  for (size_t i = 0; i < degree; i++) {
    root_list[i] = cos(2 * PI * i/(double)degree) 
									+ I*sin(2* PI * i/(double)degree);
  }

	char* a_file_name = (char*)malloc(sizeof(char)*100);
	char* c_file_name = (char*)malloc(sizeof(char)*100);
	sprintf(a_file_name, "newton_attractors_x%d.ppm", degree);
	sprintf(c_file_name, "newton_convergence_x%d.ppm", degree);

  attractors_file = fopen(a_file_name, "w");
	convergences_file = fopen(c_file_name, "w");

  fprintf(attractors_file, "P3\n%d %d\n255\n", n_lines, n_lines);
  fprintf(convergences_file, "P3\n%d %d\n255\n", n_lines, n_lines);

	free(a_file_name);
	free(c_file_name);

	write_attractors = (char*)malloc(sizeof(char)*(n_lines*n_lines*line_size+1));
	write_convergence = (char*)malloc(sizeof(char)*(n_lines*n_lines*line_size+1));



	mtx_t mtx;
	mtx_init(&mtx, mtx_plain);

	cnd_t cnd;
	cnd_init(&cnd);
	
	thrd_t threads[n_threads];
	thrd_t w_thread;
	compute_struct_t thrds_info[n_threads];
	write_struct_t thrd_info_write;
	const unsigned int n_lines_loc = n_lines/n_threads;
	unsigned int row_start;
	unsigned int n_written = 0;

	for(int i = 0, row_start = 0 ; i < n_threads; ++i, row_start += n_lines_loc){
		thrds_info[i].n_written = &n_written;
		thrds_info[i].row_start = row_start;
		thrds_info[i].row_end = i != n_threads - 1 ? row_start + n_lines_loc : n_lines;
		thrds_info[i].mtx = &mtx;
		thrds_info[i].cnd = &cnd;
		thrds_info[i].root_list = root_list;
		thrds_info[i].dx = dx;
		thrds_info[i].x_min = x_min;

		int r = thrd_create(threads+i, compute_thread, (void*)(thrds_info+i));
		if(r != thrd_success){
			printf("Failed to create thread.\n");
			exit(1);
		}
		printf("Successfully created thread %d\n", i);
	}

	{
		thrd_info_write.n_written = &n_written;
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

  fclose(attractors_file);
	fclose(convergences_file);

  free(color_palette);
	free(grey_palette);
  free(root_list);
	free(write_convergence);
	free(write_attractors);
  return 0;
}
