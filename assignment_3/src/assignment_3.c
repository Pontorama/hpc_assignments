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
//double x_min; double x_max;
//double dx;
//double root_list[degree][2];
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
		sprintf(grey_palette+line_size*i, "%3i %3i %3i", i,i,i);
	}
/*  memcpy(grey_palette, "25  25  25\n", line_size);
  memcpy(grey_palette + line_size, " 50  50  50\n", line_size);
  memcpy(grey_palette + 2 * line_size, " 75  75  75\n", line_size);
  memcpy(grey_palette + 3 * line_size, "100 100 100\n", line_size);
  memcpy(grey_palette + 4 * line_size, "125 125 125\n", line_size);
  memcpy(grey_palette + 5 * line_size, "150 150 150\n", line_size);
  memcpy(grey_palette + 6 * line_size, "175 175 175\n", line_size);
  memcpy(grey_palette + 7 * line_size, "200 200 200\n", line_size);
  memcpy(grey_palette + 8 * line_size, "225 225 225\n", line_size);
  memcpy(grey_palette + 9 * line_size, "250 250 250\n", line_size);
  memcpy(grey_palette + 10 * line_size, "  0   0   0\n", line_size);
*/
}
// In case we'd want to change to something else
static inline double calc_derivative(double x, const unsigned short degree) {
  return degree * pow(x, degree - 1);
}

static inline void cart_to_polar(const double x, const double y, double* r, double* theta){
	*r = sqrt(x*x+y*y);
	*theta = asin(y/ *r);
}


static inline int newton_iterations(double complex z, double complex* root_list,
		unsigned char* attractor, unsigned char* convergence){
	unsigned char root;
	double complex check;
//	double complex f_z;
//	double complex f_z_d;
	double complex norm;
	double complex conj_z;
	double complex div;
	*attractor = 10;
	*convergence = 0;
	for(size_t i = 0; i < MAX_ITERS; i++){
		for(size_t j = 0; j < degree; j++){
			check = (root_list[j] - z)*1000;
//				if(cabs(check) < 0.001){
			if(creal(check*conj(check)) <= 0.001){
//				printf("Found a root!(%i, %i)\n", j, i);
				*attractor = (unsigned char)j;
				*convergence = (unsigned char)i;
				return 0;
			}
		}
		if(creal(z*conj(z))*1000 <= 0.001){
				*attractor = 10;
				*convergence = 0;
				return 0;
			}
		if(creal(z) > 10000000000 || cimag(z) > 10000000000 ||
					creal(z) < -10000000000 || cimag(z) < -10000000000){
				*attractor = 10;
				*convergence = 0;
				return 0;
			}

			switch(degree){
				case 1:
					z = 1;
					break;
				case 2:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 0.5*(z + div);
					break;
				case 3:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 1/3.*(2*z+div*div);
					break;
				case 4:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 0.25*(3*z+div*div*div);
					break;
				case 5:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 0.2*(4*z+div*div*div*div);
					break;
				case 6:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 1/6.*(5*z+div*div*div*div*div);
					break;
				case 7:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 1/7.*(6*z+div*div*div*div*div*div);
					break;
				case 8:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 0.125*(7*z+div*div*div*div*div*div*div);
					break;
				case 9:
					conj_z = conj(z);
					norm = creal(z*conj_z);
					div = conj_z/norm;
					z = 1/9.*(8*z+div*div*div*div*div*div*div*div);
					break;
			}

/*			f_z = z;
			f_z_d = z;
			for(size_t j = 0; j < degree-2; j++){
				f_z_d *= f_z_d;
			}
			f_z = f_z_d*f_z_d;
*/
			//z = z - (cpow(z,degree)-1)/(degree*cpow(z,degree-1));
		
	}
	return 0;
	
}

typedef struct {
	unsigned int n_thread;
	unsigned int *n_written;
  unsigned int row_start;
	unsigned int row_end;
	unsigned int *row_index;
  unsigned char **attractors;
  unsigned char **convergences;
	double complex *root_list;
	bool *ready;
  mtx_t *mtx;
  cnd_t *cnd;
	double dx;
	double x_min;
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
	double complex* root_list = thrd_info->root_list;
	bool* ready = thrd_info->ready;
	cnd_t* cnd = thrd_info->cnd;
	mtx_t* mtx = thrd_info->mtx;
	unsigned int n_thread = thrd_info->n_thread;
	unsigned int* n_written = thrd_info->n_written;
	const double x_min = thrd_info->x_min;
	const double dx = thrd_info->dx;

	double x;
	double y;
	//double r;
	//double theta;
	double complex z;
	//char* attractor = (char*)malloc(sizeof(char)*(row_end-row_start)+1);
	for(unsigned int i = row_start; i < row_end; i++){
		for(unsigned int j = 0; j < n_lines; j++){
			x = x_min + i*dx;
			y = x_min + j*dx;
			//cart_to_polar(x,y,&r,&theta);
			z = x + I*y;
			int r;
			r = newton_iterations(z, root_list, &attractors[i][j], &convergences[i][j]);
			
//			attractors[i][j] = n_thread;
//			convergences[i][j] = n_thread;
		}
		mtx_lock(mtx);
//		printf("Done with row index %i\n");
//		*row_index = i;
//		*ready = true;
//		printf("%u\n",*n_written);
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
  //fseek(attractors_file, 15, SEEK_SET);
  for (unsigned int i = 0; i < n_lines; i++) {
    for (unsigned int j = 0; j < n_lines; j++) {
        fwrite(color_palette +
                   attractors[i][j] *
                       line_size,
               sizeof(char), line_size, attractors_file);
				fwrite(grey_palette + convergences[i][j] * line_size, sizeof(char),line_size, convergences_file);
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  // Read command lines
  degree = atoi(argv[3]);

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
  const double dx = (x_max - x_min) / n_lines;

  color_palette = (char *)malloc(sizeof(char) * max_degree * line_size + 1);
  create_color_palette(color_palette, line_size);

  grey_palette = (char *)malloc(sizeof(char) * MAX_ITERS * line_size + 1);
  create_grey_palette(grey_palette, line_size);


//	double *root_list_entries = (double*)malloc(sizeof(double)*degree*2);
//	double **root_list = (double**)malloc(sizeof(double*)*degree);
	
//	for(size_t i = 0; i < degree; i++){
//		root_list[i] = root_list_entries + 2*i;
//	}


  double complex* root_list = (double complex*)malloc(sizeof(double complex) * max_degree);
  for (size_t i = 0; i < degree; i++) {
    root_list[i] = cos(2 * PI * i/(double)degree) 
									+ I*sin(2* PI * i/(double)degree);
//		root_list[i] = sin(2* PI * i) / (double)degree;
  }

	char* a_file_name = (char*)malloc(sizeof(char)*100);
	char* c_file_name = (char*)malloc(sizeof(char)*100);
	sprintf(a_file_name, "newton_attractors_x%d.ppm", degree);
	sprintf(c_file_name, "newton_convergences_x%d.ppm", degree);

  attractors_file = fopen(a_file_name, "w");
	convergences_file = fopen(c_file_name, "w");

  fprintf(attractors_file, "P3\n%d %d\n255\n", n_lines, n_lines);
  fprintf(convergences_file, "P3\n%d %d\n255\n", n_lines, n_lines);

	free(a_file_name);
	free(c_file_name);

//	for(unsigned int i = 0; i < n_lines*n_lines; i++){
//		fwrite("           \n", sizeof(char), line_size, attractors_file);
//	}

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

	switch(degree){
		case 1:

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
	fclose(convergences_file);
  free(color_palette);
	free(grey_palette);
  free(root_list);
	//free(root_list_entries);
  free(attractors);
  free(attractors_entries);
  free(convergences);
  free(convergences_entries);
  return 0;
}
