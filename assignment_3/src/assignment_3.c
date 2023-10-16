#include <stdio.h>
#include <stdlib.h>
#include <threads.h>
#include <math.h>
#include <string.h>

#define PI 3.14

void
create_color_palette(char* color_palette, const unsigned short line_size){
	memcpy(color_palette,"158   1  66\n",line_size);
	memcpy(color_palette + line_size,"213  62  79\n",line_size);
	memcpy(color_palette + 2*line_size, "244 109  67\n",line_size);
	memcpy(color_palette + 3*line_size, "253 174  97\n", line_size);
	memcpy(color_palette + 4*line_size, "254 224 139\n",line_size);
	memcpy(color_palette + 5*line_size, "230 245 152\n",line_size);
	memcpy(color_palette + 6*line_size, "171 221 164\n",line_size);
	memcpy(color_palette + 7*line_size, "102 194 165\n", line_size);
	memcpy(color_palette + 8*line_size, " 50 136 189\n", line_size);
	memcpy(color_palette + 9*line_size, " 94  79 162\n",line_size);
}

//In case we'd want to change to something else
static inline
double
calc_derivative(double x, const unsigned short degree){
	return degree*pow(x, degree-1);
}
typedef struct{
	unsigned int start;
	unsigned int stop;
	char** attractors;
	char** convergences;
	mtx_t* mtx;
	cnd_t* cnd;
} compute_struct_t


int
compute_thread(void* args){
	
	return 0;
}

int
write_thread(void* args){
	return 0;
}


int
main(int argc, char* argv[]){
	//Read command lines
	unsigned short n_threads;
	unsigned short n_lines;
	const unsigned short degree = atoi(argv[3]);

	char first_char = argv[1][1];
	
	if(first_char == 't'){
		n_threads = atoi(argv[1]+2);
		n_lines = atoi(argv[2]+2);
	}
	else{
		n_threads = atoi(argv[2]+2);
		n_lines = atoi(argv[1]+2);
	}

	//Initialize constants
	const unsigned short max_degree = 10;
	const unsigned short line_size = 12;
	const double x_max = 2;
	const double x_min = -2;
	const double d_x = (x_max - x_min)/n_lines;
	
	char* color_palette = (char*)malloc(sizeof(char)*max_degree*line_size+1);
	create_color_palette(color_palette, line_size);
	
	double* root_list = (double*)malloc(sizeof(double)*max_degree);
	for(size_t i = 0; i < max_degree; i++){
		root_list[i] = cos(2*PI*i/(double)max_degree);
	}

	FILE* file = fopen("newton_attractors_xd.ppm", "w");

	fprintf(file, "P3\n%d %d\n255\n", n_lines, n_lines);
//	char* file_header = (char*)malloc(sizeof(char)*100);
//	sprintf(file_header, "P3\n%d %d\nM\n", n_lines, n_lines);
//	fwrite(file_header, sizeof(char), 9, file);
//	free(file_header);

	printf(color_palette);


	//Global result arrays
	unsigned char* attractors_entries = (unsigned char*)malloc(sizeof(char)*n_lines*n_lines+1);
	unsigned char** attractors = (unsigned char**)malloc(sizeof(char*)*n_lines);
	unsigned char* convergences_entries = (unsigned char*)malloc(sizeof(char)*n_lines*n_lines+1);
	unsigned char** convergences = (unsigned char**)malloc(sizeof(char*)*n_lines);

	for(unsigned int i = 0; i < n_lines; i++){
		attractors[i] = attractors_entries + i*n_lines;
		convergences[i] = convergences_entries + i*n_lines;
	}

	//Threading: Use threads to calculate one row at a time for the arrays
	//Connect with a write thread that writes row for row
	//Need to keep track of which line to set the file pointer correctly
	//Set up the framwork for computational threads, but just return a set value to check synchronization and writing







	//Test writing to file
  //Assume n_lines = 10 and write all them colors for the number of lines

	for(unsigned int i = 0; i < n_lines; i++){
		for(unsigned int j = 0; j < n_lines; j++){
			attractors[i][j] = i;
		}
	}


	for(unsigned int i = 0; i < n_lines; i++){
		for(unsigned int j = 0; j < n_lines; j++){
			fwrite(color_palette+attractors[i][j]*line_size, sizeof(char), line_size, file);
		}
	}



	fclose(file);
	free(color_palette);
	free(root_list);
	free(attractors);
	free(attractors_entries);
	free(convergences);
	free(convergences_entries);
	return 0;
}
