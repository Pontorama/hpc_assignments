#include <stdio.h>
#include <stdlib.h>

int** create_matrix(int size){
  int* asentries = (int*) malloc(sizeof(int) * size * size);
  int** as = (int**) malloc(sizeof(int*) * size);

  for(size_t i = 0, j = 0; i < size; ++i, j += size){
    as[i] = asentries + j;
  }

  for(size_t i = 0; i < size; ++i){
    for(size_t j = 0; j < size; ++j){
      as[i][j] = j*i;
    }
  }

  return as;
}

int main(){
  int size = 10;

  int** m = create_matrix(size);
  // Open file for writing
  FILE* file = fopen("matrix_data.dat", "w");
  // Check return value
  if(file == NULL){
    printf("Error opening file\n");
    return -1;
  }
  // Can write like this because matrix is contigous in memory
  // Otherwise would have to use for-loop, more expensive (?)
  fwrite((void*)m[0], sizeof(int), size * size, file); 
  fclose(file);
  // Remember to free memory
  free(m[0]);
  free(m);

  // Open file for reading
  file = fopen("matrix_data.dat", "r");
  if(file == NULL){
    printf("Unable to open file\n");
    return -1;
  }
  for(size_t i = 0; i < size; ++i){
    int temp[size];
    fread((void*) temp, sizeof(int), 10, file);
    for(size_t j = 0; j < size; j++){
      printf("%d\t", temp[j]);
      if(temp[j] != j*i){
        printf("Wrong file contents!\n");
        return -1;
      }
    }
    printf("\n");
  }
  fclose(file);
}
