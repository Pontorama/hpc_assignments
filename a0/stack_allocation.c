#include <stdio.h>

int main(){
  int size = 10;
  int as[size];
  for(size_t i = 0; i < size; ++i){
    as[i] = 0;
  }
  // Segfault if size is greater than int size 2^32(?)
  printf("%d\n", as[0]);
}
