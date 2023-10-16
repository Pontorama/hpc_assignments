#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int
main(){
	clock_t t;
	double time;
	int n_iters = 10000000;
	double x;
	double res;
	t = clock();
	for(int i = 0; i < n_iters; i++){
		x = 0.82519;
		res = pow(x,10);
	}
	t = clock()-t;
	time = (double)t/CLOCKS_PER_SEC;
	printf("Time for pow: %f s\n", time);

	t = clock();
	for(int i = 0; i < n_iters; i++){
		x = 0.82519;
		res = x*x*x*x*x*x*x*x*x*x;
	}
	t = clock()-t;
	time = (double)t/CLOCKS_PER_SEC;
	printf("Time for manual: %f s\n", time);

}
