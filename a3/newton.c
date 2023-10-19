#include <threads.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

type main(double real, double imag, unsigned short degree){

  //List of the roots where the first column is real and the second imaginary
  double root_list[degree][2];
  for(unsigned short ix = 0; ix < degree; ix++){
    root_list[ix][1] = cos(2*M_PI*ix/degree);
    root_list[ix][2] = sin(2*M_PI*ix/degree);
  }
 
double complex z = CMPLX(real, imag);
 double complex z_k = CMPLX(1, 0);
 // f(z)
 double complex f = z;
 // f'(z)
 double complex fp = z;
 bool run = true;
 unsigned short iter = 0;
 while(iter < 128){
   if( creal(z_k) < -2 || 2 < creal(z_k) || cimag < -2 || 2 < cimag(z_k)) break;  
   for(unsigned short kx = 0; kx < degree; kx++){
     if( (root[kx][1]-creal(z_k))^2 + (root[kx][2]-cimag(z_k))^2 <= 0.000001) break; //0.001^2
   }
     
     for(unsigned short jx = 0; jx < degree -1; jx++){
      z_k = z_k * z;
     }
   f = z_k * z - 1;
   fp = (degree-1)*z_k;
   z_k = z_k-f/fp;
   iter++
 }


 switch(degree){
 case 1: degree == 1
   
 case 2: degree == 2
    root(
    


 default:
   fprintf(stderr, "Unexpected degree\n");
   exit(1);

 }
   // x_k = x-(x^d-1)/(x^(d-1)) = (x-1)/d + 1/(d*x^(d-1))  
     return x_k;
     }
