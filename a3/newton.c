#include <threads.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>

double* newton(double real, double imag, unsigned short degree){

  //List of the roots where the first column is real and the second imaginary
  static double root_list[degree][2];
  for(unsigned short ix = 0; ix < degree; ix++){
    root_list[ix][1] = cos(2*M_PI*ix/degree);
    root_list[ix][2] = sin(2*M_PI*ix/degree);
  }
  // z_k
  double complex z_k = CMPLX(real, imag);
  // z_k^d-1
  double complex z = CMPLX(1, 0);
 // f(z)
 double complex f = z;
 // f'(z)
 double complex fp = z;
 unsigned short iter = 0;
 unsigned short kx;
 while(iter < 128){
   if( creal(z_k) < -1000000000 || 1000000000 < creal(z_k) || cimag < -1000000000 || 1000000000 < cimag(z_k)) break;  
   for( kx = 0; kx < degree; kx++){
     if( (root_list[kx][1]-creal(z_k))^2 + (root_list[kx][2]-cimag(z_k))^2 <= 0.000001) break; //0.001^2
   }
     
     for(unsigned short jx = 0; jx < degree -1; jx++){
      z = z_k * z;
     }
   f = z_k * z - 1;
   fp = degree*z;
   z_k = z_k-f/fp;
   iter++
 }
 if(iter == 128)
   return 0;
 return root_list[kx];
}
