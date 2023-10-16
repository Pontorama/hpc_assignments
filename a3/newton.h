#ifndef Newton_h
#define Newton_h

typedef struct {
  float re;
  float im;
} complex_float;

complex_float newton_iter(complex_float (*f)(complex_float *),
                          complex_float (*df)(complex_float *),
                          complex_float *x);

// Functions
complex_float poly_2(complex_float *x);
complex_float poly_3(complex_float *x);
complex_float poly_4(complex_float *x);
complex_float poly_5(complex_float *x);
complex_float poly_6(complex_float *x);
complex_float poly_7(complex_float *x);
complex_float poly_8(complex_float *x);
complex_float poly_9(complex_float *x);
// Derivatives
complex_float dpoly_2(complex_float *x);
complex_float dpoly_3(complex_float *x);
complex_float dpoly_4(complex_float *x);
complex_float dpoly_5(complex_float *x);
complex_float dpoly_6(complex_float *x);
complex_float dpoly_7(complex_float *x);
complex_float dpoly_8(complex_float *x);
complex_float dpoly_9(complex_float *x);
#endif
