#include "newton.h"
#include <math.h>
#include <stdio.h>

complex_float newton_iter(complex_float (*f)(complex_float *),
                          complex_float (*df)(complex_float *),
                          complex_float *x) {
  complex_float next;

  complex_float f_val = f(x);
  complex_float df_val = df(x);
  next.re = x->re + f_val.re / df_val.re;
  next.im = x->im + f_val.im / df_val.im;

  return next;
}

// All base functions x^n -1, 0<n<10
complex_float poly_2(complex_float *x) {
  complex_float out;
  out.re = x->re * x->re - x->im * x->im - 1;
  out.im = 2 * x->re * x->im;

  return out;
}

complex_float poly_3(complex_float *x) {
  complex_float out;
  out.re = x->re * x->re * x->re - 3 * x->re * x->im * x->im - 1;
  out.im = 3 * x->re * x->re - x->im * x->im * x->im;

  return out;
}
complex_float poly_4(complex_float *x) {
  complex_float out;
  out.re = x->re * x->re * x->re * x->re - 6 * x->re * x->re * x->im * x->im +
           x->im * x->im * x->im * x->im - 1;
  out.im = 4 * x->re * x->re * x->re - 4 * x->re * x->im * x->im * x->im;

  return out;
}
complex_float poly_5(complex_float *x) {
  complex_float out;
  out.re = x->re * x->re * x->re * x->re * x->re -
           10 * x->re * x->re * x->im * x->im * x->im +
           5 * x->re * x->im * x->im * x->im * x->im - 1;
  out.im = 5 * x->re * x->re * x->re * x->re * x->im - 10 * x->re * x->re +
           x->im * x->im * x->im * x->im * x->im;
  return out;
}

complex_float poly_6(complex_float *x) {
  complex_float out;
  out.re = pow(x->re, 6) - 15 * pow(x->re, 4) * pow(x->im, 2) +
           15 * pow(x->re, 2) * pow(x->im, 4) - pow(x->im, 6) - 1;
  out.im = 6 * pow(x->re, 5) * x->im - 20 * pow(x->re, 3) * pow(x->im, 3) +
           6 * x->re * pow(x->im, 5);

  return out;
}

complex_float poly_7(complex_float *x) {
  complex_float out;
  out.re = pow(x->re, 7) - 21 * pow(x->re, 5) * pow(x->im, 2) +
           35 * pow(x->re, 3) * pow(x->im, 4) - 7 * x->re * pow(x->im, 6) - 1;
  out.im = 7 * pow(x->re, 6) * x->im - 35 * pow(x->re, 4) * pow(x->im, 3) +
           21 * pow(x->re, 2) * pow(x->im, 5) - pow(x->im, 7);
  return out;
}

complex_float poly_8(complex_float *x) {
  complex_float out;
  out.re = pow(x->re, 8) - 28 * pow(x->re, 6) * pow(x->im, 2) +
           70 * pow(x->re, 4) * pow(x->im, 4) -
           28 * pow(x->re, 2) * pow(x->im, 6) + pow(x->im, 8) - 1;
  out.im = 8 * pow(x->re, 7) * x->im - 56 * pow(x->re, 5) * pow(x->im, 3) +
           56 * pow(x->re, 3) * pow(x->im, 5) - 8 * x->re * pow(x->im, 7);
  return out;
}
complex_float poly_9(complex_float *x) {
  complex_float out;
  out.re = pow(x->re, 9) - 36 * pow(x->re, 7) * pow(x->im, 2) +
           126 * pow(x->re, 5) * pow(x->im, 4) -
           84 * pow(x->re, 3) * pow(x->im, 6) + 9 * x->re * pow(x->im, 8) - 1;
  out.im = 9 * pow(x->re, 8) * x->im - 84 * pow(x->re, 6) * pow(x->im, 3) +
           126 * pow(x->re, 4) * pow(x->im, 5) -
           36 * pow(x->re, 2) * pow(x->im, 7) + pow(x->im, 9);
  return out;
}

// All derivative functions, n*x^(n-1)
complex_float dpoly_2(complex_float *x) {
  complex_float out;
  out.re = 2 * x->re;
  out.im = 2 * x->im;
  return out;
}

complex_float dpoly_3(complex_float *x) {
  complex_float out = poly_2(x);
  out.re = 3 * (out.re + 1);
  out.im *= 3;
  return out;
}

complex_float dpoly_4(complex_float *x) {
  complex_float out = poly_3(x);
  out.re = 4 * (out.re + 1);
  out.im *= 4;
  return out;
}
complex_float dpoly_5(complex_float *x) {
  complex_float out = poly_4(x);
  out.re = 5 * (out.re + 1);
  out.im *= 5;
  return out;
}
complex_float dpoly_6(complex_float *x) {
  complex_float out = poly_5(x);
  out.re = 6 * (out.re + 1);
  out.im *= 6;
  return out;
}
complex_float dpoly_7(complex_float *x) {
  complex_float out = poly_6(x);
  out.re = 7 * (out.re + 1);
  out.im *= 7;
  return out;
}
complex_float dpoly_8(complex_float *x) {
  complex_float out = poly_7(x);
  out.re = 8 * (out.re + 1);
  out.im *= 8;
  return out;
}
complex_float dpoly_9(complex_float *x) {
  complex_float out = poly_8(x);
  out.re = 9 * (out.re + 1);
  out.im *= 9;
  return out;
}
