#ifndef KERNELS_AUX_H
#define KERNELS_AUX_H
void histogram_intersection_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out);
void chi_square_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out);
void linear_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out);
void polynomial_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out, double alpha, double c, int degree);
#endif
