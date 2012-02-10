#include "kernels_aux.h"

static inline double powi(double base, int times) {
    // Adapted from libsvm 3-1 (BSD license)
    double tmp = base, ret = 1.0;
    int t;
    for (t = times; t > 0; t /= 2) {
        if (t % 2 == 1)
            ret *= tmp;
        tmp = tmp * tmp;
    }
    return ret;
}

// Same as the following in python:  return np.array([[np.sum(np.min([x0, y0], 0)) for x0 in x] for y0 in y])
void histogram_intersection_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out) {
    int i, j, k;
    double cur_sum, *cur_x, *cur_y;
    for (i = 0; i < rows_x; ++i)
        for (j = 0; j < rows_y; ++j) {
            cur_sum = 0.;
            cur_x = x + i * cols;
            cur_y = y + j * cols;
            for (k = 0; k < cols; ++k)
                cur_sum += cur_x[k] < cur_y[k] ? cur_x[k] : cur_y[k];
            out[i * rows_y + j] = cur_sum;
        }
}


void linear_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out) {
    int i, j, k;
    double cur_sum, *cur_x, *cur_y;
    for (i = 0; i < rows_x; ++i)
        for (j = 0; j < rows_y; ++j) {
            cur_sum = 0.;
            cur_x = x + i * cols;
            cur_y = y + j * cols;
            for (k = 0; k < cols; ++k)
                cur_sum += cur_x[k] * cur_y[k];
            out[i * rows_y + j] = cur_sum;
        }
}

void polynomial_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out, double alpha, double c, int degree) {
    int i, j, k;
    double cur_sum, *cur_x, *cur_y;
    for (i = 0; i < rows_x; ++i)
        for (j = 0; j < rows_y; ++j) {
            cur_sum = 0.;
            cur_x = x + i * cols;
            cur_y = y + j * cols;
            for (k = 0; k < cols; ++k)
                cur_sum += cur_x[k] * cur_y[k];
            out[i * rows_y + j] = powi(alpha * cur_sum + c, degree);
        }
}

void chi_square_fast(double *x, double *y, const int rows_x, const int rows_y, const int cols, double *out) {
    int i, j, k;
    double cur_sum, *cur_x, *cur_y, temp_diff, temp_sum;
    for (i = 0; i < rows_x; ++i)
        for (j = 0; j < rows_y; ++j) {
            cur_sum = 0.;
            cur_x = x + i * cols;
            cur_y = y + j * cols;
            for (k = 0; k < cols; ++k) {
                temp_diff = cur_x[k] - cur_y[k];
                temp_sum = cur_x[k] + cur_y[k];
                cur_sum += temp_diff * temp_diff / temp_sum;
            }
            out[i * rows_y + j] = 1. - 2 * cur_sum;
        }
}
