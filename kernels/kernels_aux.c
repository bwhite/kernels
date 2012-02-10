#include "kernels_aux.h"

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
