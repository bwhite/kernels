import numpy as np
cimport numpy as np


cdef extern from "kernels_aux.h":
    void histogram_intersection_fast(double *x, double *y, int rows_x, int rows_y, int cols, double *out)
    void linear_fast(double *x, double *y, int rows_x, int rows_y, int cols, double *out)
    void chi_square_fast(double *x, double *y, int rows_x, int rows_y, int cols, double *out)
    void polynomial_fast(double *x, double *y, int rows_x, int rows_y, int cols, double *out, double alpha, double c, int degree)


cpdef histogram_intersection(np.ndarray[np.float64_t, ndim=2, mode='c'] x, np.ndarray[np.float64_t, ndim=2, mode='c'] y):
    cdef np.ndarray out = np.zeros((x.shape[0], y.shape[0]))
    histogram_intersection_fast(<double *>x.data, <double *>y.data, x.shape[0], y.shape[0], x.shape[1], <double *>out.data)
    return out


cpdef linear(np.ndarray[np.float64_t, ndim=2, mode='c'] x, np.ndarray[np.float64_t, ndim=2, mode='c'] y):
    cdef np.ndarray out = np.zeros((x.shape[0], y.shape[0]))
    linear_fast(<double *>x.data, <double *>y.data, x.shape[0], y.shape[0], x.shape[1], <double *>out.data)
    return out


cpdef polynomial(np.ndarray[np.float64_t, ndim=2, mode='c'] x, np.ndarray[np.float64_t, ndim=2, mode='c'] y, alpha, c, degree):
    cdef np.ndarray out = np.zeros((x.shape[0], y.shape[0]))
    polynomial_fast(<double *>x.data, <double *>y.data, x.shape[0], y.shape[0], x.shape[1], <double *>out.data, alpha, c, degree)
    return out


cpdef chi_square(np.ndarray[np.float64_t, ndim=2, mode='c'] x, np.ndarray[np.float64_t, ndim=2, mode='c'] y):
    cdef np.ndarray out = np.zeros((x.shape[0], y.shape[0]))
    chi_square_fast(<double *>x.data, <double *>y.data, x.shape[0], y.shape[0], x.shape[1], <double *>out.data)
    return out
