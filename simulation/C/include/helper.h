#ifndef __HELPER_H__
#define __HELPER_H__

void print_vector(double * vec, size_t J);
void print_matrix(double * mat, size_t m);
double max_vector(double * vec, size_t J);
double min_vector(double * vec, size_t J);
double sum_vector(double * vec, size_t J);
double subsum_vector(double * vec, size_t i, size_t j);
int cmp_vector(double * vec, size_t J, double val, double eps);
int fill_vector(double * vec, size_t J, double val);
void print_progress(unsigned n, size_t N);

inline void swap_ptr(void ** one, void ** two) {
    void * tmp = *one;
    *one = *two;
    *two = tmp;

    return;
}

inline void swap_dbl(double * one, double * two) {
    double tmp = *one;
    *one = *two;
    *two = tmp;

    return;
}

#endif
