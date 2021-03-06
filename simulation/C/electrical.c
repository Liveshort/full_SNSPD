#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "omp.h"

// advances the electric model one timestep using BLAS matrix algebra.
// you need the LAPACK for this, or work out the matrix logic on your own (not recommended)
int advance_time_electrical(unsigned n, size_t numberOfIv, size_t numberOfC, double ** Iv, double ** V_c, double * XW, double ** XX, double * Y, double ** R_w, double ** R, double * Iv_b) {
    // set up matrix and vector for Ax = b calculation
    lapack_int m, mrhs, info;

    m = 3;
    mrhs = 1;

    double * A = calloc(m*m, sizeof(double));
    double ** Am = calloc(m, sizeof(double *));
    double * b = calloc(m*mrhs, sizeof(double));
    lapack_int * ipiv = calloc(m, sizeof(lapack_int));

    for (int j=0; j<m; ++j)
        Am[j] = &A[j*m];

    Am[0][0] = -(- XW[0] - R_w[0][n] + - R[0][0]);
    Am[0][2] = -1;
    Am[0][1] = - R[1][0];
    
    Am[1][0] = 1;
    Am[1][1] = 1;
    
    Am[2][1] = -Y[0];
    Am[2][2] = 1;
    
    b[0] = Iv[0][n-1]*(-(- XW[0] + R_w[0][n-1] + R[0][0])) + V_c[0][n-1] + Iv[1][n-1]*(R[1][0]);
    b[1] = Iv_b[0];
    b[2] = Y[0]*Iv[1][n-1] + V_c[0][n-1];

    // print_matrix(A, m);
    // print_vector(b, m);
    // exit(99);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, m, mrhs, A, m, ipiv, b, mrhs);

    if (info != 0) {
        printf("Error encountered in matrix calculation of electrical model (info = %d)...\nExiting with code 3.\n", info);
        exit(3);
    }

    for (unsigned j=0; j<numberOfIv; ++j)
        Iv[j][n] = b[j];
    for (unsigned j=0; j<numberOfC; ++j)
        V_c[j][n] = b[j + numberOfIv];

    free(Am);
    free(A);
    free(b);
    free(ipiv);

    //if (n>100)
    //    exit(99);

    return 0;
}
