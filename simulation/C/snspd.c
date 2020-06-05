#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "electrical.h"
#include "transmission.h"
#include "omp.h"

// print parameters
int print_parameters(SimData * data) {
    puts("    PARAMETERS:");
    printf("    Number of (T, Iv, C): %zu %zu %zu\n", data->numberOfT, data->numberOfIv, data->numberOfC);
    printf("    Resistor groups: ");
    for (unsigned j=0; j<data->groupsOfR; ++j) {
        if (j > 0) printf("                     ");
        printf("[ ");
        for (unsigned k=0; k<data->numberOfR[j]; ++k)
            printf("%4.2f ", data->R[j][k]);
        printf("]\n");
    }
    printf("    Inductor groups: ");
    for (unsigned j=0; j<data->groupsOfL; ++j) {
        if (j > 0) printf("                     ");
        printf("[ ");
        for (unsigned k=0; k<data->numberOfL[j]; ++k)
            printf("%7.1e ", data->L_x[j][k]);
        printf("]\n");
    }
    printf("    Capacitors: [ ");
    for (unsigned j=0; j<data->numberOfC; ++j)
        printf("%7.1e ", data->C[j]);
    puts("]");
    printf("    Timing (N, tMax, timeskip, ETratio, opt): %zu %7.1e %zu %zu %d\n", data->N, data->tMax, data->timeskip, data->ETratio, data->allowOpt);
    printf("    Thermal (c_p, c_e, alpha, T_ref): %7.1e %7.1e %7.1e %4.2e\n", data->c_p, data->c_e, data->alpha, data->T_ref);
    printf("    Thermal (T_sub, T_sub_eps, T_c): %4.2f %7.1e %4.2f\n", data->T_sub, data->T_sub_eps, data->T_c);
    printf("    Bias currents: [ ");
    for (unsigned j=0; j<data->numberOfIv_b; ++j)
        printf("%7.1e ", data->Iv_b[j]);
    puts("]");
    printf("    Impurities (offset spread): %7.1e %7.1e\n", data->impurityOffset, data->impuritySpread);
    printf("    Wire sptl smpls: [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7zu ", data->J[j]);
    puts("]");
    printf("    Wirelengths    : [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7.1e ", data->wireLength[j]);
    puts("]");
    printf("    Wirethicknesses: [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7.1e ", data->wireThickness[j]);
    puts("]");
    printf("    Wirewidths     : [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7.1e ", data->wireWidth[j]);
    puts("]");
    printf("    Wire inductancs: [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7.1e ", data->L_w[j]);
    puts("]");
    printf("    Wire crit currs: [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%7.1e ", data->Iv_c0[j]);
    puts("]");
    printf("    Init currs: [ ");
    for (unsigned j=0; j<data->numberOfIv; ++j)
        printf("%8.2e ", data->Iv_init[j]);
    puts("]");
    printf("    Init volts C: [ ");
    for (unsigned j=0; j<data->numberOfC; ++j)
        printf("%8.2e ", data->V_c_init[j]);
    puts("]");
    printf("    Triggers: [ ");
    for (unsigned j=0; j<data->numberOfT; ++j)
        printf("%d ", data->trigger[j]);
    puts("]");
    printf("    Trns line (type, NTL, VF, LTL): %d %zu %4.2f %4.2f\n\n", data->simTL, data->NTL, data->VF, data->LTL);

    return 0;
}

// main function that runs the simulation
int run_sim(SimRes * res, SimData * data, double * dX, double dt, size_t * J, size_t N, size_t NE, size_t NTL) {
    puts("\n    SETTING UP SIMULATION...\n");

    // print the parameters of the sim
    print_parameters(data);

    // first locally save some important parameters that we will need all the time
    double *** T = calloc(data->numberOfT, sizeof(double **));
    for (unsigned j=0; j<data->numberOfT; ++j)
        T[j] = res->T[j];

    double ** Iv = calloc(data->numberOfIv, sizeof(double *));
    for (unsigned j=0; j<data->numberOfIv; ++j)
        Iv[j] = res->Iv[j];

    double ** R_w = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j)
        R_w[j] = res->R_w[j];

    double ** V_c = calloc(data->numberOfC, sizeof(double *));
    for (unsigned j=0; j<data->numberOfC; ++j)
        V_c[j] = res->V_c[j];

    // set up vectors to temporarily save the current and next T
    // this saves a lot of memory for larger simulations
    double ** T_stash = calloc(2*data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        T_stash[2*j] = calloc(J[j], sizeof(double));
        T_stash[2*j+1] = calloc(J[j], sizeof(double));
    }
    double ** T_prev = calloc(data->numberOfT, sizeof(double *));
    double ** T_curr = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        T_prev[j] = T[j][0];
        T_curr[j] = T[j][1];
        // set up initial thermal values at t = 1, add a steady state time step at t = 0
        fill_vector(T_prev[j], J[j], data->T_sub);
        fill_vector(T_curr[j], J[j], data->T_sub);
    }

    // determine halfway point and set up a beginning hotspot at t = 1
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (data->trigger[j] > 0) {
            unsigned halfway = J[j]/2;
            unsigned initHS_segs = (unsigned) (data->HS_l_init[j]/dX[j]) + 1;
            // check if there is a nonzero number of segments
            if (initHS_segs < 2) {
                puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
                exit(2);
            }
            for (unsigned k=halfway - initHS_segs/2; k<halfway + initHS_segs/2; ++k)
                T_curr[j][k] = data->HS_T_init[j];
        }

    // set up vectors for the critical currents to simulate nanowire impurities
    double ** Iv_c = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        Iv_c[j] = calloc(J[j], sizeof(double));

        for (unsigned k=0; k<J[j]; ++k)
            Iv_c[j][k] = data->Iv_c0[j] + data->impurityOffset + data->impuritySpread*(((double) rand() / (double) RAND_MAX));

        Iv_c[j][J[j]/2] = data->Iv_c0[j];
    }

    // determine initial conditions
    for (unsigned n=0; n<=data->timeskip; ++n) {
        for (unsigned j=0; j<data->numberOfIv; ++j)
            Iv[j][n] = data->Iv_init[j];
        for (unsigned j=0; j<data->numberOfC; ++j)
            V_c[j][n] = data->V_c_init[j];
    }

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref/data->T_c)*(data->T_ref/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref, 3));

    printf("    DERIVED THERMAL PROPERTIES:\n    Delta: %e\n    A:     %e\n    gamma: %e\n    B:     %e\n\n", DeltaRef, A, gamma, B);

    // define the resistance of a segment of wire in the normal state
    double * R_seg = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; j++)
        R_seg[j] = res->rho_norm[j]*dX[j]/(data->wireWidth[j]*data->wireThickness[j]);
    // declare the nanowire resistance and current density
    for (unsigned n=0; n<=data->timeskip; ++n)
        for (unsigned j=0; j<data->numberOfT; ++j)
            R_w[j][n] = 0;
    double * currentDensity_w = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        currentDensity_w[j] = 0;

    // allocate space for the state and temperature dependent variables for each time step
    double ** alpha_n = calloc(data->numberOfT, sizeof(double *));
    double ** kappa_n = calloc(data->numberOfT, sizeof(double *));
    double ** c_n = calloc(data->numberOfT, sizeof(double *));
    double ** rho_seg_n = calloc(data->numberOfT, sizeof(double *));
    double ** R_seg_n = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        alpha_n[j] = calloc(J[j], sizeof(double));
        kappa_n[j] = calloc(J[j], sizeof(double));
        c_n[j] = calloc(J[j], sizeof(double));
        rho_seg_n[j] = calloc(J[j], sizeof(double));
        R_seg_n[j] = calloc(J[j], sizeof(double));
    }

    // set up two characteristic numbers for the electrical calculations
    double * XW = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        XW[j] = (2*data->L_w[j])/dt;
    double ** XX = calloc(data->groupsOfL, sizeof(double *));
    for (unsigned j=0; j<data->groupsOfL; ++j) {
        XX[j] = calloc(data->numberOfL[j], sizeof(double));
        for (unsigned k=0; k<data->numberOfL[j]; ++k)
            XX[j][k] = (2*data->L_x[j][k])/dt;
    }
    double * Y = calloc(data->numberOfC, sizeof(double));
    for (unsigned j=0; j<data->numberOfC; ++j)
        Y[j] = dt/(2*data->C[j]);

    puts("\n    STARTING SIMULATION...\n");

    // main time loop
    for (unsigned n=data->timeskip+1; n<NE; ++n) {
        // print progress
        print_progress(n, NE);

        // advance the thermal model to the next time step after the initial step
        if (n > data->timeskip+1 && n < N)
            #pragma omp parallel for
            for (unsigned j=0; j<data->numberOfT; ++j)
                advance_time_thermal(T_prev[j], T_curr[j], J[j], data->T_sub, alpha_n[j], c_n[j], rho_seg_n[j], kappa_n[j], data->wireThickness[j], currentDensity_w[j], dt, dX[j]);

        if (n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            #pragma omp parallel for
            for (unsigned j=0; j<data->numberOfT; ++j) {
                update_thermal_values(alpha_n[j], kappa_n[j], c_n[j], rho_seg_n[j], R_seg_n[j], T_curr[j], J[j], A, B, gamma, data->T_c, Iv[data->cor_Iv[j]][n-1], Iv_c[j], res->rho_norm[j], data->c_p, data->T_ref, R_seg[j]);
            }
            // update the current nanowire resistance
            for (unsigned j=0; j<data->numberOfT; ++j)
                R_w[j][n] = sum_vector(R_seg_n[j], J[j]);
        } else {
            for (unsigned j=0; j<data->numberOfT; ++j)
                R_w[j][n] = 0;
        }

        // update the current density through the nanowire
        for (unsigned j=0; j<data->numberOfT; ++j)
                currentDensity_w[j] = Iv[data->cor_Iv[j]][n-1]/(data->wireWidth[j]*data->wireThickness[j]);
        // update the electric values
        advance_time_electrical(n, data->numberOfIv, data->numberOfC, Iv, V_c, XW, XX, Y, R_w, data->R, data->Iv_b);

        // shuffle the T pointers around so the old and new timestep don't point to the same array
        for (unsigned j=0; j<data->numberOfT; ++j)
            T_prev[j] = T_curr[j];
        if (n % data->timeskip == 0 && n < N) {
            for (unsigned j=0; j<data->numberOfT; ++j)
                T_curr[j] = T[j][n/data->timeskip];
        } else {
            for (unsigned j=0; j<data->numberOfT; ++j) {
                if (T_prev[j] == T_stash[2*j])
                    T_curr[j] = T_stash[2*j+1];
                else
                    T_curr[j] = T_stash[2*j];
            }
        }
    }

    // free allocated space
    for (unsigned j=0; j<2*data->numberOfT; ++j)
        free(T_stash[j]);
    free(T_stash);

    free(T_prev);
    free(T_curr);

    for (unsigned j=0; j<data->numberOfT; ++j) {
        free(alpha_n[j]);
        free(kappa_n[j]);
        free(c_n[j]);
        free(rho_seg_n[j]);
        free(R_seg_n[j]);
        free(Iv_c[j]);
    }
    free(alpha_n);
    free(kappa_n);
    free(c_n);
    free(rho_seg_n);
    free(R_seg_n);
    free(Iv_c);

    free(T);
    free(Iv);
    free(R_w);
    free(V_c);

    free(R_seg);
    free(currentDensity_w);

    free(XW);
    for (unsigned j=0; j<data->groupsOfL; ++j)
        free(XX[j]);
    free(XX);
    free(Y);

    // print result
    puts("\n\n    SIMULATION COMPLETED.");
    return 0;
}
