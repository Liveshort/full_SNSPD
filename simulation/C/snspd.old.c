#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"
#include "yang.h"
#include "yang_par.h"
#include "waterfall_2s_res.h"
#include "waterfall_3s_res.h"
#include "resistive_summation_3s_res.h"
#include "binary_neural_network.h"
#include "binary_neural_network_xor.h"

// function that coordinates the overall simulation. Data comes into this function from python
//     or another C script or whatever, is processed by the library, and is then returned to
//     the user as a result struct.
// runType is used as follows:
//   - 0: the snspd standard model is simulated. it assumes nothing but the snspd and a load
//            resistor connected after a capacitor
//   - 1: the snspd standard model with a parallel high pass filter is simulated. it assumes an
//            snspd with a parallel resistor and inductor, connected to a load resistor and
//            capacitor.
SimRes * run_snspd_simulation(SimData * data, int runType) {
    printf("    Runtype %d\n\n", runType);

    // set the random seed to something random
    srand((unsigned) time(0));

    // first locally save some important parameters that we will need all the time
    size_t J0, J1, J2;
    J0 = data->J0;
    // put J1 and J2 to zero to suppress "maybe uninitialized" warning from the gcc compiler
    J1 = 0;
    J2 = 0;
    // ...and overwrite if we actually want to use it
    if (data->numberOfT > 1) J1 = data->J1;
    if (data->numberOfT > 2) J2 = data->J2;
    size_t N = data->N;
    size_t NT = data->N/data->timeskip;
    size_t NE = data->N*data->ETratio;
    size_t NTL = data->NTL;

    // create the simulation result struct and allocate sufficient memory
    //   number of time samples N
    //   number of spacial samples J
    SimRes * res = calloc(1, sizeof(SimRes));
    res->runType = data->runType;
    if (NTL > 0 && data->simTL) res->tlType = 1;
    else if (NTL > 0 && !data->simTL) res->tlType = 2;
    else res->tlType = 0;

    res->J = calloc(data->numberOfT, sizeof(size_t));
    res->J[0] = J0;
    if (data->numberOfT > 1) res->J[1] = J1;
    if (data->numberOfT > 2) res->J[2] = J2;
    if (data->runType == 10)
        for (unsigned k=3; k<data->numberOfT; ++k)
            res->J[k] = res->J[k-3];
    if (data->runType == 20 || data->runType == 21)
        for (unsigned k=1; k<data->numberOfT; ++k)
            res->J[k] = res->J[0];
    res->N = N;

    res->numberOfT = data->numberOfT;
    res->numberOfI = data->numberOfI;
    res->numberOfR = data->numberOfR;
    res->numberOfC = data->numberOfC;

    res->timeskip = data->timeskip;
    res->ETratio = data->ETratio;

    res->T = calloc(data->numberOfT, sizeof(double **));
    for (unsigned t=0; t<data->numberOfT; ++t) {
        res->T[t] = calloc(NT, sizeof(double *));
        for (unsigned n=0; n<NT; ++n)
            res->T[t][n] = calloc(res->J[t], sizeof(double));
    }

    res->I = calloc(data->numberOfI, sizeof(double *));
    for (unsigned i=0; i<data->numberOfI; ++i)
        res->I[i] = calloc(NE, sizeof(double));

    res->R = calloc(data->numberOfR, sizeof(double *));
    for (unsigned r=0; r<data->numberOfR; ++r)
        res->R[r] = calloc(NE, sizeof(double));

    res->V_c = calloc(data->numberOfC, sizeof(double *));
    for (unsigned v=0; v<data->numberOfC; ++v)
        res->V_c[v] = calloc(NE, sizeof(double));

    // allocate bias currents
    res->I_b = calloc(res->numberOfT, sizeof(double));

    // calculate delta x and delta t
    res->dX = calloc(res->numberOfT, sizeof(double));
    res->dX[0] = data->wireLength / (J0 - 1);
    if (data->numberOfT > 1) res->dX[1] = data->wireLength_1 / (J1 - 1);
    if (data->numberOfT > 2) res->dX[2] = data->wireLength_2 / (J2 - 1);
    if (data->runType == 10)
        for (unsigned k=3; k<data->numberOfT; ++k)
            res->dX[k] = res->dX[k-3];
    if (data->runType == 20 || data->runType == 21) {
        for (unsigned j=0; j<data->I_bnn; ++j)
            res->dX[j] = data->wl_bnn[0] / (res->J[j] - 1);
        if (data->L_bnn == 0) {
            for (unsigned j=0; j<data->O_bnn; ++j)
                res->dX[data->I_bnn + data->L_bnn*data->M_bnn + j] = data->wl_bnn[1 + data->L_bnn] / (res->J[data->I_bnn + data->L_bnn*data->M_bnn + j] - 1);
        } else if (data->L_bnn == 1) {
            for (unsigned k=0; k<data->M_bnn; ++k)
                res->dX[data->I_bnn + k] = data->wl_bnn[1] / (res->J[data->I_bnn + k] - 1);
            for (unsigned j=0; j<data->O_bnn; ++j)
                res->dX[data->I_bnn + data->L_bnn*data->M_bnn + j] = data->wl_bnn[1 + data->L_bnn] / (res->J[data->I_bnn + data->L_bnn*data->M_bnn + j] - 1);
        } else {
            exit(9);
        }
    }
    res->dt = data->tMax / (N - 1);

    switch(runType) {
        case 0:
            run_yang(res, data, res->dX[0], res->dt, J0, N, NE, NTL);
            break;
        case 1:
            run_yang_parallel(res, data, res->dX[0], res->dt, J0, N, NE, NTL);
            break;
        case 4:
            run_waterfall_2s_res(res, data, res->dX[0], res->dX[1], res->dt, J0, J1, N, NE, NTL);
            break;
        case 6:
            run_waterfall_3s_res(res, data, res->dX, res->dt, res->J, N, NE, NTL);
            break;
        case 10:
            run_resistive_summation_3s_res(res, data, res->dX, res->dt, res->J, N, NE, NTL);
            break;
        case 20:
            run_bnn(res, data, res->dX, res->dt, res->J, N, NE, NTL);
            break;
        case 21:
            run_bnn_xor(res, data, res->dX, res->dt, res->J, N, NE, NTL);
            break;
        default:
            printf("    Unknown runtype %d...\nExiting with error 1 (wrong runtype)...\n", runType);
            exit(1);
    }

    res->exitValue = 0;

    return res;
}
