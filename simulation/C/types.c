#include <stdlib.h>

#include "types.h"

// some global constants
const double Kb = 1.3806503E-23;      // boltzmann constant
const double Lorentz = 2.45E-8;       // Lorentz number

// function that frees SimRes struct
void free_simres(SimData * data, SimRes * res) {
    // first free all the matrix contents
    for (unsigned i=0; i<data->numberOfIv; ++i)
        free(res->Iv[i]);
    for (unsigned r=0; r<data->numberOfT; ++r)
        free(res->R_w[r]);
    for (unsigned v=0; v<data->numberOfC; ++v)
        free(res->V_c[v]);
    for (unsigned t=0; t<data->numberOfT; ++t) {
        for (unsigned n=0; n<data->N/data->timeskip; ++n)
            free(res->T[t][n]);
        free(res->T[t]);
    }
    free(res->T);
    free(res->Iv);
    free(res->R_w);
    free(res->V_c);
    free(res->dX);
    free(res->rho_norm);
    free(res);

    return;
}

// function that frees SimData struct
void free_simdata(SimData * data) {
    for (unsigned j=0; j<data->groupsOfR; ++j)
        free(data->R[j]);
    for (unsigned j=0; j<data->groupsOfL; ++j)
        free(data->L_x[j]);
    free(data->numberOfR);
    free(data->numberOfL);
    free(data->Iv_b);
    free(data->J);
    free(data->wireLength);
    free(data->wireThickness);
    free(data->wireWidth);
    free(data->L_w);
    free(data->Iv_c0);
    free(data->cor_Iv);
    free(data->Iv_init);
    free(data->C_init);
    free(data->trigger);
    free(data->initHS_l);
    free(data->initHS_T);

    free(data);

    return;
}
