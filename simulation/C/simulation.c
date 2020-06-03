#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include <sys/sysinfo.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"
#include "linalg.h"

int collect_data(char * inputPath, SimData * data) {
    // open file
    FILE * fp;

    if ((fp = fopen(inputPath, "r")) == NULL) {
        printf("Something went wrong trying to open the file \"%s\"...\nExiting...\n", inputPath);
        exit(5);
    }

    // create throwaway variable for text from the setup file we don't use
    char dump[2000];

        // skip first bogus line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan sim memory parameters
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfT, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfIv, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfC, dump) < 1) exit(6);

    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->groupsOfR, dump) < 1) exit(6);
    if (data->groupsOfR > 0) {
        data->numberOfR = calloc(data->groupsOfR, sizeof(size_t));
        for (unsigned j=0; j<data->groupsOfR; ++j)
            if (fscanf(fp, "%zu;", &data->numberOfR[j]) < 1) exit(6);
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);
    }

    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->groupsOfL, dump) < 1) exit(6);
    if (data->groupsOfL > 0) {
        data->numberOfL = calloc(data->groupsOfL, sizeof(size_t));
        for (unsigned j=0; j<data->groupsOfL; ++j)
            if (fscanf(fp, "%zu;", &data->numberOfL[j]) < 1) exit(6);
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);
    }

    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfIv_b, dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan sim timing parameters
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->N, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->tMax, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->timeskip, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->ETratio, dump) < 1) exit(6);
    if (fscanf(fp, "%d;%2000[^\n]\n", &data->allowOpt, dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan thermal parameters
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_p, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_e, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->alpha, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_ref, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_gamma, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_c, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub_eps, dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan bias current values
    data->Iv_b = calloc(data->numberOfIv_b, sizeof(double));
    for (unsigned j=0; j<data->numberOfIv_b; ++j)
        if (fscanf(fp, "%lf;", &data->Iv_b[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan nanowire parameters
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impurityOffset, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impuritySpread, dump) < 1) exit(6);

    data->J = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%zu;", &data->J[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->wireLength = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->wireLength[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->wireThickness = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->wireThickness[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->wireWidth = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->wireWidth[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->L_w = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->L_w[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->Iv_c0 = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->Iv_c0[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->cor_Iv = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%u;", &data->cor_Iv[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan other component parameters
    data->R = calloc(data->groupsOfR, sizeof(double *));
    for (unsigned j=0; j<data->groupsOfR; ++j) {
        data->R[j] = calloc(data->numberOfR[j], sizeof(double));
        for (unsigned k=0; k<data->numberOfR[j]; ++k)
            if (fscanf(fp, "%lf;", &data->R[j][k]) < 1) exit(6);
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);
    }

    data->L_x = calloc(data->groupsOfL, sizeof(double *));
    for (unsigned j=0; j<data->groupsOfL; ++j) {
        data->L_x[j] = calloc(data->numberOfL[j], sizeof(double));
        for (unsigned k=0; k<data->numberOfL[j]; ++k)
            if (fscanf(fp, "%lf;", &data->L_x[j][k]) < 1) exit(6);
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);
    }

    data->C = calloc(data->numberOfC, sizeof(double));
    for (unsigned j=0; j<data->numberOfC; ++j)
        if (fscanf(fp, "%lf;", &data->C[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan initial conditions
    data->Iv_init = calloc(data->numberOfIv, sizeof(double));
    for (unsigned j=0; j<data->numberOfIv; ++j)
        if (fscanf(fp, "%lf;", &data->Iv_init[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->C_init = calloc(data->numberOfC, sizeof(double));
    for (unsigned j=0; j<data->numberOfC; ++j)
        if (fscanf(fp, "%lf;", &data->C_init[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->trigger = calloc(data->numberOfT, sizeof(int));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%d;", &data->trigger[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->initHS_l = calloc(data->numberOfT, sizeof(int));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->initHS_l[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->initHS_T = calloc(data->numberOfT, sizeof(int));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->initHS_T[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        // skip line
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    // scan transmission line parameters
    if (fscanf(fp, "%d;%2000[^\n]\n", &data->simTL, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->NTL, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->VF, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->LTL, dump) < 1) exit(6);

    fclose(fp);

    return 0;
}

SimRes * snspd_simulation(SimData * data) {
    // set the random seed to something random
    srand((unsigned) time(0));

    size_t N = data->N;
    size_t NT = data->N/data->timeskip;
    size_t NE = data->N*data->ETratio;
    size_t NTL = data->NTL;

    // create the simulation result struct and allocate sufficient memory
    SimRes * res = calloc(1, sizeof(SimRes));

    res->T = calloc(data->numberOfT, sizeof(double **));
    for (unsigned t=0; t<data->numberOfT; ++t) {
        res->T[t] = calloc(NT, sizeof(double *));
        for (unsigned n=0; n<NT; ++n)
            res->T[t][n] = calloc(data->J[t], sizeof(double));
    }

    res->Iv = calloc(data->numberOfIv, sizeof(double *));
    for (unsigned i=0; i<data->numberOfIv; ++i)
        res->Iv[i] = calloc(NE, sizeof(double));

    res->R_w = calloc(data->numberOfT, sizeof(double *));
    for (unsigned r=0; r<data->numberOfT; ++r)
        res->R_w[r] = calloc(NE, sizeof(double));

    res->V_c = calloc(data->numberOfC, sizeof(double *));
    for (unsigned v=0; v<data->numberOfC; ++v)
        res->V_c[v] = calloc(NE, sizeof(double));

    // calculate delta x and delta t
    res->dX = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        res->dX[j] = data->wireLength[j] / (data->J[j] - 1);
    res->dt = data->tMax / (N - 1);

    // calculate rho_norm
    res->rho_norm = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        res->rho_norm[j] = data->R_gamma*data->wireThickness[j];

    run_sim(res, data, res->dX, res->dt, data->J, N, NE, NTL);

    return res;
}

// int write_results(char * outputPath, FILE * fp, SimData * data, SimRes * res) {
//     // generate correct filenames
//     const char Tbin[] = "T.bin";
//     const char Ibin[] = "I.bin";
//     const char Rbin[] = "R.bin";
//     const char Cbin[] = "V_c.bin";
//     const char paramInfo[] = "param.info";
//     char *TFilename = calloc(strlen(outputPath) + strlen(Tbin) + 1, sizeof(char));
//     char *IFilename = calloc(strlen(outputPath) + strlen(Ibin) + 1, sizeof(char));
//     char *RFilename = calloc(strlen(outputPath) + strlen(Rbin) + 1, sizeof(char));
//     char *CFilename = calloc(strlen(outputPath) + strlen(Cbin) + 1, sizeof(char));
//     char *paramInfoFilename = calloc(strlen(outputPath) + strlen(paramInfo) + 1, sizeof(char));
//     snprintf(TFilename, strlen(outputPath) + strlen(Tbin) + 1, "%s%s", outputPath, Tbin);
//     snprintf(IFilename, strlen(outputPath) + strlen(Ibin) + 1, "%s%s", outputPath, Ibin);
//     snprintf(RFilename, strlen(outputPath) + strlen(Rbin) + 1, "%s%s", outputPath, Rbin);
//     snprintf(CFilename, strlen(outputPath) + strlen(Cbin) + 1, "%s%s", outputPath, Cbin);
//     snprintf(paramInfoFilename, strlen(outputPath) + strlen(paramInfo) + 1, "%s%s", outputPath, paramInfo);
//
//     // write data to binary files
//     FILE * fp;
//     fp = fopen(TFilename, "wb");
//     for (unsigned q=0; q<res->numberOfT; ++q) {
//         for (unsigned n=0; n<res->N/res->timeskip; ++n)
//         fwrite(res->T[q][n], sizeof(double), res->J[q], fp);
//     }
//     fclose(fp);
//
//     fp = fopen(IFilename, "wb");
//     for (unsigned q=0; q<res->numberOfI; ++q) {
//         for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
//             fwrite(&res->I[q][n], sizeof(double), 1, fp);
//     }
//     fclose(fp);
//
//     fp = fopen(RFilename, "wb");
//     for (unsigned q=0; q<res->numberOfR; ++q) {
//         for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
//             fwrite(&res->R[q][n], sizeof(double), 1, fp);
//     }
//     fclose(fp);
//
//     fp = fopen(CFilename, "wb");
//     for (unsigned q=0; q<res->numberOfC; ++q) {
//         for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
//             fwrite(&res->V_c[q][n], sizeof(double), 1, fp);
//     }
//     fclose(fp);
//
//     // write simulation parameters to file for readout
//     fp = fopen(paramInfoFilename, "w");
//     fprintf(fp, "%40s; %d\n", "runtype of the simulation", res->runType);
//     fprintf(fp, "%40s; %d\n", "transmission simulation type", res->tlType);
//     fprintf(fp, "%40s; ", "J (# of spatial elements)");
//     for (unsigned i=0; i<res->numberOfT; ++i) {
//         if (i > 0) fprintf(fp, "; ");
//         fprintf(fp, "%zu", res->J[i]);
//     }
//     fprintf(fp, "\n%40s; %zu\n", "N (# of temporal elements)", res->N);
//     fprintf(fp, "%40s; %zu\n", "timeskip factor", res->timeskip);
//     fprintf(fp, "%40s; %zu\n", "electrical / thermal time ratio", res->ETratio);
//     fprintf(fp, "%40s; %zu\n", "# of nanowires", res->numberOfT);
//     fprintf(fp, "%40s; %zu\n", "# of currents", res->numberOfI);
//     fprintf(fp, "%40s; %zu\n", "# of resistances", res->numberOfR);
//     fprintf(fp, "%40s; %zu\n", "# of capacitor voltages", res->numberOfC);
//     fprintf(fp, "%40s; ", "I_b [A]");
//     for (unsigned i=0; i<res->numberOfT; ++i) {
//         if (i > 0) fprintf(fp, "; ");
//         fprintf(fp, "%8.6e", res->I_b[i]);
//     }
//     fprintf(fp, "\n%40s; ", "dX [m]");
//     for (unsigned i=0; i<res->numberOfT; ++i) {
//         if (i > 0) fprintf(fp, "; ");
//         fprintf(fp, "%8.6e", res->dX[i]);
//     }
//     fprintf(fp, "\n%40s; %8.6e\n", "dt [s]", res->dt);
//     fclose(fp);
//
//     free(TFilename);
//     free(IFilename);
//     free(RFilename);
//     free(CFilename);
//     free(paramInfoFilename);
// }

int main(int argc, char * argv[]) {
    puts("########################################################################################\n################################### SNSPD SIMULATION ###################################\n########################################################################################\n");

    // check and open setup file
    if (argc != 3) {
        puts("Wrong number of arguments...\nUsage: ./simulation inputFile outputFolder\nExiting...");
        exit(4);
    }

    printf("    %d processors detected. %d processors available.\n", get_nprocs_conf(), get_nprocs());
    printf("    OMP_NUM_THREADS set to %s, using that amount of threads.\n\n", getenv("OMP_NUM_THREADS"));

    printf("    Input file:       %s\n", argv[1]);
    printf("    Output folder:    %s\n\n", argv[2]);

    // put in the experiment data
    SimData * data = calloc(1, sizeof(SimData));
    collect_data(argv[1], data);

    // check if some of the parameters make sense
    if (data->N % data->timeskip != 0 || data->N <= 0) {
        printf("Error: invalid N (%ld), not larger than 0 or N is not an integer multiple of timeskip (%ld)...\nExiting with code 8.\n", data->N, data->timeskip);
        exit(8);
    }

    if (data->timeskip % 10 != 0 || data->timeskip > 100000) {
        printf("Error: timeskip is not a multiple of 10...\nExiting with code 9.\n");
        exit(9);
    }

    // run simulation
    SimRes * res = snspd_simulation(data);

    puts("    Writing files...");
    // write_results(argv[2], fp, data, res);
    puts("    Done.\n");

    free_simres(data, res);
    free_simdata(data);

    puts("########################################################################################\n#################################### END SIMULATION ####################################\n########################################################################################");

    exit(0);
}
