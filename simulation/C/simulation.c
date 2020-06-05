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

    data->V_c_init = calloc(data->numberOfC, sizeof(double));
    for (unsigned j=0; j<data->numberOfC; ++j)
        if (fscanf(fp, "%lf;", &data->V_c_init[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->trigger = calloc(data->numberOfT, sizeof(int));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%d;", &data->trigger[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->HS_l_init = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->HS_l_init[j]) < 1) exit(6);
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    data->HS_T_init = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        if (fscanf(fp, "%lf;", &data->HS_T_init[j]) < 1) exit(6);
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

int write_results(char * outputPath, SimData * data, SimRes * res) {
    // generate correct filenames
    const char Tbin[] = "T.bin";
    const char Ibin[] = "Iv.bin";
    const char Rbin[] = "R_w.bin";
    const char Cbin[] = "V_c.bin";
    const char paramInfo[] = "param.info";
    char *TFilename = calloc(strlen(outputPath) + strlen(Tbin) + 1, sizeof(char));
    char *IFilename = calloc(strlen(outputPath) + strlen(Ibin) + 1, sizeof(char));
    char *RFilename = calloc(strlen(outputPath) + strlen(Rbin) + 1, sizeof(char));
    char *CFilename = calloc(strlen(outputPath) + strlen(Cbin) + 1, sizeof(char));
    char *paramInfoFilename = calloc(strlen(outputPath) + strlen(paramInfo) + 1, sizeof(char));
    snprintf(TFilename, strlen(outputPath) + strlen(Tbin) + 1, "%s%s", outputPath, Tbin);
    snprintf(IFilename, strlen(outputPath) + strlen(Ibin) + 1, "%s%s", outputPath, Ibin);
    snprintf(RFilename, strlen(outputPath) + strlen(Rbin) + 1, "%s%s", outputPath, Rbin);
    snprintf(CFilename, strlen(outputPath) + strlen(Cbin) + 1, "%s%s", outputPath, Cbin);
    snprintf(paramInfoFilename, strlen(outputPath) + strlen(paramInfo) + 1, "%s%s", outputPath, paramInfo);

    // write data to binary files
    FILE * fp;
    fp = fopen(TFilename, "wb");
    for (unsigned q=0; q<data->numberOfT; ++q) {
        for (unsigned n=0; n<data->N/data->timeskip; ++n)
        fwrite(res->T[q][n], sizeof(double), data->J[q], fp);
    }
    fclose(fp);

    fp = fopen(IFilename, "wb");
    for (unsigned q=0; q<data->numberOfIv; ++q) {
        for (unsigned n=0; n<data->N*data->ETratio; n += data->timeskip/10)
            fwrite(&res->Iv[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    fp = fopen(RFilename, "wb");
    for (unsigned q=0; q<data->numberOfT; ++q) {
        for (unsigned n=0; n<data->N*data->ETratio; n += data->timeskip/10)
            fwrite(&res->R_w[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    fp = fopen(CFilename, "wb");
    for (unsigned q=0; q<data->numberOfC; ++q) {
        for (unsigned n=0; n<data->N*data->ETratio; n += data->timeskip/10)
            fwrite(&res->V_c[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    // write simulation parameters to file for readout
    fp = fopen(paramInfoFilename, "w");
    fprintf(fp, "%40s\n", "SIM MEMORY PARAMETERS");
    fprintf(fp, "%40s; %zu\n", "number of nanowires", data->numberOfT);
    fprintf(fp, "%40s; %zu\n", "number of currents", data->numberOfIv);
    fprintf(fp, "%40s; %zu\n", "number of capacitors", data->numberOfC);
    fprintf(fp, "%40s; %zu\n", "number of resistor groups", data->groupsOfR);
    fprintf(fp, "%40s; ", "number of resistors per group, in order");
    for (unsigned i=0; i<data->groupsOfR; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%zu", data->numberOfR[i]);
    }
    fprintf(fp, "\n%40s; %zu\n", "number of inductor groups", data->groupsOfL);
    fprintf(fp, "%40s; ", "number of inductors per group, in order");
    for (unsigned i=0; i<data->groupsOfL; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%zu", data->numberOfL[i]);
    }
    fprintf(fp, "\n%40s; %zu\n", "number of bias currents", data->numberOfIv_b);

    fprintf(fp, "%40s\n", "SIM TIMING PARAMETERS");
    fprintf(fp, "%40s; %zu\n", "N (# of temporal elements)", data->N);
    fprintf(fp, "%40s; %10.3e\n", "tMax (sim time for thermal model)", data->tMax);
    fprintf(fp, "%40s; %10.3e\n", "delta t", res->dt);
    fprintf(fp, "%40s; %zu\n", "timeskip factor", data->timeskip);
    fprintf(fp, "%40s; %zu\n", "electrical / thermal time ratio", data->ETratio);
    fprintf(fp, "%40s; %d\n", "allow thermal optimization", data->allowOpt);

    fprintf(fp, "%40s\n", "THERMAL PARAMETERS");
    fprintf(fp, "%40s; %10.3e\n", "phonon specific heat", data->c_p);
    fprintf(fp, "%40s; %10.3e\n", "electron specific heat", data->c_e);
    fprintf(fp, "%40s; %10.3e\n", "thermal boundary conductivity", data->alpha);
    fprintf(fp, "%40s; %10.3e\n", "reference temperature for thermal", data->T_ref);
    fprintf(fp, "%40s; %10.3e\n", "sheet resistance", data->R_gamma);
    fprintf(fp, "%40s; %10.3e\n", "critical temperature", data->T_c);
    fprintf(fp, "%40s; %10.3e\n", "substrate temperature", data->T_sub);
    fprintf(fp, "%40s; %10.3e\n", "sub temp epsilon, optimization strategy", data->T_sub_eps);

    fprintf(fp, "%40s\n", "BIAS CURRENTS");
    fprintf(fp, "%40s; ", "bias currents, in order");
    for (unsigned i=0; i<data->numberOfIv_b; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->Iv_b[i]);
    }

    fprintf(fp, "\n%40s\n", "NANOWIRE PARAMETERS");
    fprintf(fp, "%40s; %10.3e\n", "impurity offset", data->impurityOffset);
    fprintf(fp, "%40s; %10.3e\n", "impurity spread", data->impuritySpread);
    fprintf(fp, "%40s; ", "J (# of spatial elements per wire)");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9zu", data->J[i]);
    }
    fprintf(fp, "\n%40s; ", "wire lengths");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->wireLength[i]);
    }
    fprintf(fp, "\n%40s; ", "wire thicknesses");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->wireThickness[i]);
    }
    fprintf(fp, "\n%40s; ", "wire widths");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->wireWidth[i]);
    }
    fprintf(fp, "\n%40s; ", "delta X's");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", res->dX[i]);
    }
    fprintf(fp, "\n%40s; ", "wire inductances");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->L_w[i]);
    }
    fprintf(fp, "\n%40s; ", "wire critical currents");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->Iv_c0[i]);
    }
    fprintf(fp, "\n%40s; ", "corresponding current per wire");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9d", data->cor_Iv[i]);
    }

    fprintf(fp, "\n%40s\n", "INITIAL CONDITIONS");
    fprintf(fp, "%40s; ", "initial currents");
    for (unsigned i=0; i<data->numberOfIv; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->Iv_init[i]);
    }
    fprintf(fp, "\n%40s; ", "initial voltages");
    for (unsigned i=0; i<data->numberOfC; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->V_c_init[i]);
    }
    fprintf(fp, "\n%40s; ", "triggered nanowires");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9d", data->trigger[i]);
    }
    fprintf(fp, "\n%40s; ", "initial hot-spot size");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->HS_l_init[i]);
    }
    fprintf(fp, "\n%40s; ", "initial hot-spot temperature");
    for (unsigned i=0; i<data->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%9.2e", data->HS_T_init[i]);
    }

    fprintf(fp, "\n%40s\n", "TRANSMISSION LINE PARAMETERS");
    fprintf(fp, "%40s; %d\n", "transmission type", data->simTL);
    fprintf(fp, "%40s; %zu\n", "NTL (# of transmission line elements)", data->NTL);
    fprintf(fp, "%40s; %10.3e\n", "velocity factor", data->VF);
    fprintf(fp, "%40s; %10.3e\n", "transmission line length", data->LTL);

    fclose(fp);

    free(TFilename);
    free(IFilename);
    free(RFilename);
    free(CFilename);
    free(paramInfoFilename);
}

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

    puts("\n    WRITING RESULTS...");
    write_results(argv[2], data, res);
    puts("\n    DONE.\n");

    free_simres(data, res);
    free_simdata(data);

    puts("########################################################################################\n#################################### END SIMULATION ####################################\n########################################################################################");

    exit(0);
}
