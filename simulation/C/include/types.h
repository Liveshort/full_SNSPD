#ifndef __TYPES_H__
#define __TYPES_H__

// some global constants
const double Kb;            // boltzmann constant
const double Lorentz;       // Lorentz number

enum transmissionType{NOTRNS, DELAYTRNS, VARRPLCTRNS, CONSTRPLCTRNS, NORPLCTRNS};

// structure that contains the result of the simulation of an snspd
typedef struct _simres {
    double * dX;            // delta X
    double dt;              // delta t
    double * rho_norm;      // resistivity of the nanowire in the normal state
    double *** T;           // temperature matrices T[number][time][space]
    double ** Iv;           // current vectors Iv[number][time]
    double ** R_w;          // resistance vectors R_w[number][time]
    double ** V_c;          // voltage vectors for capacitors V_c[number][time]
} SimRes;

// structure that contains the (input) data of the simulation of an snspd
typedef struct _simdata {
    // sim parameters
    size_t numberOfT;       // number of nanowires
    size_t numberOfIv;      // number of currents
    size_t numberOfC;       // number of capacitors
    size_t groupsOfR;       // number of resistor groups
    size_t * numberOfR;     // number of resistors per group, in order
    size_t groupsOfL;       // number of inductor groups
    size_t * numberOfL;     // number of inductors per group, in order
    size_t numberOfIv_b;    // number of bias currents
    // timing parameters
    size_t N;               // number of time samples
    double tMax;            // maximum time, sim will run from 0 to tMax
    size_t timeskip;        // saved time step reduction
    size_t ETratio;         // ratio to which electrical is calculated more than thermal
    int allowOpt;           // optimizes thermal when wire has cooled to within T_sub_eps
    // thermal parameters
    double c_p;             // phonon specific heat
    double c_e;             // electron specific heat
    double alpha;           // thermal boundary conductivity
    double T_ref;           // reference temperature for thermal values
    double R_gamma;         // sheet resistance of the nanowire in normal state
    double T_c;             // nominal critical temperature
    double T_sub;           // substrate temperature
    double T_sub_eps;       // sub temp epsilon, optimization strategy (detect steady state)
    // bias currents
    double * Iv_b;          // bias currents
    // nanowire parameters
    double impurityOffset;  // minimum impurity of Iv_c
    double impuritySpread;  // spread in impurity of Iv_c
    size_t * J;             // spatial elements per nanowire, in order
    double * wireLength;    // length of the nanowires, in order
    double * wireThickness; // thickness of the nanowires, in order
    double * wireWidth;     // width of the nanowires, in order
    double * L_w;           // inductances of the nanowires, in order
    double * Iv_c0;         // critical currents of the nanowires, in order
    unsigned * cor_Iv;      // corresponding current per nanowire, in order
    // other component parameters
    double ** R;            // resistor values, per group, in order
    double ** L_x;          // inductor values, per group, in order
    double * C;             // capacitor values, in order
    // initial conditions
    double * Iv_init;       // initial currents, in order
    double * C_init;        // initial voltages over the capacitors, in order
    int * trigger;          // which wires to trigger
    double * initHS_l;      // initial hotspot size per nanowire, in order
    double * initHS_T;      // initial hotspot temperature per nanowire, in order
    // transmission line parameters
    int simTL;              // 0: none, 1-4: delay, implicit, explicit par R, explicit
    size_t NTL;             // number of transmission line elements
    double VF;              // velocity factor transmission line. Calculated as 1/sqrt(eps_r)
    double LTL;             // length of transmission line in [m]
} SimData;

void free_simres(SimData * data, SimRes * res);
void free_simdata(SimData * data);

#endif
