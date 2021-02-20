#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "../kissfft/kiss_fft.h"


typedef struct {
  kiss_fft_cpx *cx_in, *cx_out;
  kiss_fft_cfg cfg, icfg;
}kissfft_struct;

/* error of O(tau^2) */
void euler_method(double complex *in,double tau);

/* Function used by UCN_method */
void operator(double complex *in, double complex *out);

/* error of O(tau^3) */
void UCN_method(double complex *in,double tau);

/*******************************************************************************
Strang splitting and necessary functions
*******************************************************************************/
void strangsplitting_method(double complex *in, double tau);
// calucalates one time iteration step
// uses kiss_fft

#endif //INTEGRATOR_H
