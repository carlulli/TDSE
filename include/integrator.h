#ifndef INTEGRATOR_H
#define INTEGRATOR_H


#include "../kissfft/kiss_fft.h"

void integrator(double complex* in, double tau, int integ_choice);

/* takes extra parameter intg_choice to choose the type of integrator it uses */
void integrator(double complex* in, double tau, int integ_choice);

/* error of O(tau^2) */
void euler_method(double complex *in, double tau);

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
