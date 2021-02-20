#ifndef INTEGRATOR_H
#define INTEGRATOR_H


#include "../kissfft/kiss_fft.h"

/*******************************************************************************
Function that calls the right integrator depending on the integ_choice
0 = Euler ; 1 = UCM ; 2 = SSM
*******************************************************************************/
void integrator(double complex* in, double tau, int integ_choice);

/*******************************************************************************
Calculates one wavefunction time step by euler method
*******************************************************************************/
/* error of O(tau^2) */
void euler_method(double complex *in, double tau);


/* Function used by UCN_method */
void operator(double complex *in, double complex *out);

/*******************************************************************************
Calculates one wavefunction time step by Unitary Crank Nicolson method
-> uses conjugate_gradient() method
*******************************************************************************/
void UCN_method(double complex *in,double tau);

/*******************************************************************************
Calculates one wavefunction time step by Strang Splittin Method
-> see integrator.c for more
-> uses kissfft library for FFT/IFFT
*******************************************************************************/
void strangsplitting_method(double complex *in, double tau);

#endif //INTEGRATOR_H
