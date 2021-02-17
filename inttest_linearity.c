#include <stdio.h>

#include "wavefunction.h"
#include "integrator.h"

/*******************************************************************************
// Linearity test

// goal:
//   time evolution for ||H(a*psi + b*theta) - a*H(psi) + b*H(theta) ||
calculate this norm for every time step

create random wavefunctions psi and theta
create random coefficients

for (int q=0; q<Q; q++) { //q is the current multiplication of tau, Q*tau = t
    \\ integrator(a*psi + b*theta, q) - a*integrator(psi,q) + b*integrator(theta,q) ||

    print( calculated value to some file and/or to terminal)
}
*******************************************************************************/
void integrator(double complex* in, double tau, int integ_choice) {
  if (integ_choice == 0) {
    euler_method(in, tau);
  }
  else if (integ_choice == 1) {
    UCN_method(in, tau);
  }
}


int main(int argc, char const *argv[]) {
/*
argv[0] = N
argv[1] = tau
argv[2] = integrator_choice
*/


  int N = atoi(argv[0]);
  double tau = atoi(argv[1]);
  int integrator_choice = atoi(argv[2]);
  double complex *psia, *psib;
  double complex alpha,  beta;
  double complex *in_left, *in_right, *out_left, *out_right, *delta;

  /* dynamic memory allication (remember to free at the end) */
  psia = (double*) malloc(sizeof(psia) * N);
  psib = (double*) malloc(sizeof(psib) * N);
  in_left = (double*) malloc(sizeof(in_left) * N);
  in_right = (double*) malloc(sizeof(in_right) * N);
  out_left= (double*) malloc(sizeof(out_left) * N);
  out_right = (double*) malloc(sizeof(out_right) * N);
  delta = (double*) malloc(sizeof(delta) * N);

  /* generate random normalized wavefunctions psia and psib */
  set_random_wavefunction(psia, N);
  set_random_wavefunction(psib, N);

  /* generate random coomplex coefficients */
  random_complex_coefficient(alpha);
  random_complex_coefficient(beta);

  /* Calculate int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib) */
  for(int n = 0; n < N; n++) {
    in_left[n] = alpha*psi[n] + beta*theta[n];
  }

  if ( integrator_choice == 0) {

  }
H(in_left,out_left);//computes the left hand side

H(psi,out_psi);
H(theta,out_theta);

  for(int i = 0; i < N; i++) {
    out_right[i] = alpha*out_psi[i] + beta*out_theta[i];
    delta[i] = out_left[i] - out_right[i];
  }

  return 0;
}
