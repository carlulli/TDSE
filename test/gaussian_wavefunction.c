//gaussian_wavefunction.case
//goal : test that set_gaussian_wavefunction correctly generates it


#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
//#include "hamiltonian.h"
#include "wavefunction.h"
//#include "integrator.h"


#define _FILE_NAME_ "test/inttest_convergence.c"

/*********************************************
1) Generates a gaussian wavepacket with a set of hardcored parameters
2) Print them on terminal and checks norm;
3) Can eventually be printed on file and plotted with python
***********************************************/

int main(int argc, char *argv[]) {

  int N = 65;
  double complex psi[N];
  double time = 100;
  double nsteps = 100;
  double tau = time/nsteps;
  double mu = 20;
  double sigma = 5;
  double pbar = 2;
  set_gaussian_wavefunction(psi,mu,sigma,pbar,N);

  for(int i = 0; i < N; i++) {
  printf("psi[i] = %f  (I) %f\n", creal(psi[i]), cimag(psi[i]));
  }
  printf("norm = %f",norm(psi,N));

  return 0;
}
