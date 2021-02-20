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

  /* Here the parameters are passed in the set_params function, which checks the validity */
  set_params(argc, (char**) argv);
  /* mass is always hardcoded */
  double mass = 2.3512;
  N = get_N();
  double time = get_time();
  double nsteps = get_nsteps();
  double tau = time/nsteps;
  int integrator_choice = get_integ_choice();
  int pot_choice = get_pot_choice();
  double mu = atof(argv[6]);
  double dx = atof(argv[7]);
  double dp = atof(argv[8]);
  set_kinetic_params(mass);
  set_potential(pot_choice);
  print_hamiltonian_info();
  set_gaussian_wavefunction(psi,mu,dx,dp,N);

  for(int i = 0; i < N; i++) {
  printf("psi[i] = %f  (I) %f\n", creal(psi[i]), cimag(psi[i]));
  }
  printf("norm = %f",norm(psi,N));

  return 0;
}
