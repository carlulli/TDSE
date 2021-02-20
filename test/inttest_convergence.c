//inttest_convergence.c
// goal:  time evolution for the time difference
// calculate the norm and write it on a file as a function of time paramter discretized via the time step tau

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"


#define _FILE_NAME_ "test/inttest_convergence.c"
static int N = 0;
/***************************************************************

1) generate a guassian wavefunction
2) compute time evolution with the chosen potential and integrator
3)
*******************************************************************/


int main( int argc, char *argv[]) {

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
  set_kinetic_params(mass);
  set_potential(pot_choice);
  print_hamiltonian_info();



  printf("Test prams are: N=%d, mass=%f, time= %f nsteps=%d, integrator_choice=%d, potential = %d mu = %d sigma= %f\n\n",
  N, mass, time, nsteps, integrator_choice,pot_choice,mu,sigma);
  printf("In this test we generate a gaussian wavepacket with the parameters mu and sigma passed by the user\n");
  printf("Then we compute the ")
  double complex psi[N];
  /* the wavefunction generated is not normalized by default */
  set_random_wavefunction_NN(psi,N);
  return 0;
}
