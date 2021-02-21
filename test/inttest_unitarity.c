
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"
#include "assert.h"

#define _FILE_NAME_ "test/inttest_unitarity.c"

static int N = 0;

/******************************************************************
//unitarity_test.c

// goal:  time evolution for the time difference
// calculate the norm and write it on a file as a function of time paramter discretized via the time step tau


0) take as input     N    TIME    NSTEP    INTEGRATOR    POTENTIAL
1) generate a random normalized initial wavefunction
2) compute the norm of psi(t) and check that it is still 1
3)print on terminal and or on a file
*******************************************************************/

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
  set_kinetic_params(mass);
  set_potential(pot_choice);
  print_hamiltonian_info();



  printf("This program is tasked with testing that the various integration methods implement respect the unitarity of the time evolution operator (only those who uphold this property should be ran) \n");
  printf("\nThe program generates a random NON normalized wavefunction. Then it applies the chosen integrator to approximate the time evolution of the wavefunction.\nAt each step (in regard to the parameters passed by the user) the norm is computed and printed on the screen. Here the user can check that the norm is maintained through time (up to a certain tolerance)\n");
  printf("The test passes is all these numbers are < 1e-15\n\nTest params are: N=%d, mass=%f, tau=%.e, integrator_choice=%d\n\n",N, mass, tau,integrator_choice);
  double complex psi[N];
  /* the wavefunction generated is not normalized by default */
  set_random_wavefunction_NN(psi,N);
  printf("Computing the norm of the vector psi \n ||psi|| = %f \n", norm(psi,N));
  /* Now we apply the integrator chosen by the user,
  we want to check that the time operator is unitary and that psi keeps the same norm!  (not Euler)*/
  for(int i = 0; i < nsteps; i++)
    if (integrator_choice == 1) {
      UCN_method(psi,tau);
      printf("Norm of vector at time %f  with UCN_method  = %f\n", i*tau, norm(psi,N));
  }
  /*  else if( int_type == 2) {
      splitting_method(psi,tau);
      printf("Norm of vector at time %f   with splitting_method  = %f\n",i*tau, norm(psi,N));

  }
*/
/* Printing test infos to text file */

  return 0;
}
