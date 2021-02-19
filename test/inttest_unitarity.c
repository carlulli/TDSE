/******************************************************************
//unitarity_test.c

// goal:  time evolution for the time difference
// calculate the norm and write it on a file as a function of time paramter discretized via the time step tau


0) take as input     N    TIME    ITER    INTEGRATOR    POTENTIAL
1) generate a random normalized initial wavefunction
2) compute the norm of psi(t) and check that it is still 1
3)print on terminal and or on a file

*******************************************************************/

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



int main( int argc, char *argv[]) {

  /* checks that the input parameters are correct */
  assert(argc == 7, _FILE_NAME_, "inttest_anal","Usage: Please input: [N]  [mass]  [time]  [N-steps]  [integrator]  [potential]\n");
  /* sets N as the number of lattice points */
  set_params(argc, argv);
  set_kinetic_params(atof(argv[2]));
  /* set potential */
  set_potential(atoi(argv[5]));
  /* set the integrator type */
  int int_type = atoi(argv[4]);
  /* define tau the time step as time/Nsteps */
  double time = atof(argv[3]);
  int nsteps = atoi(argv[4]);
  double tau = time/nsteps;
  /* sets the N size of the wavefunction psi[N] */
  int N = get_N();
  double complex psi[N];
  /* the wavefunction generated is not normalized by default */
  set_random_wavefunction_NN(psi,N);
  printf("Computing the norm of the vector psi \n ||psi|| = %f \n",norm(psi,N));
  /* Now we apply the integrator chosen by the user,
  we want to check that the time operator is unitary and that psi keeps the same norm!  (not Euler)*/
  for(int i = 0; i < nsteps; i++)
    if (int_type == 1) {
      UCN_method(psi,tau);
      printf("Norm of vector at time %f   with UCN_method  = %f\n",i*tau, norm(psi,N));
  }
  /*  else if( int_type == 2) {
      splitting_method(psi,tau);
      printf("Norm of vector at time %f   with splitting_method  = %f\n",i*tau, norm(psi,N));

  }
*/
  else {
    printf("WRONG INTEGRATOR INPUT PLEASE INSERT ONLY [1] for UCN OR [2] for Strang Splitting\n");
    return(-1);
  }





  return 0;
}
