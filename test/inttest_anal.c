#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"
#include "assert.h"

#define _FILE_NAME_ "test/inttest_anal.c"

/* This test compares the TDSE obtained integrating with different methods against a SE with V = 0
  In this test a randomwavefunction is generated and the Hamiltonian is set with V = 0
  The user inputs the total simulation time t, and the time step tau.
  The energy of the wavefunction is printed on a text file for all 3 integration methods the time step is increased
  up to the value of time. The number of steps to get there is passed as a parameter by the user
  Different runs are made at different time steps tau and the continuum limit is then extrapolated via linear fit */


int main (int argc, char *argv[]) {
  /*checks that the input parameters are correct */
  assert(argc== 5,_FILE_NAME_,"inttest_anal","Usage: Please input [N]  [time]  [tau] [Iteration]\n");
  /* sets N as the number of lattice points */
  set_params(argc, argv);
  /* sets the potential in H as V = 0 for this test */
  set_zero_potential();
  double time = atof(argv[2]);
  double tau = atof(argv[3]);
  int iter = atoi(argv[4]);
  int N = get_N();
  /* generates initial random normalized wavefunction */
  double complex psi[N];


  set_random_wavefunction(psi,N);

  FILE *fptr;


  return 0;
}
