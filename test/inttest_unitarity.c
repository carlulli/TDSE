
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"


#define _FILE_NAME_ "test/inttest_unitarity.c"

static int N = 0;

/******************************************************************
//unitarity_test.c

// goal:  time evolution for the time difference
// calculate the norm and write it on a file as a function of time paramter discretized via the time step tau


0) take as input     N    TIME    ITER    INTEGRATOR    POTENTIAL
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
  printf("\nThe program generates a random NON normalized wavefunction. Then it applies the chosen integrator to approximate the time evolution of the wavefunction.\n\n");
  printf("The Test params are: N=%d, mass=%f, tau=%.e, integrator_choice=%d\n\n",N, mass, tau, integrator_choice);
  double complex psi[N],psi_in[N];
  /* the wavefunction generated is not normalized by default */
  set_random_wavefunction_NN(psi,N);
  for(int i = 0; i < N; i++) {
    psi_in[i] = psi[i];
  }
  //printf("Computing the norm of the vector psi \n ||psi|| = %f \n", norm(psi,N));
  /* Now we apply the integrator chosen by the user,
  we want to check that the time operator is unitary and that psi keeps the same norm!  (not Euler)*/
  for(int i = 0; i < N; i++) {
  integrator(psi,tau,integrator_choice);
  }
  printf("| ||psi(tf)|| - ||psi(ti)|| | = %.12e\ttau = %f\tint_method = %d\n",abs(norm(psi,N)-norm(psi_in,N)),tau,integrator_choice);
  }


  /* Printing test infos to text file */
  /*
  FILE *fp;
  int namesize = 60;
  for (int i=1; i<=8; i++) { namesize += strlen(argv[i]); }
  char filename[namesize];

  snprintf(
    filename, sizeof(filename),
    "data/int_uni_test_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);

  fp = fopen(filename, "w");
  fprintf(fp, "n\tREAL(psi[n])\tIMAG(psi[n])\ttau\taveren\taverx\tdeltax\taverp\tdeltap\n");
  for (int q = 0; q < nsteps; j++) {
    integrator(psi,tau,integrator_choice);
    for (int i = 0; i < N; i++) {
    fprintf(fp, "%.e\t%.e\t%.e\t%.e\t%.e\t%.e\t%.e\t%.e\n", creal(psi[i]),cimag(psi[i]),tau*q,average_state_energy(psi),get_avgx(psi),get_deltax(psi),get_avgp(psi),get_deltap(psi));
  }
  fclose(fp);
  */
  return 0;
}
