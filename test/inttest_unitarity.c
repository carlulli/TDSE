
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <time.h>

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
  srand(time(NULL));
  double mass = 2.3512;
  N = get_N();
  double ttime = get_time();
  double nsteps = get_nsteps();
  double tau = ttime/nsteps;
  int integrator_choice = get_integ_choice();
  int pot_choice = get_pot_choice();
  set_kinetic_params(mass);
  set_potential(pot_choice);
  print_hamiltonian_info();



  printf("This program tests if the integration methods implemented respect the unitarity of the time evolution operator\n"
  "(only those who uphold this property should be ran) \n");
  printf("\nThe program generates a random NON normalized wavefunction.\n"
  "Then it applies the chosen integrator to approximate the time evolution of the wavefunction.\n"
  "The norm is computed and printed on the screen and to a file.\n"
  "Here the user can check that the norm is maintained through time (up to a certain tolerance)\n");
  printf("The test passes if all these numbers are < 1e-15\n\n"
  "Test params are: N=%d, mass=%f, tau=%.e, integrator_choice=%d\n\n",N, mass, tau,integrator_choice);

  double complex psi[N], psi_cp[N];
  /* the wavefunction generated is not normalized by default */
  set_random_wavefunction(psi,N);
  copy_wf(psi, psi_cp, N);
  printf("Computing the norm of the vector psi \n ||psi|| = %f \n\n", norm(psi,N));
  /* Now we apply the integrator chosen by the user,
  we want to check that the time operator is unitary and that psi keeps the same norm!  (not Euler)*/

  /* isnt one time step enough? */
  integrator(psi, tau, integrator_choice);
  printf("Norm of vector at time 1 step  with integrator:%d \n ||psi|| = %f\n", integrator_choice, norm(psi,N));
  printf("\n||psiold|| - ||psinew|| = %.6e\n", fabs(norm(psi_cp, N)-norm(psi,N)));

  /* Printing test infos to text file */
  FILE *fp;
  int namesize = 45;
  for (int i=1; i<=5; i++) { namesize += strlen(argv[i]); }
  char filename[namesize];

  snprintf(
    filename, sizeof(filename),
    "data/int_uni_test_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4] ,argv[5]);

  fp = fopen(filename, "w");
  fprintf(fp, "n\tREAL(psi_cp[n])\tIMAG(psi_cp[n])\n");
  for (int i=0; i<N; i++) { fprintf(fp, "%d\t%.6e\t%.6e\n", i, creal(psi_cp[i]), cimag(psi_cp[i])); }
  fprintf(fp, "\nn\tREAL(int(psi[n]))\tIMAG(int(psi[n]))\n");
  for (int i=0; i<N; i++) {  fprintf(fp, "%d\t%.6e\t%.6e\n", i, creal(psi[i]), cimag(psi[i])); }

  if (pot_choice == 0) {fprintf(fp, "\nPotential used:\tZERO potential\n");}
  else if (pot_choice == 1) {fprintf(fp, "\nPotential used:\tHARMONIC potential\n");}
  else if (pot_choice == 2) {fprintf(fp, "\nPotential used:\tWELL potential\n");}
  else if (pot_choice == 3) {fprintf(fp, "\nPotential used:\tWALL potential\n");}
  else {fprintf(fp, "\nPotential used:\tERROR\n");}

  if (integrator_choice == 0) {fprintf(fp, "\nIntegrator used:\tEuler Method\n");}
  else if (integrator_choice == 1) {fprintf(fp, "\nIntegrator used:\tUCM Method\n");}
  else if (integrator_choice == 2) {fprintf(fp, "\nIntegrator used:\tSSM Method\n");}
  else {fprintf(fp, "\nIntegrator used:\tERROR\n");}

  fprintf(fp, "\nNorm before\tNorm after int\n" "%.6e\t%.6e\n", norm(psi_cp,N), norm(psi,N));
  fprintf(fp, "\n||psiold|| - ||psinew|| = %.6e\n", fabs(norm(psi_cp, N)-norm(psi,N)));


  fclose(fp);

  printf("\nUNITARITY test finished.\n\n");
  return 0;
}
