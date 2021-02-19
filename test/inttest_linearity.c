#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "wavefunction.h"
#include "integrator.h"
#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h" //shouldnt be necessary if imported by the integrator


static int N = 0;

/*******************************************************************************
// Linearity test

// goal:
//  1x time evolution for ||int(a*psi + b*theta) - a*int(psi) + b*int(theta) ||
calculate this norm for every time step (here 1 time step is enough)

create random wavefunctions psi and theta
create random coefficients

for (int q=0; q<Q; q++) { //q is the current multiplication of tau, Q*tau = t
    \\ integrator(a*psi + b*theta, q) - a*integrator(psi,q) + b*integrator(theta,q) ||

    print( calculated value to some file and/or to terminal)
}
*******************************************************************************/
void integrator(double complex* in, double tau, int integ_choice) {
    /***************************************************************
   function that calculates the time evolution of input wavefunciton
   for a chosen integrator method
   ****************************************************************/
  if (integ_choice == 0) {
    euler_method(in, tau);
    printf("Integrator used: Euler Method!\n");
  }
  else if (integ_choice == 1) {
    UCN_method(in, tau);
    printf("Integrator used: Unitary Crank Nicolson Method!\n");
  }
  else if (integ_choice == 2) {
    // init_strangsplitting();
    strangsplitting_method(in, tau);
    printf("Integrator used: Strang Splitting Method!\n");
    // finished_strangsplitting();
  }
  else {
    printf("[inttest_linearity.c | integrator()] Error! Choice of integrator is out of range!\n"
  "Remember: Integrator choice is 3rd input when calling inttest_linearity.c.\n "
  "Euler Method = 1, Unitary Crank-Nicolson Method = 2, Strang Splitting Method = 2\n" );

  exit(-1);
  }
}


int main(int argc, char const *argv[]) {
/****************************************************************
argv[1] = N
argv[2] = tau
argv[3] = integrator_choice
argv[4] = potential
****************************************************************/

  // assert(argc==5,_FILE_NAME_,"main","ERROR. Necessary number of input parameters 4!\n
  // Usage: inttest_linearity {N} {tau} {integrator_choice} {potential_choice}\n
  // Remember: The executable Name is the first parameter.\n");

  printf("This program calculates: norm(delta) = norm( ||int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|| )\n\n" );
  // if (argc == 4) {
  //   printf(
  //     "\n[ inttest_linearity.c | input parameters ] ERROR. Necessary number of input parameters 4!\n"
  //     "You have %d!\n" "Remember: The executable Name is the first parameter.\n", argc
  //   );
  //   exit(-1);
  // }

  srand(time(NULL)); // is called in beginning of the main.c
  set_params(argc, (char**) argv);
  N = get_N();
  double tau = get_tau();
  int integrator_choice = get_integ_choice();
  double mass = 2.3512;
  set_kinetic_params(mass);
  int pot = get_pot_choice();
  set_potential(pot);
  print_hamiltonian_info();

  double complex *psia, *psib;
  double complex alpha,  beta;
  double complex *left, *right;
  double maxdev, dev;

  /* dynamic memory allication (remember to free at the end) */
  psia = (double complex*) malloc(sizeof(double) * N);
  psib = (double complex*) malloc(sizeof(double) * N);
  left = (double complex*) malloc(sizeof(double) * N);
  right = (double complex*) malloc(sizeof(double) * N);
  assert(psia!=NULL); // check of condition is true else -> aboort?
  assert(psib!=NULL);
  assert(left!=NULL);
  assert(right!=NULL);

  /* generate random normalized wavefunctions psia and psib */
  set_random_wavefunction(psia, N);
  set_random_wavefunction(psib, N);

  printf(
    "Linearity test: wavefuncitons:\n"
    "\tPsia[0] = %.e + %.e * i \n" "\tPsib[0] = %.e + %.e * i\n",
    creal(psia[0]), cimag(psia[0]), creal(psib[0]), cimag(psib[0])
  );

  /* generate random coomplex coefficients */
  alpha = random_complex_coefficient();
  beta = random_complex_coefficient();

  printf(
    "Linearity test: coefficients:\n"
    "\talpha= %.e + %.e * i \n" "\tPsib= %.e + %.e * i\n",
    creal(alpha), cimag(alpha), creal(beta), cimag(beta)
  );

  /* Calculate int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib) */
  /* Left hand side */
  if (integrator_choice==2) {init_strangsplitting();}
  for(int n = 0; n < N; n++) {
    left[n] = alpha*psia[n] + beta*psib[n];
  }

  integrator(left, tau, integrator_choice);

  /* Right hand side */
  integrator(psia, tau, integrator_choice);
  integrator(psib, tau, integrator_choice);
  if (integrator_choice==2) {finished_strangsplitting();}
  /* ...and difference */
  maxdev = 0.0;
  for(int n = 0; n < N; n++) {
    right[n] = alpha*psia[n] + beta*psib[n];
    dev = cabs(left[n] - right[n]);
    if(dev>maxdev)  maxdev=dev;
  }

  printf("Calcualted difference for every n: |int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|");
  printf(
    "Maximum Deviation is=\t%.e\n"
    "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n", maxdev);

  /* Printing test infos to text file */
  FILE *fp;
  int namesize = 45;
  for (int i=1; i<=4; i++) { namesize += strlen(argv[i]); }
  char filename[namesize];

  snprintf(
    filename, sizeof(filename),
    "data/int_lin_test_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4]
  );

  fp = fopen(filename, "w");
  fprintf(fp, "n\tREAL(psia[n])\tIMAG(psia[n])\tREAL(psib[n])\tIMAG(psib[n])\n");
  for (int i=0; i<N; i++) {
    fprintf(fp, "%d\t%.e\t%.e\t%.e\t%.e\n", i, creal(psia[i]), cimag(psia[i]), creal(psib[i]), cimag(psib[i]));
  }
  fprintf(
    fp, "Coefficients:\n" "alpha = %.e + %.e * i\n" "beta = %.e + %.e * i\n",
    creal(alpha), cimag(alpha), creal(beta), cimag(beta)
  );
  if (pot == 0) {fprintf(fp, "Potential used:\tZERO potential\n");}
  else if (pot == 1) {fprintf(fp, "Potential used:\tHARMONIC potential\n");}
  else if (pot == 2) {fprintf(fp, "Potential used:\tWELL potential\n");}
  else if (pot == 3) {fprintf(fp, "Potential used:\tWALL potential\n");}
  else {fprintf(fp, "Potential used:\tERROR\n");}

  // fprintf(
  //   fp, "norm(delta) = norm( ||int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|| )\n"
  //   "norm(delta) = %.e\n", norm(delta, N)
  // );
  fprintf(fp, "MAXIMUM DEVIATION =\t%.e\n" "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n", maxdev);

  fclose(fp);

  /* Free allicated wavefunctions */
  free(psia);
  free(psib);
  free(left);
  free(right);

  return 0;
}
