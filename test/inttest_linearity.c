#include <stdio.h>
#include <time.h>

#include "wavefunction.h"
#include "integrator.h"
#include "geometry.h"
#include "linearalgebra.h"
#include "assert.h"

#define _FILE_NAME_ "test/inttest_linearity.c"

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
    strangsplitting_method(in, tau);
    printf("Integrator used: Strang Splitting Method!\n");
  }
  else {
    printf("[inttest_linearity.c | integrator()] Error! Choice of integrator is out of range!\n"
  "Remember: Integrator choice is 3rd input when calling inttest_linearity.c.\n ",
  "Euler Method = 1, Unitary Crank-Nicolson Method = 2, Strang Splitting Method = 2\n" );

  exit(-1);
  }
}


int main(int argc, char const *argv[]) {
/****************************************************************
argv[0] = N
argv[1] = tau
argv[2] = integrator_choice
****************************************************************/

  assert(argc==4,_FILE_NAME_,"main","ERROR. Necessary number of input parameters 4!\n
  Usage: inttest_linearity {N} {tau} {integrator_choice}\n Remember: The executable Name is the first parameter.\n");

  // if (argc == 4) {
  //   printf(
  //     "\n[ inttest_linearity.c | input parameters ] ERROR. Necessary number of input parameters 4!\n"
  //     "You have %d!\n" "Remember: The executable Name is the first parameter.\n", argc
  //   );
  //   exit(-1);
  // }

  srand(time(NULL)); // is called in beginning of the main.c
  set_params(argc, argv);
  int N = get_N();
  double tau = atoi(argv[2]);
  int integrator_choice = atoi(argv[3]);
  double complex *psia, *psib;
  double complex alpha,  beta;
  double complex *left, *right, *delta;

  /* dynamic memory allication (remember to free at the end) */
  psia = (double*) malloc(sizeof(psia) * N);
  psib = (double*) malloc(sizeof(psib) * N);
  left = (double*) malloc(sizeof(left) * N);
  right = (double*) malloc(sizeof(right) * N);
  delta = (double*) malloc(sizeof(delta) * N);

  /* generate random normalized wavefunctions psia and psib */
  set_random_wavefunction(psia, N);
  set_random_wavefunction(psib, N);

  printf(
    "Linearity test: wavefuncitons:\n"
    "\tPsia[0] = %.e + %.e * i \n" "\tPsib[0] = %.e + %.e * i\n",
    creal(psia[0]), cimag(psia[0]), creal(psib[0]), cimag(psib[0])
  );

  /* generate random coomplex coefficients */
  random_complex_coefficient(alpha);
  random_complex_coefficient(beta);

  printf(
    "Linearity test: coefficients:\n"
    "\talpha= %.e + %.e * i \n" "\tPsib= %.e + %.e * i\n",
    creal(psia), cimag(psia), creal(psib), cimag(psib)
  );

  /* Calculate int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib) */
  /* Left hand side */
  for(int n = 0; n < N; n++) {
    left[n] = alpha*psi[n] + beta*theta[n];
  }

  integrator(left, tau, integrator_choice);

  /* Right hand side */
  integrator(psia, tau, integrator_choice);
  integrator(psib, tau, integrator_choice);

  /* ...and difference */
  for(int n = 0; n < N; n++) {
    right[n] = alpha*psia[n] + beta*psib[i];
    delta[i] = left[i] - right[i];
  }

  printf("Calcualted difference ||int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|| = %.e \n", norm(delta, N));
  printf("Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n");

  /* Printing test infos to text file */
  FILE *fp;
  int namesize = 19;
  for (int i=1; i<=3; i++) { namesize += strlen(argv[i]); }
  char filename[namesize];

  snprintf(
    filename, sizeof(filename),
    "data/int_lin_test_%s_%s_%s.txt", argv[1], argv[2], argv[3]
  );

  fp = fopen(filename, "w");
  fprintf(fp, "n\tREAL(psia[n])\tIMAG(psia[n])\tREAL(psib[n])\tIMAG(psib[n])\n");
  for (int i=0; i<N; i++) {
    fprintf(fp, "%d\t%.e\t%.e\t%.e\t%.e\n", i, creal(psia[i]), cimag(psia[i]), creal(psib[i]), cimag(psib[i]));
  }
  fprintf(
    fp, "Coefficients:\n" "alpha = %.e + %.e * I\n" "beta = %.e + %.e * I\n",
    creal(alpha), cimag(alpha), creal(beta), cimag(beta)
  );
  fprintf(
    fp, "norm(delta) = norm( ||int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|| )\n"
    "norm(delta) = %.e", norm(delta, N)
  );
  fprintf("Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n");

  fclose(fp);

  /* Free allicated wavefunctions */
  free(psia);
  free(psib);
  free(left);
  free(right);
  free(delta);

  return 0;
}
