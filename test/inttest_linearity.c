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

int main(int argc, char const *argv[]) {
/****************************************************************
argv[1] = N
argv[2] = time
argv[3] = nsteps
argv[4] = integrator_choice
argv[5] = potential

****************************************************************/


  // printf("This program calculates: maxdev = int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib) )\n"
  //   "For random wavefunctions psia and psib and random complex coefficients alpha and beta\n"
  //   "If the deviation is small (<10e-14?) the test was successful˜\n"
  // );


  srand(time(NULL)); // is called in beginning of the main.c
  set_params(argc, (char**) argv);
  N = get_N();
  double tau = get_time()/get_nsteps();
  int integrator_choice = get_integ_choice();
  double mass = 2.3512;
  set_kinetic_params(mass);
  int pot = get_pot_choice();
  set_potential(pot);
  // print_hamiltonian_info();

  double complex *psia, *psib, *psia_cp, *psib_cp;
  double complex alpha,  beta;
  double complex *left, *right;
  double maxdev, dev;

  /* dynamic memory allication (remember to free at the end) */
  psia = (double complex*) malloc(sizeof(double complex) * N);
  psib = (double complex*) malloc(sizeof(double complex) * N);
  psia_cp = (double complex*) malloc(sizeof(double complex) * N);
  psib_cp = (double complex*) malloc(sizeof(double complex) * N);
  left = (double complex*) malloc(sizeof(double complex) * N);
  right = (double complex*) malloc(sizeof(double complex) * N);
  assert(psia!=NULL); // check of condition is true else -> aboort?
  assert(psib!=NULL);
  assert(psia_cp!=NULL);
  assert(psib_cp!=NULL);
  assert(left!=NULL);
  assert(right!=NULL);

  /* generate random normalized wavefunctions psia and psib */
  set_random_wavefunction(psia, N);
  set_random_wavefunction(psib, N);

  // printf(
  //   "Linearity test: wavefuncitons:\n"
  //   "\tPsia[0] = %.e + %.e * i \n" "\tPsib[0] = %.e + %.e * i\n",
  //   creal(psia[0]), cimag(psia[0]), creal(psib[0]), cimag(psib[0])
  // );

  /* generate random coomplex coefficients */
  alpha = random_complex_coefficient();
  beta = random_complex_coefficient();

  // printf(
  //   "Linearity test: coefficients:\n"
  //   "\talpha= %.e + %.e * i \n" "\tPsib= %.e + %.e * i\n\n",
  //   creal(alpha), cimag(alpha), creal(beta), cimag(beta)
  // );

  /* Calculate int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib) */
  /* Left hand side */
  // if (integrator_choice==2) {init_strangsplitting();}

  for(int n = 0; n < N; n++) {
    left[n] = alpha*psia[n] + beta*psib[n];
  }
  copy_wf(psia, psia_cp, N);
  copy_wf(psib, psib_cp, N);

  integrator(left, tau, integrator_choice);
  // printf("DEBUGGING after int use\n");
  // exit(-1);
  /* Right hand side */
  integrator(psia, tau, integrator_choice);
  integrator(psib, tau, integrator_choice);
  // if (integrator_choice==2) {finished_strangsplitting();}
  /* ...and difference */
  maxdev = 0.0;
  for(int n = 0; n < N; n++) {
    right[n] = alpha*psia[n] + beta*psib[n];
    // printf("left[%d]= %.4e + %.4e * i\t\t right[%d]= %.4e + %.4e * i\n", n,creal(left[n]), cimag(left[n]), n,creal(right[n]), cimag(right[n]));
    dev = cabs(left[n] - right[n]);
    if(dev>maxdev)  {maxdev=dev;}
  }

  // printf(
  //   "\nMaximum Deviation is=\t%.6e\n"
  //   "Tolerances for to check for success:\t" "\teuler method = \t" "\tUCM method = \t" "\tstrang splittin method = 10e^-16\n", maxdev);
  printf("LINTEST with NUM=%d, time=%f, nsteps=%d, integ_choice=%d, pot_choice=%d\n", N, get_time(), get_nsteps(), integrator_choice, pot);
  printf("Maximum Deviation is =\t%.12e\n", maxdev);

  /* Printing test infos to text file */
  FILE *fp;
  int namesize = 45;
  for (int i=1; i<=5; i++) { namesize += strlen(argv[i]); }
  char filename[namesize];

  snprintf(
    filename, sizeof(filename),
    "data/int_lin_test_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4] ,argv[5]);

  fp = fopen(filename, "w");
  fprintf(fp, "n\tREAL(psia_cp[n]))\tIMAG(psia_cp[n]))\tREAL(psib_cp[n]))\tIMAG(psib_cp[n]))\n");
  for (int i=0; i<N; i++) {
    fprintf(fp, "%d\t%.6e\t%.6e\t%.6e\t%.6e\n", i, creal(psia_cp[i]), cimag(psia_cp[i]), creal(psib_cp[i]), cimag(psib_cp[i]));
  }
  fprintf(fp, "\nn\tREAL(int(psia[n]))\tIMAG(int(psia[n]))\tREAL(int(psib[n]))\tIMAG(int(psib[n]))\tREAL(alpha)\tIMAG(alpha)\tREAL(beta)\tIMAG(beta)\n");
  for (int i=0; i<N; i++) {
    fprintf(fp, "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", i, creal(psia[i]), cimag(psia[i]), creal(psib[i]), cimag(psib[i]),creal(alpha), cimag(alpha), creal(beta), cimag(beta));
  }
  fprintf(fp, "\nn\tREAL(int(left[n]))\tIMAG(int(left[n]))\tREAL(int(right[n]))\tIMAG(int(right[n]))\n");
  for (int i=0; i<N; i++) {
    fprintf(fp, "%d\t%.6e\t%.6e\t%.6e\t%.6e\n", i, creal(left[i]), cimag(left[i]), creal(right[i]), cimag(right[i]));
  }
  if (pot == 0) {fprintf(fp, "Potential used:\tZERO potential\n");}
  else if (pot == 1) {fprintf(fp, "Potential used:\tHARMONIC potential\n");}
  else if (pot == 2) {fprintf(fp, "Potential used:\tWELL potential\n");}
  else if (pot == 3) {fprintf(fp, "Potential used:\tWALL potential\n");}
  else {fprintf(fp, "Potential used:\tERROR\n");}

  // fprintf(
  //   fp, "norm(delta) = norm( ||int(alpha*psia + beta*psib) - (alpha*int(psia) + beta*int(psib)|| )\n"
  //   "norm(delta) = %.e\n", norm(delta, N)
  // );
  fprintf(fp, "MAXIMUM DEVIATION =\t%.16e\n" "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n", maxdev);

  fclose(fp);

  /* Free allicated wavefunctions */
  free(psia);
  free(psib);
  free(psia_cp);
  free(psib_cp);
  free(left);
  free(right);

  // printf("\n Linearity Test FINISHED! \n\n");

  return 0;
}
