#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>

#include "wavefunction.h"
#include "integrator.h"
#include "geometry.h"
#include "linearalgebra.h"
// #include "assert.h"
#include "hamiltonian.h"

#define _FILE_NAME_ "test/inttest_analytical.c"

static int N = 0;
// static int *ssmcount=NULL;

/*******************************************************************************
// Linearity test

// goal:
//  1x time evolution for || f(H)psi - f(E)psi||
calculate this norm for 1 time step

Use V=0 -> eigenvlue and eigenvector of H=laplace are known
eigenvector(n) = sqrt(2/(N+1)) sin(pi*(n+1)*k/(N+1))
eigenvalue(k) = -4sin^2(pi*k/2(N+1))
k = 0,...,N ???
H = K + V = K = - 1/2m_hat * laplace
f(H) = exp(-i*tau_hat*H) -> integrator applied
f(E) = exp(i*tau_hat/2m_hat*eigenvalue)

create eigenvector psi of laplace as above
  -> for k = 0,...,N

for every k:
  calculate psi, int(psi), fofE
  for every n:
    calcualte dev=cabs(left[n] - fofE*psi[n]);
  printf k maxdev (out and to file )
*******************************************************************************/



void set_eigenfunction(double complex *psi,int k) {
  for(int n=0;n<N;n++)
    // psi[n]=sin((M_PI*((n+1)*k))/(N+1));
    psi[n] = sqrt(2./(N+1))*sin((M_PI*(n+1)*k)/(N+1));
}

void copy_wf(double complex *in, double complex *out) {
  for (int n=0; n<N; n++) {
    out[n] = in[n];
  }
}

int main(int argc, char const *argv[]) {

      /* HERE THE MASS IS HARDCODED */
    set_params(argc, (char**) argv);
    double mass = 2.3512;
    N = get_N();
    double tau = get_time()/get_nsteps();
    int integrator_choice = get_integ_choice();
    int pot = 0;
    set_kinetic_params(mass);
    set_potential(pot);
    print_hamiltonian_info();



    printf("\nThis programs checks the eigenvalues and eigenvectors of the Hamiltonian"
     "with V=0.\n In this case one can calculate analytically that the eigenvectors"
     "are labeled by an integer index k.\n The k-th eigenvector is given by\n"
     "   psi(n) = sqrt(2/(N+1)) * sin(pi*(n+1)*k/(N+1))     for n=0,1,2,...,N-1\n"
     "and the corresponding eigenvalue is given by\n"
     "   E = -4 * sin( pi*k/(2*(N+1)) )^2\n\n"
     "For each k=0,1,2,...,N this program prints\n"
     "   maxdev = max_n | H.psi(n) - E psi(n) |\n"
     "The test passes is all these numbers are < 1e-15\n\n"
   "Test prams are: N=%d, mass=%f, tau=%.e, integrator_choice=%d\n\n",
   N, mass, tau, integrator_choice);



    double complex *psi, *left;
    // double complex *delta;
    double fofE, maxdev, dev;

    /* dynamic memory allication (remember to free at the end) */
    psi = (double complex*) malloc(sizeof(double complex) * N);
    left = (double complex*) malloc(sizeof(double complex) * N);
    assert(left!=NULL);
    assert(psi!=NULL);
    // delta = (double*) malloc(sizeof(delta) * N);

    /* set psi as eigenvector of laplace
    -> k=1
    */
    // for (int n=0; n<N; n++) {
    //   psi[n] = sqrt(2./(N+1))*sin(M_PI*(n+1)*k/(N+1));
    //   left[n] = psi[n]; // copy of psi
    // }
    // for (int n=0; n<2; n++) {
    //     printf("First three psi values:\t psi[%d] = %.e + %e * i\n", n, creal(psi[n]), cimag(psi[n]);
    // }

    FILE *fp;
    int namesize = 30;
    for (int i=1; i<=3; i++) { namesize += strlen(argv[i]); }
    char filename[namesize];

    snprintf(
      filename, sizeof(filename),
      "data/int_anal_test_%s_%s_%s.txt", argv[1], argv[2], argv[3]
    );

    fp = fopen(filename, "w");
    /* file header  */
    if (pot == 0) {fprintf(fp, "Potential used:\tZERO potential\n");}
    else if (pot == 1) {fprintf(fp, "Potential used:\tHARMONIC potential\n");}
    else if (pot == 2) {fprintf(fp, "Potential used:\tWELL potential\n");}
    else if (pot == 3) {fprintf(fp, "Potential used:\tWALL potential\n");}
    else {fprintf(fp, "Potential used:\tERROR\n");}
    fprintf(fp, "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splitting method = 10e^-16\n");
    // fprintf(fp, "n\tREAL(psi[n])\tIMAG(psi[n])\tmaxdeviation\n");
    fprintf(fp, "k\tMAXDEVIATION\tn\treal(int(psi))\timag(int(psi))\treal(f(E)*psi)\treal(f(E)*psi)\n");

    if (integrator_choice==2) {init_strangsplitting();}
    for (int k=0;k<=N;k++) {
          set_eigenfunction(psi,k);
          for (int n=0; n<N; n++) {
            printf("DEBUGGING inttest_analytical psi[%d]= %.12e + %.12e * i\n", n, creal(psi[n]), cimag(psi[n]));
          }
          copy_wf(psi, left);
          integrator(left, tau, integrator_choice);
          // ev = 2.*sin((M_PI*k)/(2*(N+1)))*sin((M_PI*k)/(2*(N+1)))/mass;
          fofE = exp( (double) (I*tau*(-4)*sin(M_PI*k/(2*(N+1)))*sin(M_PI*k/(2*(N+1)))/mass) );
          printf("fofE=%.e\n", fofE);
          maxdev=0.0;
          printf("maxdev = %.e\n", maxdev);

          for (int n=0;n<N;n++) {
             dev=cabs(left[n] - fofE*psi[n]);
             if(dev>maxdev) maxdev=dev;
             fprintf(fp, "%d\t%.e\t%d\t%.e\t%.e\t%.e\t%.e\n", k, maxdev, n, creal(left[n]), cimag(left[n]), creal(fofE*psi[n]), cimag(fofE*psi[n]));
          }
          printf("Calculating...\tk= %d\tmaxdev= %.2e\n",k,maxdev);
          // fprintf(fp, "%d\t%.e\n", k, maxdev);
       }
      if (integrator_choice==2) {finished_strangsplitting();}



    // printf("Calcualted difference || integrator(psi) - fofE*psi || = %.e \n", norm(delta, N));
    // printf("Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n");
    fclose(fp);

    /* Free allicated wavefunctions */
    free(psi);
    free(left);
    // free(delta);

  return 0;
}
