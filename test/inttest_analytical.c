#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "wavefunction.h"
#include "integrator.h"
#include "geometry.h"
#include "linearalgebra.h"
#include "assert.h"
#include "hamiltonian.h"

#define _FILE_NAME_ "test/inttest_analytical.c"

static int N = 0;

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

create eigenvector phi of laplace as above
  -> for k = 0,...,N

for every k:
  calculate phi, int(phi), fofE
  for every n:
    calcualte dev=cabs(left[n] - fofE*phi[n]);
  printf k maxdev (out and to file )
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
  "Remember: Integrator choice is 3rd input when calling inttest_linearity.c.\n "
  "Euler Method = 1, Unitary Crank-Nicolson Method = 2, Strang Splitting Method = 2\n" );

    exit(-1);
  }
}

void set_sinwave(double complex *psi,int k)
{

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
  /****************************************************************
  argv[1] = N
  argv[2] = tau
  argv[3] = integrator_choice
  ****************************************************************/

    // assert(argc==5,_FILE_NAME_,"main","ERROR. Necessary number of input parameters 4!\n
    // Usage: inttest_analytical {N} {tau} {integrator_choice}\n
    // Remember: The executable Name is the first parameter.\n");

      /* HERE THE MASS IS HARDCODED */
    set_params(argc, argv);
    double mass = 2.3512;
    int N = get_N();
    double tau = atoi(argv[2]);
    int integrator_choice = atoi(argv[3]);
    int pot = 0;
    set_kinetic_params(mass);
    set_potential(pot);
    print_hamiltonian_info();

    printf("\nThis programs checks the eigenvalues and eigenvectors of the Hamiltonian\n"
     "with V=0. In this case one can calculate analytically that the eigenvectors\n"
     "are labeled by an integer index k. The k-th eigenvector is given by\n"
     "   psi(n) = sqrt(2/(N+1)) * sin(pi*(n+1)*k/(N+1))     for n=0,1,2,...,N-1\n"
     "and the corresponding eigenvalue is given by\n"
     "   E = -4 * sin( pi*k/(2*(N+1)) )^2\n\n"
     "For each k=0,1,2,...,N this program prints\n"
     "   maxdev = max_n | H.psi(n) - E psi(n) |\n"
     "The test passes is all these numbers are < 1e-15\n\n");

    double complex *phi, *left;
    // double complex *delta;
    double fofE, maxdev, dev;

    /* dynamic memory allication (remember to free at the end) */
    phi = (double complex*) malloc(sizeof(phi) * N);
    left = (double complex*) malloc(sizeof(double) * N);
    // delta = (double*) malloc(sizeof(delta) * N);

    /* set phi as eigenvector of laplace
    -> k=1
    */
    // for (int n=0; n<N; n++) {
    //   phi[n] = sqrt(2./(N+1))*sin(M_PI*(n+1)*k/(N+1));
    //   left[n] = phi[n]; // copy of phi
    // }
    // for (int n=0; n<2; n++) {
    //     printf("First three phi values:\t phi[%d] = %.e + %e * i\n", n, creal(phi[n]), cimag(phi[n]);
    // }

    FILE *fp;
    int namesize = 19;
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
    fprintf(fp, "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n");
    // fprintf(fp, "n\tREAL(phi[n])\tIMAG(phi[n])\tmaxdeviation\n");
    fprintf(fp, "k\tMAXDEVIATION\n");

    for (int k=0;k<=N;k++) {
          set_sinwave(phi,k);
          copy_wf(phi, left);
          integrator(left, tau, integrator_choice);
          // ev = 2.*sin((M_PI*k)/(2*(N+1)))*sin((M_PI*k)/(2*(N+1)))/mass;
          fofE = exp( I*tau*(-4)*sin(M_PI*k/(2*(N+1)))*sin(M_PI*k/(2*(N+1)))/mass );
          maxdev=0.0;
          for (int n=0;n<N;n++) {
             dev=cabs(left[n] - fofE*phi[n]);
             if(dev>maxdev) maxdev=dev;
             // fprintf(fp, "%d\t%.e\t%.e\t%.e\n", i, creal(psia[i]), cimag(psia[i]), creal(psib[i]), cimag(psib[i]));
          }
          printf("k= %d\tmaxdev= %.2e\n",k,maxdev);
          fprintf(fp, "%d\t%.e\n", k, maxdev);
       }


    // printf("Calcualted difference || integrator(phi) - fofE*phi || = %.e \n", norm(delta, N));
    // printf("Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splittin method = 10e^-16\n");


    fclose(fp);

    /* Free allicated wavefunctions */
    free(phi);
    free(left);
    // free(delta);

  return 0;
}
