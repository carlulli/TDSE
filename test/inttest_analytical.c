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



void set_eigenfunction(double complex *psi, int k) {
  for(int n=0;n<N;n++)
    // psi[n]=sin((M_PI*((n+1)*k))/(N+1));
    psi[n] = sqrt(2./(N+1))*sin((M_PI*(n+1)*k)/(N+1));
}


/*******************************************************************************
Function you run in bash with:
$ ./inttest_analytical.exe [N] [TIME] [NSTEPS] [INTEGRATOR] [POTENTIAL]
*******************************************************************************/
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
     "   maxdev = max_n | int(H).psi(n) - int(EV).psi(n) |\n"
     "The test passes is all these numbers are < 1e-15\n\n"
   "Test prams are: N=%d, mass=%f, tau=%.e, integrator_choice=%d\n\n",
   N, mass, tau, integrator_choice);



    double complex *psi, *left, *right;
    double maxdev, dev;

    /* dynamic memory allication (remember to free at the end) */
    psi = (double complex*) malloc(sizeof(double complex) * N);
    left = (double complex*) malloc(sizeof(double complex) * N);
    right = (double complex*) malloc(sizeof(double complex) * N);
    assert(left!=NULL);
    assert(psi!=NULL);
    assert(right!=NULL);

    /* opening the file for the output */
    FILE *fp;
    int namesize = 50;
    for (int i=1; i<=5; i++) { namesize += strlen(argv[i]); }
    char filename[namesize];

    snprintf(
      filename, sizeof(filename),
      "data/int_anal_test_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4], argv[5]
    );

    fp = fopen(filename, "w");
    /* file header  */
    if (pot == 0) {fprintf(fp, "Potential used:\tZERO potential\n");}
    else if (pot == 1) {fprintf(fp, "Potential used:\tHARMONIC potential\n");}
    else if (pot == 2) {fprintf(fp, "Potential used:\tWELL potential\n");}
    else if (pot == 3) {fprintf(fp, "Potential used:\tWALL potential\n");}
    else {fprintf(fp, "Potential used:\tERROR\n");}
    fprintf(fp, "Tolerances for to check for success:\n" "\teuler method = \n" "\tUCM method = \n" "\tstrang splitting method = 10e^-16\n");
    fprintf(fp, "\nk\tDEVIATION\tn\treal(int(psi))\timag(int(psi))\treal(f(E)*psi)\timag(f(E)*psi)\tcreal(psi)\tcimag(psi)\n");

/*******************************************************************************
Actual running the test
*******************************************************************************/
    // if (integrator_choice==2) {init_strangsplitting();}
    double complex arg, fofE, rdummy, idummy;
    double finaldev=0.0;
    for (int k=0;k<=N;k++) {
          set_eigenfunction(psi,k);
          // for (int i=0;i<N;i++) { psi[i]=1.+I; }
          copy_wf(psi, left, N);
          integrator(left, tau, integrator_choice);

          /*
           * arg = H => K = -EV/2m = -( -4*sin^2( pi*k/(2(N+1)) )/2m )
           * EigenValue of laplace operator
           * the same for all integrators
           */
          arg = (2 * sin(M_PI*k/(2*(N+1))) * sin(M_PI*k/(2*(N+1)))) / mass;

          if (integrator_choice==0) { fofE = 1. - I*tau*arg; } // check complex number assigning
          else if (integrator_choice==1) { fofE = (1. - I*0.5*tau*arg) / (1. + I*0.5*tau*arg); }
          else if (integrator_choice==2) { fofE = cos(tau*arg) - I * sin(tau*arg); }

          // printf("fofE (%d) = %.4e + %.4e * I\n", integrator_choice , creal(fofE), cimag(fofE));
          maxdev=0.0;
          for (int n=0;n<N;n++) {
            rdummy = creal(fofE)*creal(psi[n])-cimag(fofE)*cimag(psi[n]);
        		idummy = creal(fofE)*cimag(psi[n])+cimag(fofE)*creal(psi[n]);
        		right[n] = rdummy + idummy * I;
            // printf("DEBUGGING 1 psi[%d] = %.4e + i * %.4e\n", n, creal(psi[n]), cimag(psi[n]));
            // printf("DEBUGGING 1 right[%d] = %.4e + i * %.4e\n", n, creal(right[n]), cimag(right[n]));
            // multply_dcx_element(&fofE, &psi[n], &right[n]); // doesnt work with the pointers correctly
            dev=cabs(left[n] - right[n]);
            if(dev>maxdev) maxdev=dev;
            fprintf(fp, "%d\t%.6e\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
            k, dev, n, creal(left[n]), cimag(left[n]), creal(right[n]), cimag(right[n]), creal(psi[n]), cimag(psi[n]) );
          }
          // printf("...Calculating\tk= %d\tmaxdev= %.6e\n",k,maxdev);
          if(maxdev>finaldev) finaldev=maxdev;
       }
      // if (integrator_choice==2) {finished_strangsplitting();}
    fprintf(fp, "\n\nFINAL MAXDEV = %.12e\t\n", finaldev);
    fclose(fp);
    printf( "\n\nFINAL MAXDEV = %.12e\t\n", finaldev);
    /* Free allicated wavefunctions */
    free(psi);
    free(left);
    free(right);

    printf("\n Analytical Test FINISHED! \n\n");

  return 0;
}
