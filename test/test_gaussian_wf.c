//gaussian_wavefunction.case
//goal : test that set_gaussian_wavefunction correctly generates it


#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "wavefunction.h"

#define _FILE_NAME_ "test/test_gaussion_wf.c"

/*********************************************
1) Generates a gaussian wavepacket with a set of hardcored parameters
2) Print them on terminal and checks norm;
3) Can eventually be printed on file and plotted with python
***********************************************/

int main(int argc, char *argv[]) {

  /* Here the parameters are passed in the set_params function, which checks the validity */
  /* some params are not used and not called (get function) */
  set_params(argc, (char**) argv);

  int N = get_N();
  double mu = atof(argv[6]);
  double dx = atof(argv[7]);
  double dp = atof(argv[8]);

  /* wavefunction */
  double complex *psi;
  psi = (double complex*) malloc(sizeof(double complex)*N);

  set_gaussian_wavefunction(psi,mu,dx,dp,N);

  for(int i = 0; i < N; i++) {
  printf("psi[i] = %f  (I) %f\n", creal(psi[i]), cimag(psi[i]));
  }
  printf("norm = %f",norm(psi,N));

  return 0;
}
