#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <math.h>


#include "linearalgebra.h"

void set_random_wavefunction(double complex* psi, int N) {

  // srand(time(NULL)); is called in beginning of the main.c

  for(int  i = 0; i < N; i++) {
    psi[i] = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
  }
  double psi_norm = norm(psi,N);
  for(int  i = 0; i < N; i++) {
    psi[i]  = psi[i]/psi_norm;

  }
}

void set_random_wavefunction_NN(double complex* psi, int N) {

  for(int i = 0; i < N; i++) {
      psi[i] = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
    }
  }


void set_gaussian_wavefunction(double complex* psi, double mu, double sigma, int N) {


  for(int i = 0; i < N; i++) {
    psi[i] = 1/(sqrt(2*M_PI)*sigma)*cexp(-0.5*(i-mu)*(i-mu)/(sigma*sigma));

  }
  for (int i = 0; i < N; i ++) {
    psi[i]/= norm(psi,N);
  }
}
