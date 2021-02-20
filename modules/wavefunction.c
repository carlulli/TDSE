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


void set_gaussian_wavefunction(double complex* psi, double mu, double dx,double dp, int N) {

  for(int i = 0; i < N; i++) {
    psi[i] = /*1/(sqrt(2*M_PI)*sigma)*/cexp(-1/4*(i-mu)*(i-mu)/(dx*dx) + I*i*dp);

  }
  double normalizer = norm(psi,N);
  for (int i = 0; i < N; i ++) {
    psi[i]/= normalizer;
  }
}

double complex random_complex_coefficient() {
  return 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
}
