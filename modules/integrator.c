#include <complex.h>
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "geometry.h"


static double time_step;

void euler_method(double complex *in, double tau) {

  int N = get_N();
  double complex Hpsi[N];

  H(in,Hpsi);

  for(int i = 0; i < N; i++) {
    in[i] -= tau*Hpsi[i]*I;
  }
}

void operator(double complex *in, double complex *out) {
  int N = get_N();
  double complex Hpsi[N];
  H(in,Hpsi);
  H(Hpsi,Hpsi);
  for(int i = 0; i < N; i ++) {
    out[i] = 1 + 1/4*time_step*time_step*Hpsi[i];
  }
}

void UCN_method(double complex *in,double tau) {
  time_step = tau;
  int N = get_N();
  double complex Hnu[N],HHnu[N], nu[N];

  for(int i = 0; i < N; i++) {

  }

  conjugate_gradient(in,nu,&operator);

  H(nu,Hnu);
  H(Hnu,HHnu);
  for(int i = 0; i < N; i++) {
    in[i] = nu[i] - tau*Hnu[i]*I -1/4*tau*tau*HHnu[i];
  }
}
