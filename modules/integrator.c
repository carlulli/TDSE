#include <complex.h>

#include "conjugategradient.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "geometry.h"

/* time step of the integration method */
static double time_step;

/* takes array and integration step , modifies the array input*/
void euler_method(double complex *in, double tau) {
  /* uses the geometry.h library to obtain the N number of lattice points */
  int N = get_N();
  double complex Hpsi[N];
  H(in,Hpsi);

  for(int i = 0; i < N; i++) {
    in[i] = in[i] - tau*Hpsi[i]*I;
  }
}



/* function to use with conjugate_gradient() , applies the operator (1 +t*t/4*H*H) psi and gives back out */
void operator(double complex *in, double complex *out) {
  int N = get_N();
  //double complex Hpsi[N];
  H(in,in);
  H(in,in);
  for(int i = 0; i < N; i ++) {
    out[i] = 1 + 1/4*time_step*time_step*in[i];
  }
}

/* uses conjugategradient to calculate inverse of the operator it needs */
void UCN_method(double complex *in,double tau) {
  //time_step = tau;
  int N = get_N();
  double complex Hnu[N],HHnu[N], nu[N];
  void (*op_ptr)(double complex *, double complex *);
  op_ptr = &operator;

  for(int i = 0; i < N; i++) {

  }

  conj_grad(in,nu,(*op_ptr));

  H(nu,Hnu);
  H(Hnu,HHnu);
  
  for(int i = 0; i < N; i++) {
    in[i] = nu[i] - tau*Hnu[i]*I -1/4*tau*tau*HHnu[i];
  }
}
