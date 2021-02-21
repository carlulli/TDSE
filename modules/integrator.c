#include <complex.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "integrator.h"
#include "conjugategradient.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "geometry.h"
#include "../kissfft/kiss_fft.h"
// in kiss_fft.h in line 83: changed default from float to double


/* time step of the integration method */
static double time_step;

void integrator(double complex* in, double tau, int integ_choice) {
    /***************************************************************
   function that calculates the time evolution of input wavefunciton
   for a chosen integrator method
   remeber to initialze and finish strangsplitting if used
   ****************************************************************/
  if (integ_choice == 0) {
    euler_method(in, tau);
    printf("...Integrator used: Euler Method!\n");
  }
  else if (integ_choice == 1) {
    UCN_method(in, tau);
    printf("...Integrator used: Unitary Crank Nicolson Method!\n");
  }
  else if (integ_choice == 2) {
    // if (ssmcount==NULL) {
    //   init_strangsplitting();
    //   ssmcount=1;
    // }
    strangsplitting_method(in, tau);
    printf("...Integrator used: Strang Splitting Method!\n");
  }
  else {
    printf("[inttest_linearity.c | integrator()] Error! Choice of integrator is out of range!\n"
  "Remember: Integrator choice is 3rd input when calling inttest_linearity.c.\n "
  "Euler Method = 1, Unitary Crank-Nicolson Method = 2, Strang Splitting Method = 2\n" );

    exit(-1);
  }
}


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
  double complex Hin[N], HHin[N];
  H(in,Hin);
  H(Hin,HHin);
  for(int i = 0; i < N; i ++) {
    out[i] = in[i] + 1/4*time_step*time_step*HHin[i];
  }
}

/* uses conjugategradient to calculate inverse of the operator it needs */
void UCN_method(double complex *in, double tau) {
  //sets time_step
  time_step = tau;
  int N = get_N();
  double complex Hnu[N],HHnu[N], nu[N];
  void (*op_ptr)(double complex *, double complex *);
  op_ptr = &operator;

  conj_grad(in,nu,(*op_ptr));

  H(nu,Hnu);
  H(Hnu,HHnu);

  for(int i = 0; i < N; i++) {
    in[i] = nu[i] - tau*Hnu[i]*I -1/4*tau*tau*HHnu[i];
  }
}


/*******************************************************************************
Strang Splitting method
6 Steps

Using Kiss_FFT as FFT
https://github.com/mborgerding/kissfft
downloaded from github and ran
make KISSFFT_DATATYPE=double
to install
needed one more package: libpng
i got this through homebrew (for apple)
dont know if really necessary

0. Before first iteration.
initialize strangsplitting to initialize the kiss_fft parameters
with init_strangsplitting

1. calculates eta_q(n) = exp(-i*tau*V(n)/2)psi_q(n)
    V(n) is a scalar field
    n = {0, ... ,N-1} -> N

2. extending eta_q to dim 2N+2 to eta_ext_q
    but eta_ext_q has to be periodic woth period 2N+2 -> eta(n+(2N+2))=eta(n)
          (not really necessary as extension doesnt go further than 2N+2)
    and antisymmetric eta(2N-n)=-eta(n)
      -> this way Dirichlet bc are considered implicitely
    algorithm: for n<2N+2
      if n<N: eta_ext_q(n)=eta_q(n)
      if else n>N (maybe || n<2N+1): eta_ext_q=-eta_q(2N-n)
      else: eta_ext_q=0
  Q: is the dimension thing considere right?
        eta_ext_q has dim 2N+1. from n=0 to 2N with eta_ext_q(2N/2)=0
        additionally form Dbc: eta_ext_q(-1)=eta_ext_q(2N+1)=0

3. use Kiss_FFT for forward DFT
    DFT(eta_ext_q) calculates sum_n=0^2N+1 eta_ext_q * e_k(n) = eta_hat_q
    with e_k(n)=sin(pi(n+1)k/(N+1))
    Kiss_FFT calculates without norm
    (doesnt matter after reverse DFT) or yes.
      specifics on how to use Kiss_FFT below

4. calculates chi_hat_q(k) = Noramlization * exp((i*tau/2m)*laplace(k))eta_hat_q
    laplace(k) is a diagonal matrix with eigenvalues of laplace op.
    ev(k) = -4sin^2(pi*k/2(N+1))
    for k = 0,1,...,2N+1 (len(2N+2))
    normierung von davor mit drauf multiplizieren

5. calculates reverse DFT^-1
    chi_q = DFT-^1(chi_hat_q) for 2N+2

6. calculate psi_q+1 for only n = 0,1,...,N-1
    psi_q+1 = exp(-i*tau*V(n)/2)chi_q(n)

7. After last interation free kiss_fft parameters
*******************************************************************************/
void strangsplitting_method(double complex *in, double tau) {
  /* in = psi_q and out = psi_q+1 */
  int N = get_N();
  double mass;
  mass = get_m(); // function that needs to be defined

  double *V;
  V = (double*) malloc(sizeof(double)*N);
  assert(V!=0);

  /* wavefunctions and configuration for fft */
  kiss_fft_cpx *cx_in, *cx_out;
  kiss_fft_cfg cfg, icfg;

  cx_in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx)*(2*N+2));
  cx_out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx)*(2*N+2));
  cfg = kiss_fft_alloc( N, 0, 0, 0);
  icfg = kiss_fft_alloc( N, 1, 0, 0);
  assert(cx_in!=NULL); // checks if condition is true, otherwise assert
  assert(cx_out!=NULL);

  /* 1. part */
  double complex exparg1, rdummy1, idummy1;
  for (int n=0; n<N; n++) {
    V[n] = return_V(n);
    exparg1 = cos(0.5*tau*V[n])-I*sin(0.5*tau*V[n]); // is exp(-I*0.5*V[n]);
    // multply_dcx_element(in[n], exparg1, in[n]);
    rdummy1 = creal(in[n])*creal(exparg1)-cimag(in[n])*cimag(exparg1);
		idummy1 = creal(in[n])*cimag(exparg1)+cimag(in[n])*creal(exparg1);
		in[n] = rdummy1 + idummy1 * I; // c = dummy should be equivalent
  }

  /* 2. part */
  for (int n=0; n<2*N+2; n++) {
    if (n>0 && n<N+1) {
      cx_in[n].r = creal(in[n-1]);
      cx_in[n].i = cimag(in[n-1]);
    }
    else if (n>N+1) {
      cx_in[n].r = (-1)*cx_in[(2*N+2)-n].r;
      cx_in[n].i = (-1)*cx_in[(2*N+2)-n].i;
  }
    else {
      cx_in[n].r = 0.0;
      cx_in[n].i = 0.0;
    }
  }

  /* 3. part */
  if (cx_in != NULL) {  kiss_fft(cfg, cx_in, cx_out); }
  else {
    printf("[integrator.c | strangsplitting_method()] ERROR! FFT Plan not prepared.\n");
    exit(0);
  }

  /* 4. part */
  double complex arg,  exparg2, rdummy2, idummy2; // multiple dummys to be sure (even if more memory allocated)
  for (int k=0; k<2*N+2; k++) {
    arg = ((-2)*sin(M_PI*k/(2*N+2))*sin(M_PI*k/(2*N+2))*tau) / mass; // without I
    exparg2 = ( cos(arg) + I * sin(arg) ) / (2*N+2); // 1 / (2*N+2) for Noramlization (tested in helloworld program)
    // can't use multply_dcx_element as cx_in[k] is type kiss_fft_cpx not double complex
    // but cx_in[k].r and cx_in[k].i are type double
    rdummy2 = cx_out[k].r * creal(exparg2) - cx_out[k].i * cimag(exparg2);
		idummy2 = cx_out[k].r * cimag(exparg2) + cx_out[k].i * creal(exparg2);
		cx_out[k].r = rdummy2;
    cx_out[k].i = idummy2;
  }

  /* 5. part */
  kiss_fft(icfg, cx_out, cx_in);

  /* 6. part */
  // only look at N (or N+1?) values of chi_q with and "moving 1 step back to -1"
  // should it be V[n] or V[n+1]?
  double complex exparg3, rdummy3, idummy3;
  for (int n=0; n<N; n++) {
    V[n]=return_V(n);
    exparg3 = cos(0.5*tau*V[n])-I*sin(0.5*tau*V[n]); // is exp(-I*0.5*V[n]);
    rdummy3 = cx_in[n].r * creal(exparg3) - cx_in[n].i * cimag(exparg3);
		idummy3 = cx_in[n].r * cimag(exparg3) + cx_in[n].i * creal(exparg3);
		in[n] = rdummy3 + idummy3 * I; // in[n] = dummy3[n] should be equivalent
  }
  /* Free allocated memory */
  free(cx_in);
  free(cx_out);
  kiss_fft_free(cfg);
  kiss_fft_free(icfg);

  printf("... Strang Splitting Method ran once!\n");
}
