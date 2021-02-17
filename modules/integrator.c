#include <complex.h>
#include <math.h>

#include "conjugategradient.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "geometry.h"
#include "kiss_fft.h"

/* time step of the integration method */
static double time_step;
static int initcount; // should be null pointer???!!!

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

7. After last interation free kiss_fft parameters with
    strangsplitting_finished
*******************************************************************************/


/* initialize the kiss_fftw parameters necessary for strang splitting */
void init_strangsplitting() {
  int N = get_N();
  kiss_fft_cpx *cx_in, *cx_out;
  cx_in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx)*N); // was malloc(sizeof(kiss_fft_cpx*N)) which should not
  kiss_fft_cfg cfg = kiss_fft_alloc( N, 0, NULL, NULL);
  kiss_fft_cfg icfg = kiss_fft_alloc( N, 1, NULL, NULL);
  initcount = 1; // or something with the ponter?
}

void strangsplitting_finished() {
  // maybe also free cx_in and cx_out???
  kiss_fft_free(cfg);
  kiss_fft_free(icfg);
  initcount = NULL;
}


void strangsplitting_method(double complex *in, double complex *out, double tau) {
  /* in = psi_q and out = psi_q+1 */

  int N = get_N();
  double mass = get_m(); // function that needs to be defined
  double complex eta_q[N], eta_ext_q[2*N+2], chi_q[N], chi_ext_q[2*N+2]; // safer and cleaner with dynamic allicating
  fftw_complex chi_hat_q;

  chi_hat_q = fftw_malloc(sizeof(fftw_complex)*(2*N+2));

  //
  // eta_ext_q will be tranformed to fftw_complex *
  // chi_q will be retransformed from fftw_complex *
  // chi_hat_q, eta_hat_q wil only be fftw_complex

  /* 1. part */
  for (int n=0; n<N; n++) {
    eta_q[n] = exp(-0.5*I*tau*get_V(n))*in[n]; // V(n) is some function that caluclates V(n) from hamiltonian module
  }

  /* 2. part */
  for (int n=0; n<2*N+2; n++) {
    if (n>0 || n<N+1) { eta_ext_q[n] = eta_q[n-1]; }
    if else (n>N+1) { eta_ext_q[n] = -eta_ext_q[(2*N+2)-n]; }
    else { eta_ext_q[n] = 0; }
  }

  /* 3. part */
  if (initcount != NULL) {
    for (int n=0; n<2*N+2; n++) {
        cx_in[n].r = creal(eta_ext_q[n]);
        cx_in[n].i = cimag(eta_ext_q[n]);
      }
    kiss_fft(cfg, cx_in, cx_out);
  }
  else {
    printf("[integrator.c | strangsplitting_method()] ERROR! FFTW Plan not prepared.\n",
    "init_strangsplitting was probably not called!\n");
    exit(0);
  }

/* 4. part */
  for (int k=0; k<2*N+2; k++) {
    fftw.in[k] = (double) (2*N+2)^(-1)*exp((I*tau/2*mass)*(-4)*sin^2(M_PI*k/(2*N+2)))*fftw.out[k];
  }

  /* 5. part */
  kiss_fft( icfg, cx_out, cx_in);

  for (int n=0; n<2*N+2; n++) {
    creal(chi_q[n]) = cx_in[n].r;
    cimag(chi_q[n]) = cx_in[n].i;
  }

  /* 6. part */
  for (int n=0; n<N+1; n++) {
    out[n] = exp(-0.5*I*tau*get_V(n))*chi_q[n+1];
  }

}
