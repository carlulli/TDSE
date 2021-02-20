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
// static int *initcount=NULL; // should be null pointer???!!!
static kissfft_struct *kissfft=NULL;

void integrator(double complex* in, double tau, int integ_choice) {
    /***************************************************************
   function that calculates the time evolution of input wavefunciton
   for a chosen integrator method

   remeber to initialze and finish strangsplitting if used
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
    // if (ssmcount==NULL) {
    //   init_strangsplitting();
    //   ssmcount=1;
    // }
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

7. After last interation free kiss_fft parameters with
    strangsplitting_finished
*******************************************************************************/
/* function to allocate memory for kissfft_struct elements */
void alloc_kissfft(kissfft_struct *newfft_ptr, int N) {
  newfft_ptr->cx_in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx)*N);
  newfft_ptr->cx_out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx)*N);
  newfft_ptr->cfg = kiss_fft_alloc( N, 0, NULL, NULL);
  newfft_ptr->icfg = kiss_fft_alloc( N, 1, NULL, NULL);
  assert(newfft_ptr->cx_in!=NULL); // checks if condition is true, otherwise assert
  assert(newfft_ptr->cx_out!=NULL);
}

void free_kissfft(kissfft_struct *newfft_ptr) {
  // maybe also free cx_in and cx_out???
  free(newfft_ptr->cx_in);
  free(newfft_ptr->cx_out);
  kiss_fft_free(newfft_ptr->cfg);
  kiss_fft_free(newfft_ptr->icfg);
}

/* initialize the kiss_fftw parameters necessary for strang splitting */
void init_strangsplitting() {
  if (kissfft == NULL) {
    int N = get_N();
    kissfft = malloc(sizeof(kissfft_struct));
    assert(kissfft!=NULL);
    alloc_kissfft(kissfft, 2*N+2);
  }
  else {
    printf("[integrator.c | init_strangsplitting()] ERROR! KissFFT Struct already allocated!\n");

    exit(-1);
  }
}

void finished_strangsplitting() {
  if (kissfft != NULL) {
    free_kissfft(kissfft);
    kissfft = NULL;
  }
  else {
    printf("[integrator.c | init_strangsplitting()] KissFFT struct already freed!!\n");

    exit(-1);
  }
}

/* functions to convert double complex to kiss_fft_cpx */
void double_to_kissfft_cpx(double complex *in, kiss_fft_cpx *out, int N) {
  for (int n=0; n<N; n++) {
      out[n].r = creal(in[n]);
      out[n].i = cimag(in[n]);
      printf("DEBUGGING WHAT\n");
      exit(-1);
    }
}

void kissfft_cpx_to_double(kiss_fft_cpx *in, double complex *out, int N) {
  for (int n=0; n<N; n++) {
    out[n] = (double) (in[n].r) + (double) (in[n].i) * I;
  }
}


void strangsplitting_method(double complex *in, double tau) {
  /* in = psi_q and out = psi_q+1 */
  // struct fft fft;
  int N = get_N();
  double mass;
  mass = get_m(); // function that needs to be defined
  double complex *eta_q, *eta_ext_q, *chi_q, *chi_hat_q; // safer and cleaner with dynamic allicating

  /* dynamic allication of wavefunctions */
  eta_q = (double complex*) malloc(sizeof(double complex)*N);
  eta_ext_q = (double complex*) malloc(sizeof(double complex)*2*N+2);
  chi_q = (double complex*) malloc(sizeof(double complex)*N);
  chi_hat_q = (double complex*) malloc(sizeof(double complex)*N);
  assert(eta_q!=NULL);
  assert(eta_ext_q!=NULL);
  assert(chi_q!=NULL);
  assert(chi_hat_q!=NULL);

  double *V;
  V = (double*) malloc(sizeof(double)*N);
  assert(V!=0);
  /* 1. part */
  for (int n=0; n<N; n++) {
    V[n] = return_V(n);
    eta_q[n] = exp(-0.5*I*tau*V[n])*in[n]; // V(n) is some function that caluclates V(n) from hamiltonian module
    printf("DEBUGGING integrator eta_q[%d]= %.12e + %.12e * i\n", n, creal(eta_q[n]), cimag(eta_q[n]));
  }
  // printf("DEBUGGING integrator eta_q[1]= %.12e + %.12e * i\n", creal(eta_q[1]), cimag(eta_q[1]));
  /* 2. part */
  for (int n=0; n<2*N+2; n++) {
    if (n>0 && n<N+1) {eta_ext_q[n] = eta_q[n-1];
      printf("DEBUGGING [inside extension loop] (first if) eta_ext_q[%d]= %.5e + %.5e * i\n", n, creal(eta_ext_q[n]), cimag(eta_ext_q[n]) );
    }
    else if (n>N+1) {eta_ext_q[n] = (-1)*eta_ext_q[(2*N+2)-n];
    printf("DEBUGGING [inside extension loop] (if else) eta_ext_q[%d]= %.5e + %.5e * i\n", n, creal(eta_ext_q[n]), cimag(eta_ext_q[n]) );
  }
    else {
      eta_ext_q[n] = 0.0 + I * 0.0;
      printf("DEBUGGING [inside extension loop] (else) eta_ext_q[%d]= %.5e + %.5e * i\n", n, creal(eta_ext_q[n]), cimag(eta_ext_q[n]) );
    }
  }

  printf("DEBUGGING integrator eta_ext_q[1]= %.12e + %.12e * i\n", creal(eta_ext_q[1]), cimag(eta_ext_q[1]));
  /* 3. part */
  if (kissfft != NULL) {
    double_to_kissfft_cpx(eta_ext_q, kissfft->cx_in, 2*N+2);

    }
  else {
    printf("[integrator.c | strangsplitting_method()] ERROR! FFT Plan not prepared.\n"
    "init_strangsplitting was probably not called!\n");
    exit(0);
  }
  printf("DEBUGGING WHAT WHAT\n");
  exit(-1);

  for (int n=0; n<2*N+2; n++) {
    printf("DEBUGGING integrator kissfft->cx_in[%d]= %.5e + %.5e * i\n", n, (double) kissfft->cx_in[n].r, (double) kissfft->cx_in[n].i);
  }
  // printf("DEBUGGING integrator kissfft->cx_in[1]= %.5e + %.5e * i\n", (double) kissfft->cx_in[1].r, (double) kissfft->cx_in[1].i);
  printf("DEBUGGING integrator kissfft->cx_in[1] without (double) = %.5e + %.5e * i\n", (double) kissfft->cx_in[1].r, (double) kissfft->cx_in[1].i);
  printf("DEBUGGING now FFT.\n");
  kiss_fft(kissfft->cfg, kissfft->cx_in, kissfft->cx_out);
  printf("DEBUGGING integrator kissfft->cx_in[1]= %.5e + %.5e * i\n", (double) kissfft->cx_in[1].r, (double) kissfft->cx_in[1].i);
  printf("DEBUGGING integrator kissfft->cx_out[1]= %.5e + %.5e * i\n", (double) kissfft->cx_out[1].r, (double) kissfft->cx_out[1].i);
  printf("DEBUGGING integrator kissfft->cx_out[1] without (double) = %.5e + %.5e * i\n", kissfft->cx_out[1].r, kissfft->cx_out[1].i);
  kissfft_cpx_to_double(kissfft->cx_out, chi_hat_q, 2*N+2);
  printf("DEBUGGING integrator chi_hat_q[1]= %.12e + %.12e * i\n", creal(chi_hat_q[1]), cimag(chi_hat_q[1]));
  for (int k=0; k<2*N+2; k++) {
    // chi_hat_q[k] *= 1./norm(chi_hat_q, 2*N+2);
    // cx_in[k] = (double) (2*N+2)^(-1)*exp((I*tau/2*mass)*(-4)*sin^2(M_PI*k/(2*N+2)))*cx_out[k]; //trouble with datatype??
    // kissfft->cx_in[k] = (kiss_fft_cpx) (2*N+2)^(-1)*exp((I*tau/(2*mass))*(-4)*sin(M_PI*k/(2*N+2))*sin(M_PI*k/(2*N+2)))*kissfft->cx_out[k];
    /* kissfft->cx_in[k] is of datatype kiss_fft_cpx while kissfft->cx_in[k].r is float (or hopefully if successful: double) so calculation is easier */
    chi_hat_q[k] *= (double) (1./(2*N+2)*exp((I*tau/2*mass)*(-4)*sin(M_PI*k/(2*N+2))*sin(M_PI*k/(2*N+2))));
    // kissfft->cx_in[k].r = (double) 1./(2*N+2)*exp((I*tau/(2*mass))*(-4)*sin(M_PI*k/(2*N+2))*sin(M_PI*k/(2*N+2))) * kissfft->cx_out[k].r;
    // kissfft->cx_in[k].i = (double) 1./(2*N+2)*exp((I*tau/(2*mass))*(-4)*sin(M_PI*k/(2*N+2))*sin(M_PI*k/(2*N+2))) * kissfft->cx_out[k].i;
  }
  printf("DEBUGGING integrator chi_hat_q[1]= %.6e + %.6e * i\n", creal(chi_hat_q[1]), cimag(chi_hat_q[1]));
  double_to_kissfft_cpx(chi_hat_q, kissfft->cx_in, 2*N+2);
  printf("DEBUGGING integrator kissfft->cx_in[1]= %.6e + %.6e * i\n", (double) kissfft->cx_in[1].r, (double) kissfft->cx_in[1].i);
  /* 5. part */
  printf("DEBUGGING now iFFT.\n");
  kiss_fft(kissfft->icfg, kissfft->cx_in, kissfft->cx_out);
  printf("DEBUGGING integrator kissfft->cx_out[1]= %.6e + %.6e * i\n", (double) kissfft->cx_out[1].r, (double) kissfft->cx_out[1].i);

  kissfft_cpx_to_double(kissfft->cx_out, chi_q, N);
  printf("DEBUGGING integrator chi_q[1]= %.6e + %.6e * i\n", creal(chi_q[1]), cimag(chi_q[1]));


  /* 6. part */
  // only look at N (or N+1?) values of chi_q with and "moving 1 step back to -1"
  for (int n=0; n<N+1; n++) {
    in[n] = exp(-0.5*I*tau*V[n])*chi_q[n+1];
  }

}
