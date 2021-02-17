#include <complex.h>


/* error of O(tau^2) */
void euler_method(double complex *in,double tau);

/* Function used by UCN_method */
void operator(double complex *in, double complex *out);

/* error of O(tau^3) */
void UCN_method(double complex *in,double tau);

void splitting_method(double complex *in,double tau, int N);

/* Strang Splitting method */
void init_strangsplitting();
// allocates memory for fft parameters
// call once before starting interations of q

void strangsplitting_finished();
// frees memory for FFT
// use after finishing interations

void double_to_kissfft_cpx(double complex* in, kiss_fft_cpx *out, int N);
// assigns real and imag of out by real and imag of in for all N values

void kissfft_cpx_to_double(kiss_fft_cpx *in, double complex* out, int N);
// assigns out[n] = in[n].r + i * in[n].i  for all N values

void strangsplitting_method(double complex *in, double tau);
// calucalates one time iteration step
// uses kiss_fft
