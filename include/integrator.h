#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "../kissfft/kiss_fft.h"


typedef struct {
  kiss_fft_cpx *cx_in, *cx_out;
  kiss_fft_cfg cfg, icfg;
}kissfft_struct;

/* error of O(tau^2) */
void euler_method(double complex *in,double tau);

/* Function used by UCN_method */
void operator(double complex *in, double complex *out);

/* error of O(tau^3) */
void UCN_method(double complex *in,double tau);

/*******************************************************************************
Strang splitting and necessary functions
*******************************************************************************/
/* Functions that allocate and free memory for any kissfft_strcut */
void alloc_kissfft(kissfft_struct *newfft_ptr, int N);
void free_kissfft(kissfft_struct *newfft_ptr);
/* Functions that allocate and free specific struct in integrator.c called kissfft */
void init_strangsplitting();
void finished_strangsplitting();
/* functions that convert data types for the FFT */
void double_to_kissfft_cpx(double complex* in, kiss_fft_cpx *out, int N);
// assigns real and imag of out by real and imag of in for all N values
void kissfft_cpx_to_double(kiss_fft_cpx *in, double complex* out, int N);
// assigns out[n] = in[n].r + i * in[n].i  for all N values
/* actual stranfsplitting method */
void strangsplitting_method(double complex *in, double tau);
// calucalates one time iteration step
// uses kiss_fft

#endif //INTEGRATOR_H
