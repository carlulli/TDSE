#include <complex.h>

void euler_method(double complex *in,double tau);
void operator(double complex *in, double complex *out);
void UCN_method(double complex *in,double tau);
void splitting_method(double complex *in,double tau, int N);
