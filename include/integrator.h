#include <complex.h>



/* error of O(tau^2) */
void euler_method(double complex *in, double tau);

/* Function used by UCN_method */
void operator(double complex *in, double complex *out);

/* error of O(tau^3) */
void UCN_method(double complex *in, double tau);
/* error of O(tau^3) */
void splitting_method(double complex *in, double tau, int N);
