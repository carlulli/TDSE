#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <complex.h>


double complex scalar_product(double complex* a, double complex* b,int N);
/* takes two complex vectors and the lenght of the vectors
and retruns one complex number */

double norm(double complex* a,int N);
/*analogous but only computes the norm*/

/*Function for acting with matrix on vector*/
void multiplyvec(double complex **mat, double complex *vec, double complex *res);

/* Function for creating a semi definite positive random Matrix*/
void multAtimesv( double complex* in, double complex* out);

/* function for random vector */
void randvec(double complex* vec, int M);

/*******************************************************************************
function to correctly multiply two complex vectors c = a + b for all n<N
for only one compontent pass N=1
*******************************************************************************/
void multply_dcx_wf(double complex *a, double complex *b, double complex *c, int N);

/*******************************************************************************
function to correctly multiply two complex vectors c = a + b element wise
*******************************************************************************/
// void multply_dcx_element(double complex *a, double complex *b, double complex *c);

#endif // !LinearAlgebra_h
