#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

void set_random_wavefunction(double complex* psi, int N);

/* sets a gaussian wave packet with mu mean value and sigma standard deviation */
void set_gaussian_wavefunction(double complex* psi, double mu, double sigma, int N);

/* set a non normalized wavefunction  */
void set_random_wavefunction_NN(double complex* psi, int N);

double complex random_complex_coefficient();
/* returns a random complex number */



#endif //WAVEFUNCTION_H
