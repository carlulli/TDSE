#include <complex.h>

void set_random_wavefunction(double complex* psi, int N);
/* sets the wavefunction but doesn't normalize it    */
void set_random_wavefunction_NN(double complex* psi, int N);

/* sets psi to a gaussian wavepacket with random mu and the desired variance */
void set_gaussian_wavefunction(double complex* psi, double var, int N);
