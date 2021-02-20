//main.c
// goal: simulate the time evolution of a gaussian wave packet

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"


#define _FILE_NAME_ "test/inttest_convergence.c"
static int N = 0;

int main(int argc, char const *argv[]) {
/****************************************************************
argv[1] = N
argv[2] = time
argv[3] = nsteps
argv[4] = integrator_choice
argv[5] = potential
argv[6] = mu
argv[7] = dx
argv[8] = dp
****************************************************************/
/* Here the parameters are passed in the set_params function, which checks the validity */
set_params(argc, (char**) argv);
/* mass is always hardcoded */
double mass = 2.3512;
N = get_N();
double time = get_time();
double nsteps = get_nsteps();
double tau = time/nsteps;
int integrator_choice = get_integ_choice();
int pot_choice = get_pot_choice();
double mu = atof(argv[6]);
double dx = atof(argv[7]);
double dp = atof(argv[8]);
set_kinetic_params(mass);
set_potential(pot_choice);
print_hamiltonian_info();

printf("The program will calculate the TDSE of an initial wave packet with the following parameters input by the user:\n\nLattice size N = %d\nTotal run time t = %f\nNumber of steps nsteps = %d\nChoice of integrator (Euler [1] UCN [2] Strang-Splitting[3]) = %f\nChoice of Potential = %f\nParameters of the gaussian wavepacket mu = %f\t dx = %f\t dp = %f\n",N,time,nsteps,integrator_choice,pot_choice,mu,dx,dp);

/* generates the gaussian wave packet */
double complex psi[N];
set_gaussian_wavefunction(psi,mu,dx,dp,N);
printf("\n**********************************\n");
printf("Computing time evolution and printing on text file\n");
/* Printing test infos to text file */
FILE *fp;
int namesize = 60;
for (int i=1; i<=8; i++) { namesize += strlen(argv[i]); }
char filename[namesize];

snprintf(
  filename, sizeof(filename),
  "data/int_lin_test_%s_%s_%s_%s_%s_%s_%s_%s.txt", argv[1], argv[2], argv[3], argv[4] ,argv[5],argv[6],argv[7]);

fp = fopen(filename, "w");
fprintf(fp, "n\tREAL(psi[n])\tIMAG(psi[n])\ttau\taverx\tdeltax\taverp\tdeltap\n");
for (int j = 0; i < nsteps; j++) {
  integrator(psi,tau,integrator_choice);
  for (int i = 0; i < N; i++) {
  fprintf(fp, "%.e\t%.e\t%.e\t%.e\t%.e\t%.e\t%.e\n", creal(psi[i]),cimag(psi[i]),tau*j,get_avgx(psi),get_deltax(psi),get_avgp(psi),get_deltap(psi));
}
fclose(fp);


return 0;
}
