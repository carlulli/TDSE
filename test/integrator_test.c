#include <stdio.h>
#include <stdlib.h>

#include "geometry.h"
#include "linearalgebra.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include "integrator.h"
#include "assert.h"
#include "conjugategradient.h"

/* in this routine the correct functioning of the various integration methods implemented */

int main  (int argc, char *argv[]) {
  init_strangsplitting();

  strangsplitting_finished();
  return 0;
}
