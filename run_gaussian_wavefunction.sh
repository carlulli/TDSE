#!/usr/bin/bash
gcc -Wall -o gauss.exe -I include  modules/geometry.c modules/linearalgebra.c modules/wavefunction.c test/gaussian_wavefunction.c -lm
