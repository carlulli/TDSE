# TDSE
Code to numerically solve the time-dependent Schrödinger equation

# Building the code
- kissftt (https://github.com/mborgerding/kissfft) is needed as a directory as shown in this codes tree
- Note: building the library correctly didn't work for everyone and the "hacking" version might require deleting the downloaded kissftt directory and downloading it again. Thus, if you don't want to take any chances do the 2. step right away 
  1. to build the kissfft static library run `make KISSFFT_DATAYPE=double KISSFFT_STATIC=1`(for linux/ unix) 
This should create a `libkissfft-double.a` (.a for mac and maybe linux?) file 
  2. if it doesn't, go into `kiss_fft.h` and change the default `kiss_fft_scalar` from float to double in line 83 (unfortunatly for me this only worked before the first make, so i head to reclone the whole library)

- run `make` from TDSE directory to create inttest_linearity.exe, inttest_analytical.exe and inttest_unitarity.exe

# Tests

So far both these tests take slightly different inputs. (this will be changed)
- `inttest-linearity.exe [N] [time] [ntimestep] [integ_choice] [potential_choice]`
- `inttest-analyitcal.exe [N] [time] [ntimestep] [integ_choice] [potential_choice]`
- `inttest-unitarity.exe [N] [time] [ntimestep] [integ_choice] [potential_choice]`
- `test_gaussian_wf.exe [N] [time] [ntimestep] [integ_choice] [potential_choice] [mu] [dx] [dp]`

## Options
1 N = > 0 and odd
2 time = > 0
3 ntimestep = integer
4 integ_choice = 0 - Euler Method, 1 - Uintary Crank Nicolson Method, 2 - Strang Splitting Method
5 potential_choie = 0 - Zero Potential, 1 - Harmonic Oscillator, 2 - Well Potential, 3 - Wall Potential
6 mu = double
7 dx = double
8 dp = double
