# TDSE
Code to numerically solve the time-dependent Schrödinger equation

# Building the code
- kissftt (https://github.com/mborgerding/kissfft) is needed as a directory as shown in this codes tree
  - to build the kissfft static library run `make KISSFFT_DATAYPE=double KISSFFT_STATIC=1`(for linux/ unix) 
  - this should create a libkissfft-double.a (.a for mac and maybe linux?) file 
  - if it doesn't go into kiss_fft.h and change the default KISSFFT_SCALAR from float to double in line 83 (unfortunatly for me this only worked before the first make, so i head to reclone the whole library)

- run `make` from TDSE directory to create inttest_linearity.exe and inttest_analytical.exe

# Tests

So far both these tests take slightly different inputs. (this will be changed)
- `inttest-linearity.exe [N] [tau] [integ_choice] [potential_choice]`
- `inttest-analyitcal.exe [N] [tau] [integ_choice]`

Strang splitting method is integ_choice = 2
