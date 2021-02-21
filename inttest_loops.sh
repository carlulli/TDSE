#!/usr/bin/bash
# N TIME NSTEPS INTEGRATOR POTENTIAL

nsteps=("12","25","50","100")
integrators=("0","1","2")
filename=("unitarity","linearity", "analytical")
potential=("0","1")
time="100"
num="129"

for p in "${potential[@]}"
do
  for i in "${integrators[@]}"
  do
    for f in "${filename[@]}"
    do
      for n in "${nsteps[@]}"
      do
        echo ./inttest_${f}.exe ${num} ${time} ${n} ${i} ${p}
      done
    done
  done
done
