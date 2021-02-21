#!/usr/bin/bash
time="1"
num=("129" "25")
nsteps=("10" "20" "50" "100" "1000" "5000" "10000")
test=("analytical" "linearity")
integrator=("0" "1" "2")
potential=("0" "1")
for nstep in "${nsteps[@]}"
do
  for int in "${integrator[@]}"
  do
    for pot in "${potential[@]}"
    do
      for tes in "${test[@]}"
      do
        for n in "${num[@]}"
        do
          ./inttest_$tes.exe $n ${time} $nstep $int $pot
          ./inttest_unitarity.exe $n ${time} $nstep 1 $pot
          ./inttest_unitarity.exe $n ${time} $nstep 2 $pot
        done
      done
    done
  done
done
