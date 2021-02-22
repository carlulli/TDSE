#!/usr/bin/bash
time="1"
num=("149" "129" "55" "25")
nsteps=("5" "10" "50" "100" "500" "1000")
test=("unitarity")
integrator=("1" "2")
potential=("0" "3")
for tes in "${test[@]}"
do
  for pot in "${potential[@]}"
  do
    for int in "${integrator[@]}"
    do
      for n in "${num[@]}"
      do
        for nstep in "${nsteps[@]}"
        do
          ./inttest_$tes.exe $n ${time} $nstep $int $pot
        done
      done
    done
  done
done
