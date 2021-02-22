#!/usr/bin/bash
time="1"
num=("149" "129" "55" "25")
nsteps=("5" "10" "50" "100" "500" "1000")
test=("analytical")
integrator=("1" "2")
potential=("0")
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
          # ./inttest_unitarity.exe $n ${time} $nstep 1 $pot
          # ./inttest_unitarity.exe $n ${time} $nstep 2 $pot
        done
      done
    done
  done
done


# potential=("0" "3")

# test=("analytical" "linearity" "unitarity")
