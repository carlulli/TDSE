#!/usr/bin/bash
time="1"
num="129"
nsteps=("10" "20" "50" "100" "200" "1000")
integrator=("0" "1" "2")
potential="3"
#size of wall (N+1)/2 -- (N+5)/2  height  10
for nstep in "${nsteps[@]}"
do
  for int in "${integrator[@]}"
  do
    for pot in "${potential[@]}"
    do
      for n in "${num[@]}"
      do
        ./main.exe $n ${time} $nstep $int $pot 30 6 1
      done
    done
  done
done
