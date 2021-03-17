#!/usr/bin/bash
time="1"
num="129"
nsteps=("20000" "15000" "12000" "10000")
integrator=("1" "2")
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

#!/usr/bin/bash
time="1"
num="129"
nsteps=("200000" "150000" "120000" "100000")
integrator=("0")
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
