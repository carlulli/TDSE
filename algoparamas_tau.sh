#!/usr/bin/bash

#Investigate the tau parameter dependence of the integration methods

time=10
nsteps=10
while [$nsteps < 1]
do
  echo $ ./main.exe 129 $time $nsteps 10 0 0 56 20 20
  ((nsteps -1 ))
done
echo Simulation Ended
