#!/bin/bash

#CODEDIR=/home/user/Desktop/FVM-Solver
#cd $CODEDIR

gfortran cavity-tdma.f90 -lblas -llapack # -mcmodel=large  -fbounds-check
./a.out
rm ./a.out