#!/bin/bash

#CODEDIR=/home/tushar/Desktop/CFD/fvm-codes
#cd $CODEDIR

gfortran cavity-tdma.f90 -lblas -llapack # -mcmodel=large  -fbounds-check
./a.out
rm ./a.out
