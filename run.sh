#!/bin/bash

CODEDIR=/mnt/c/Users/Tushar/Desktop/CFD/fvm-codes
cd $CODEDIR

gfortran cavity.f90 -lblas -llapack # -mcmodel=large  -fbounds-check
./a.out
rm ./a.out
