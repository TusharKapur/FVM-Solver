#!/bin/bash

#CODEDIR=/home/tushar/Desktop/CFD/fvm-codes
#cd $CODEDIR
#rm output
gfortran cavity.f90 -lblas -llapack  -mcmodel=large  -fbounds-check
./a.out
rm ./a.out
