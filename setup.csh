#!/bin/csh

rm jktebop_lib.so
gfortran jktebop_lib.f -ffree-form -fpic -shared -o jktebop_lib.so

cd oblateTransit
make clean
make
