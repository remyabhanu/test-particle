#!/bin/bash
MODULE := weimFort

all:
	f2py --overwrite-signature w2k.f -m ${MODULE} -h w2k.pyf
	gfortran -w -O2 -fbacktrace -fno-automatic -fPIC -c w2k.f
	f2py --f77flags="-w" -c w2k.pyf w2k.f

clean:
	find . -type f -name "*.o" | xargs rm -f
	find . -type f -name "*.so" | xargs rm -f
	find  -type f -name "*.pyf" | xargs rm -f
