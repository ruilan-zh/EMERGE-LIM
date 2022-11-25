#!/bin/sh
set -e

boxsize=250 #Mpc/h

lambda_obs=1.5 #micrometres

z_Ha=1.3
z_O3=2.0
z_O2=3.0

make clean
make -f makefile-luminosityfunc
echo "\n"
./luminosity_function ${lambda_obs} ${boxsize}

linenames="Ha \
	OIII \
	OII \
	noise"
names="Ha \
	OIII \
	noise"



