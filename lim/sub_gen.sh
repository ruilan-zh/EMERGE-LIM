#!/bin/sh
set -e

pixlen=6.2 #$1 #arcseconds
Nmeshz=1 #$2
boxsize=500 #Mpc/h


dir=./figures_lambda=1.5_sphereX_pixlen=${pixlen}

names="Ha \
	OIII \
	OII \
	noise"

sum_indices="0 \
	1 \
	2 \
	NLINES"

make clean -f makefile-generate_stats
make -f makefile-generate_stats
./generate_stats ${dir} ${pixlen} ${boxsize} ${Nmeshz}
