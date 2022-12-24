#!/bin/sh
set -e

boxsize=500 #Mpc/h
pixlen=100 #arcsec
I_sigma=1.912202

lambda_obs=1.5 #micrometres

z_Ha=1.47
#z_O3=2.0
#z_O2=3.0

idir=./figures_lambda=${lambda_obs}_sphereX_pixlen=${pixlen}
odir=${idir}/plots

mkdir -p $idir
mkdir -p $odir

make clean
make
echo "\n"
./lim_pinocchio ${lambda_obs} ${I_sigma} ${boxsize} ${pixlen} ${idir}

linenames="Ha \
	OIII \
	OII \
	noise"
names="Ha \
	OIII \
	noise"

#python3 projmap.py ${linenames} \
#	-idir:${idir} \
#	-fout:${odir}/projmap.png

#python3 projmap2.py ${names} \
#	-idir:${idir} \
#	-fout:${odir}/projmap2.png

#python3 stats.py ${linenames} \
#	-idir:${idir} \
#	-odir:${odir}



