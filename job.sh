#!/bin/sh

today=`date '+%Y-%m-%d'`
nohup python3 sfr_cent.py > ./tmp/out_${today}_sfr_mean.log 2> ./tmp/err_${today}.log &
