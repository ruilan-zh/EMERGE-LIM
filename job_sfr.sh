#!/bin/sh

today=`date '+%Y-%m-%d'`
nohup ./sfr > ./tmp/out_${today}_sfr.log 2> ./tmp/err_${today}.log &
