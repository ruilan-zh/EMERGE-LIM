#!/bin/sh

today=`date '+%Y-%m-%d'`
nohup python3 sfr_mean.py > ./tmp/plot_data/out_${today}_sfr_mean.log 2> ./tmp/plot_data/err_${today}_sfr_mean.log &
nohup python3 sfrd.py > ./tmp/plot_data/out_${today}_sfrd.log 2> ./tmp/plot_data/err_${today}_sfrd.log &
nohup python3 sfrd_bins.py > ./tmp/plot_data/out_${today}_sfrd_bins.log 2> ./tmp/plot_data/err_${today}_sfrd_bins.log &
