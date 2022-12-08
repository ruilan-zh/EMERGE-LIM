#!/bin/sh

today=`date '+%Y-%m-%d'`
pixlen=6.2
nohup ./sub_gen.sh > ./tmp/out_${today}_sub_gen_${pixlen}.log 2> ./tmp/err_${today}_sub_gen_${pixlen}.log &
#nohup ./sub_gen.sh > ./tmp/out_${today}_sub_gen.log 2> ./tmp/err_${today}.log &
