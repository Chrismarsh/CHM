#!/bin/bash

# ncore to use

ncore=(1 2 4 8 16 32 48)

retry=1

printf "ncore,run,duration\n"
for i in "${ncore[@]}"
do
   t=9999
   # do whatever on $i
   for j in {1..5}; do 
   		 start=$SECONDS

   		 OMP_NUM_THREADS=$i LD_LIBRARY_PATH=~/build-release/lib/gsl/lib:$LD_LIBRARY_PATH ~/build-release/bin/Release/CHM -f Flood_0p25_nc.json &> /dev/null

   		 duration=$(( SECONDS - start ))
   		 printf "%i,%i,%f\n" "$i" "$j" "$duration"
   done; 


done