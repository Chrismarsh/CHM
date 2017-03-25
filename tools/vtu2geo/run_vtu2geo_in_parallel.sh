#!/bin/bash

configlist=$1
Nthreads=$2

cat $configlist | parallel --citation -j $Nthreads

echo "Finished."
