#!/bin/bash

configlist=$1
Nthreads=$2

cat $configlist | parallel -j $Nthreads

echo "Finished."
