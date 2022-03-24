#!/bin/bash 

make clean 
make 

# Remove any time arrays in the snapshots: 
rm snapshots/time

./xdiffusion

python3 plot.py