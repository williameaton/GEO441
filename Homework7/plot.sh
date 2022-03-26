#!/bin/bash 

clean 
make clean 
make 



./xdiffusion

echo "Finished calculations. Starting plots:"
python3 plot.py