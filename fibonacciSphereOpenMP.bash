#!/bin/sh
#BSUB -J Teddy-OpenMP
#BSUB -o project_output_file
#BSUB -e project_error_file
#BSUB -n 1
#BSUB -q ht-10g
#BSUB cwd /home/stoddard.t/FibonacciSphere/
work=/home/stoddard.t/FibonacciSphere/

cd $work

export OMP_NUM_THREADS=32
bin/fibonacciSphereOpenMP 1000000 10000000 32 123456789