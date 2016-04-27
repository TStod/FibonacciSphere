#!/bin/sh
#BSUB -J Teddy-OpenACC
#BSUB -o project_output_file
#BSUB -e project_error_file
#BSUB -n 1
#BSUB -q par-gpu-2
#BSUB -R "span[ptile=32]"
#BSUB cwd /home/stoddard.t/FibonacciSphere/
work=/home/stoddard.t/FibonacciSphere/

cd $work

bin/fibonacciSphereOpenACC 1000000 10000000 4 123456789