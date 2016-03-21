#!/bin/sh
#BSUB -J Teddy-Serial
#BSUB -o project_output_file
#BSUB -e project_error_file
#BSUB -n 1
#BSUB -q ht-10g
#BSUB cwd /home/stoddard.t/FibonacciSphere/
work=/home/stoddard.t/FibonacciSphere/

cd $work
bin/fibonacciSphere 100000 1000000 123456789