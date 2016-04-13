#!/bin/sh
#BSUB -J Teddy-MPI
#BSUB -o project_output_file
#BSUB -e project_error_file
#BSUB -n 32
#BSUB -q ser-par-10g-3
#BSUB cwd /home/stoddard.t/FibonacciSphere/
#BSUB -R "span[ptile=4]"
work=/home/stoddard.t/FibonacciSphere/

cd $work
tempfile1=hostlistrun
tempfile2=hostlist-tcp
echo $LSB_MCPU_HOSTS > $tempfile1
declare -a hosts
read -a hosts < ${tempfile1}
for ((i=0; i<${#hosts[@]}; i += 2)) ;
do
 HOST=${hosts[$i]}
 CORE=${hosts[(($i+1))]}
 echo $HOST:$CORE >> $tempfile2
done

mpirun -np 32 -prot -TCP -lsf ./bin/fibonacciSphereMPI 100 10000 123456789

rm $work/$tempfile1
rm $work/$tempfile2
 