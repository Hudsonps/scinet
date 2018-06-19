#!/bin/bash
#PBS -l nodes=4:ppn=8
#PBS -l walltime=1:00:00
#PBS -N jobassign8final
cd $PBS_O_WORKDIR
module load gcc/4.8.1 openmpi/gcc/1.6.4 rarray
for i in {1..32}
do
    time mpirun -np $i ./assign2
done
# My submission script.
