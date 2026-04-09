#!/bin/bash
#PBS -N Salmon_RT_Response_highaccuracy
#PBS -l nodes=1:ppn=32
#PBS -l walltime=48:00:00

ulimit -s unlimited

module purge
module load SALMON/v.2.0.1/openmpi-3.1.3_intel-2019.3.199

cd $PBS_O_WORKDIR

mpirun -np ${PBS_NUM_PPN} salmon < ./Si_rt_response.inp > RT_response.log
