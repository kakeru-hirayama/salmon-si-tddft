#!/bin/bash
#PBS -N Salmon_RT_ex06
#PBS -l nodes=1:ppn=32
#PBS -l walltime=12:00:00

ulimit -s unlimited

module purge
module load SALMON/v.2.0.1/openmpi-3.1.3_intel-2019.3.199

cd $PBS_O_WORKDIR

# restart/ ディレクトリにGS計算の data_for_restart をコピーしてから実行
mpirun -np ${PBS_NUM_PPN} salmon < ./Si_rt_pulse.inp > Si_rt_pulse.out
