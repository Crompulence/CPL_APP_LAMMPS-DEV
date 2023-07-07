#!/bin/bash

#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:10:0
#SBATCH --exclusive
#SBATCH --export=none
#SBATCH --account=y23

#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --nodes=2

# single thread export overriders any declaration in srun
export OMP_NUM_THREADS=1

module load other-software
module load cpl-lammps

# using your own installation: remove the last three lines and use these three 'module' lines and two 'source' lines instead remembering to update the path to your two SOURCEME.sh files
#module load openfoam/com/v2106
#module load cray-fftw
#module load cray-python
#source /work/y23/y23/gavincpl/cpl-library/SOURCEME.sh
#source /work/y23/y23/gavincpl/CPL_APP_OPENFOAM/SOURCEME.sh

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=1 lmp_cpl -in ./one_wall_imaginary_sliding_coupled.in :  --het-group=1 --nodes=1 --tasks-per-node=1 python ./python_dummy/CFD_test_vs_Couette.py
