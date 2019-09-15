#!/bin/bash

# Clean directory
./clean.sh

#Run simulation
cplexec -c 1 "./CFD_dummy_suzuki.py" -m 1 "lmp_cpl < lammps/suzuki.in" 
# cplexec -c 1 "./CFD_dummy_suzuki.py" -m 1 "lmp_cpl < lammps/suzuki.in" 