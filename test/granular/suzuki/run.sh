#!/bin/bash

# Clean directory
./clean.sh

#Run simulation
cplexec -c 1 "./CFD_dummy_column.py" -m 1 "lmp_cpl < lammps/column.in" 
# cplexec -c 1 "./CFD_dummy_column.py" -m 1 "lmp_cpl < lammps/column.in" 