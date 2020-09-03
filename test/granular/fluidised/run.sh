#!/bin/bash

# Clean directory
./clean.sh

#Run simulation
cplexec -c 1 "./CFD_dummy_fluidised.py" -m 1 "lmp_cpl < lammps/fluidised.in" 