#!/bin/bash
rm -f log.lammps &> /dev/null
rm -f lammps_forces.dat &> /dev/null
rm -f debug.vels &> /dev/null
rm -r cpl/map_CFD cpl/map_MD &> /dev/null
rm -f cpl/coupler_header &> /dev/null
rm -f cfd_forces.dat &> /dev/null
