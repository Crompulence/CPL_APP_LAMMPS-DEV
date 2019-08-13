#!/bin/bash
rm -f log.lammps &> /dev/null
rm -f lammps_bc.dat &> /dev/null
rm -r cpl/map_CFD cpl/map_MD &> /dev/null
rm -f cpl/coupler_header &> /dev/null
rm -f cfd_bc[0-9]*.dat &> /dev/null
