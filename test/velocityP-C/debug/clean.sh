#!/bin/bash
rm -f log.lammps &> /dev/null
rm -f lammps_vels.dat &> /dev/null
rm -f debug.vels &> /dev/null
rm -r cpl/map_CFD cpl/map_MD &> /dev/null
rm -f cpl/coupler_header &> /dev/null
rm -f cfd_vels[0-9]*.dat &> /dev/null
