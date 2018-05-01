#!/bin/sh
set -e
if [ -f "gitlammps" ]; then
    echo "Using cached LAMMPS"
else
    case $1 in
        lammps) set -x;
            #Get LAMMPS
            git clone https://github.com/lammps/lammps.git gitlammps
        granlammps) set -x;
            #Get GranLAMMPS first as this requires password
            git clone https://edwardsmith999@bitbucket.org/granlammps/gitlammps.git gitlammps
            cd ./gitlammps
            git checkout common
            cd ../
    esac
fi

