#!/bin/sh
set -e
if [ -d "gitlammps" ]; then
    echo "Using cached LAMMPS"
else
    case $1 in
      lammps) set -x;
        #Get LAMMPS
        git clone https://github.com/lammps/lammps.git gitlammps
        cd gitlammps
        git checkout stable
        cd ../;;
      granlammps) set -x;
        #Get GranLAMMPS first as this requires password
        git clone https://edwardsmith999@bitbucket.org/granlammps/gitlammps.git gitlammps
        cd ./gitlammps
        git checkout common
        cd ../;;
      *)
        echo "Unknown lammps:" $1; exit 1;;
    esac
fi

