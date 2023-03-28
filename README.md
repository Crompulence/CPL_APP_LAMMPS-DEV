# CPL APP for LAMMPS [![CI](https://github.com/Crompulence/CPL_APP_LAMMPS-DEV/actions/workflows/main.yml/badge.svg)](https://github.com/Crompulence/CPL_APP_LAMMPS-DEV/actions/workflows/main.yml) [![Build Status](https://img.shields.io/docker/cloud/build/cpllibrary/cpl-lammps)](https://hub.docker.com/r/cpllibrary/cpl-lammps)

1 ) Pre-requisites for compilation
=================================

- A compiled version of CPL library (https://github.com/Crompulence/cpl-library) which has been added to the user's path using "source SOURCEME.sh"
- A copy of LAMMPS
- The same GCC and MPICH version used to build CPL library

2 ) Install
===========

This packages contains USER-CPL, which contains all the code to link CPL library to LAMMPS. This is a LAMMPS package which must be copied into the main src code of LAMMPS and "installed" which normally proceeds with calls "make yes-user-cpl". However, for simplicity, this is managed as follows:

First, change directory to CPL_APP_LAMMPS-DEV,

    cd CPL_APP_LAMMPS-DEV

 and create a file called CODE_INST_DIR which will tell the APP which LAMMPS to install USER-CPL for, 

    echo "/path/to/lammps/directory/" > CODE_INST_DIR

Next, the packages of interest can be tweaked in the config file, 

    vi ./config/lammps_packages.in
    
    
Then to install these packages, run

    cd config
    ./enable-packages.sh make   

which should turn on packages from lammps_packages.in.
Once you are happy with this, and you can use MPI_port to link the two codes, simply type

    make

in the top level should copy USER-CPL to the lammps directory and build LAMMPS with you packages, copying the coupled excutable back to the bin folder in CPL_APP_LAMMPS-DEV.

Note if you cannot use MPI_port to link the two codes, then you will need to run using MPMD (i.e. mpiexec -n 1 lmp_cpl : -n 1 cfd) with both codes sharing an MPI_COMM_WORLD, then you need to patch LAMMPS to remove any reference to MPI_COMM_WORLD. This can be done using'make patch-lammps' before calling make, as follows:

    make patch-lammps
    make

Note that this requires that the latest version of LAMMPS is manually checked to see how the patch should work, a constantly evolving process which may mean the patch will not be applied successfully. The key changes can be applied manually, by going into /path/to/lammps/directory/src and adjusting the main.cpp file. The aim is to  simply replace MPI_COMM_WORLD with the communicator returned by CPL::init(CPL::md_realm, comm). An example of this patch is included at the bottom of this README as an example of the key changes.

For a full discussion of this APP and CPL library itself, plese see the website www.cpl-library.org as well as the associated [wiki](http://www.cpl-library.org/wiki/index.php/Main_Page).




    diff --git a/src/main.cpp b/src/main.cpp
    index 82dac5a..23590ca 100644
    --- a/src/main.cpp
    +++ b/src/main.cpp
    @@ -17,6 +17,7 @@
     #include "error.h"
     #include <stdio.h>
     #include <stdlib.h>
    +#include "cpl.h"
     
     #if defined(LAMMPS_TRAP_FPE) && defined(_GNU_SOURCE)
     #include <fenv.h>
    @@ -48,9 +49,12 @@ int main(int argc, char **argv)
       feenableexcept(FE_OVERFLOW);
     #endif
     
    +MPI_Comm comm;
    +CPL::init(CPL::md_realm, comm);
    +
     #ifdef LAMMPS_EXCEPTIONS
       try {
    -    LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
    +    LAMMPS *lammps = new LAMMPS(argc,argv,comm);
         lammps->input->file();
         delete lammps;
       } catch(LAMMPSAbortException & ae) {
    @@ -60,7 +64,7 @@ int main(int argc, char **argv)
         exit(1);
       }
     #else
    -  LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
    +  LAMMPS *lammps = new LAMMPS(argc,argv,comm);
       lammps->input->file();
       delete lammps;
     #endif


3 ) License
==========

CPL_APP_LAMMPS is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.


4 ) Directory Structure
=========================

This application repository is structured as follows:

 - src - source files which include the USER-CPL package to copy into LAMMPS
 - test - a range of test cases run automatically on Travis CI
 - config - scripts to patch LAMMPS for MPMD mode, add the USER-CPL package and build coupled LAMMPS.

