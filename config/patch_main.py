#!/usr/bin/env python
import sys
from shutil import copyfile

lammpsdir = sys.argv[1]
mainfile = lammpsdir + "/main.cpp"

with open(mainfile, "r") as old:
    with open("./main_prepatch.cpp", "w+") as new:
        for l in old:
            if '#include "cpl.h' in l:
                sys.exit("File already patched")
            else:
                new.write(l)

with open(mainfile, "r") as old:
    with open("./main_patched.cpp", "w+") as new:
        for l in old:
            if "<mpi.h>" in l:
                new.write(l)
                new.write('#include "cpl.h"\n')
            elif "MPI_Init" in l:
                new.write(l)
                new.write('\n')
                new.write('  MPI_Comm comm;\n')
                new.write('  CPL::init(CPL::md_realm, comm);\n')
            elif "MPI_COMM_WORLD" in l:
                new.write(l.replace("MPI_COMM_WORLD","comm"))

            elif "MPI_Finalize" in l:
                new.write('CPL::finalize();\n')
                new.write(l)
            else:
                new.write(l)


#Replace main.cpp with patched file
copyfile("./main_patched.cpp", mainfile)
