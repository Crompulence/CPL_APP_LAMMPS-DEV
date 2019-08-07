mpiexec -n 1 python2 dummyCFD_bc.py : -n 3 lmp_cpl -in lammps_bc.in 
#mpiexec -n 27 lmp_cpl -in lammps_bc.in & PID=$!; mpiexec -n 27 python2 dummyCFD_bc.py ; wait $PID
