mpiexec -n 1 lmp_cpl -in lammps_forces.in : -n 1 python2 dummyCFD_forces.py
#mpiexec -n 27 lmp_cpl -in lammps_forces.in & PID=$!; mpiexec -n 27 python2 dummyCFD_forces.py ; wait $PID
