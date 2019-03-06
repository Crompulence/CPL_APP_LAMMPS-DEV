mpiexec -n 1 lmp_cpl -in lammps_vels.in : -n 1 python2 dummyCFD_vels.py
#mpiexec -n 27 lmp_cpl -in lammps_vels.in & PID=$!; mpiexec -n 27 python2 dummyCFD_vels.py ; wait $PID
