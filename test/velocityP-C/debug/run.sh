mpiexec -n 27 python2 dummyCFD_vels.py : -n 27 lmp_cpl -in lammps_vels.in 
#mpiexec -n 27 lmp_cpl -in lammps_vels.in & PID=$!; mpiexec -n 27 python2 dummyCFD_vels.py ; wait $PID
