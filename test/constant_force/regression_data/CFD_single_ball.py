import numpy as np
from mpi4py import MPI
from cplpy import CPL

g = 9.81
mi = -5.9490638385009208e-08

print("After import")

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)

print("After CPL init")

## Parameters of the cpu topology (cartesian grid)
npxyz = [1, 1, 1]
xyzL = [1.5E-003, 1.5E-003, 2.5E-003]
xyz_orig = [0.0, 0.0, 0.0]
ncxyz = [8, 8, 8]

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=3)

print("After CPL setup")


ft = True
for time in range(100):

    # send data to update
    send_array[2,:,:,:] = mi*g
    CPL.send(send_array)

    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array)

    print(("CFD time = ", time))

comm.Barrier()
CPL.finalize()
MPI.Finalize()




