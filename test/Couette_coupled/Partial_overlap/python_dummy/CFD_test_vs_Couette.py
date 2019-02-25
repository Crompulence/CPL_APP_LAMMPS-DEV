############################################################################
#
#    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
#     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
#      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
#       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
#        _\/\\\_____________\/\\\/////////____\/\\\_____________
#         _\//\\\____________\/\\\_____________\/\\\_____________
#          __\///\\\__________\/\\\_____________\/\\\_____________
#           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
#            _______\/////////__\///______________\///////////////__
#
#
#                         C P L  -  L I B R A R Y
#
#           Copyright (C) 2012-2018 Edward Smith & David Trevelyan
#
# License
#
#    This file is part of CPL-Library.
#
#    CPL-Library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CPL-Library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.
#
# Description
#
#    Minimal CFD MOCK Script for debugging
#
# Author(s)
#
#    Edward Smith
#


import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

from CouetteAnalytical import CouetteAnalytical as CA

from cplpy import CPL

plotstuff=False

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)

## Parameters of the cpu topology (cartesian grid)
npxyz = [1, 1, 1]
# FCC lattice, plus wall width plus buffer region
wallwidth = 1.
bufsize = 2.
xyzL = [16.795961913825074, 45.349097, 16.795961913825074]
xyz_orig = [0.0, 0.0, 0.0]
ncxyz = [8, 8, 8]

dy = xyzL[1]/float(ncxyz[1])

Uwall = 1.
dt = 0.005*50
nu = 1.7
liquidstart = xyz_orig[1]+2.*wallwidth+2.*bufsize
liquidend = xyzL[1]-2.*wallwidth-2.*bufsize
liquidregion = liquidend-liquidstart
liquidbins = int(liquidregion/dy)
Re = liquidregion/nu   #Note Reynolds in independent of velocity in analytical fn
y = np.linspace(xyz_orig[1],  xyzL[1], ncxyz[1])

CAObj = CA(Re=Re, U=Uwall, Lmin=liquidstart, Lmax=liquidend, npoints=liquidbins, nmodes=100*liquidbins)

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=3)

if plotstuff:
    fig, ax = plt.subplots(1,1)
    plt.ion()
    plt.show()

for time in range(500):

    print("CFD time = ", time)

    # send data to update
    send_array[...] = 0.0
    send_array[0,:,:,:] = Uwall
    CPL.send(send_array)
        
    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array)

    #y_anal, u_anal = CAObj.get_vprofile(10*(time+1)*10)

    UBC = np.mean(recv_array,(1,2,3))
    #ax.imshow(recv_array[3,:,0,:])
    #ax.plot(y, np.mean(recv_array[3,:,:,:], (0,2))/1e2, '-or')
    #ax.plot(y, np.mean(recv_array[0,:,:,:], (0,2)), '-ob')
    #ax.plot(y_anal, u_anal, '-k')
    #ax.set_ylim([-0.1,1.4])
    #plt.pause(0.002)
    #plt.cla()

comm.Barrier()
CPL.finalize()
MPI.Finalize()




