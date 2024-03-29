import numpy as np
from scipy.signal import argrelextrema

def analytical_gravity(z0, v0, t0, t, g=9.81):
    return z0 + v0*(t-t0) - 0.5*g*(t-t0)**2


def read_data(logfile='./log.lammps', 
              datafile='./thermo_output.txt'):

    #Get data from file
    import os
    print((os.path.dirname(os.path.realpath(__file__)), os.getcwd()))
    with open(logfile) as f:
        for l in f.readlines():
            if l.find("timestep") != -1:
                dt = float(l.strip('timestep'))
                break
        
    #print("timestep = ", dt)
    data = np.genfromtxt(datafile)
    t = data[:,0]*dt
    z = data[:,1]
    v = data[:,2]
    f = data[:,3]

    return t, z, v, f


def check_bouncing_error_vs_gravity(D=3.5e-4, g=9.81, 
                                    logfile='./log.lammps', 
                                    datafile='./thermo_output.txt',
                                    plot=False):

    t, z, v, f = read_data(logfile, datafile)

    #Get section between bounces
    zp = np.copy(z)
    zp[z<D/2.]=1000.
    mins = argrelextrema(zp, np.less)[0]
    m = mins[0]
    mp1 = mins[1]
    ta = t[m:mp1]

    #Error vs gravity
    error = []
    for i in range(1,len(mins)-1,2):
        m = mins[i]
        mp1 = mins[i+1]
        ta = t[m:mp1]
        za = analytical_gravity(z[m], v[m], t[m], ta, g=g)
        error.append([ta, (z[m:mp1]-za)/za])
        if plot:
            print((i, mins))
            plt.plot(ta, za, 'ro')
            plt.plot(t, z, 'k-')
            plt.show()

    return np.array(error)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    t, z, v, f = read_data()
    zp = np.copy(z)
    mins = argrelextrema(zp, np.less)[0]
    plt.plot(zp)
    plt.show()

    
    error = check_bouncing_error_vs_gravity(plot=True)
    plt.plot(error[0,0,:],error[0,1,:])
    plt.show()
