
import numpy as np

from CouetteAnalytical import CouetteAnalytical as CA

def find_nearest_indx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest(array, value):
    idx = find_nearest_indx(array, value)
    return array[idx]

def read_rec(body, no):
    step, Ny, r = [int(a) for a in body[0].split()]
    u = []
    for i in range(1+no*(Ny+1), Ny+1+no*(Ny+1)):
        u.append([float(a) for a in body[i].split()])
    u = np.array(u)
    return u

def check_LAMMPS_vs_Analytical(fdir, uwall=1., plotevolve=False):

    with open(fdir + "/profile.wall.2d", "r") as f:
        filestr = f.read()

    header = filestr.split("\n")[:3]
    body = filestr.split("\n")[3:]

    u = read_rec(body, 0)
    y = u[:,1]

    nsteps = 25
    dt = 0.005
    nu = 1.7
    Ly = 40.3103
    liquidstart = 0.13*Ly
    liquidend = 0.77*Ly
    liquidregion = liquidend-liquidstart
    liquidbins = y.shape[0]
    Re = liquidregion/nu 
    CAObj = CA(Re=Re, U=uwall, Lmin=liquidstart, Lmax=liquidend, 
               npoints=10*liquidbins, nmodes=100*liquidbins)

    if "dynamic" in plotevolve:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        plt.ion()
        plt.show()
        recds = range(nsteps)
    elif "summary" in plotevolve:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        recds = [1, 3, 7, 24]
    else:
        raise IOError("plotevolve="+plotevolve+ 
                      " which is unknown, should be dynamic or summary")

    ft = True
    timeerror = []
    for rec in range(nsteps):

        try:
            u = read_rec(body, rec)
        except IndexError:
            print("Last value read, exiting loop")
            break

        if plotevolve:
            # 0=Chunk 1=Coord1 2=Ncount 3=density/mass
            # 4=vx 5=vy 6=vz 7=temperature
            if rec in recds:
                ax.plot(u[:,4], u[:,1]*Ly, 'o', label="vx")
 
        y_anal, u_anal = CAObj.get_vprofile(1000*(rec+0.5)*dt)

        indxs = []
        mn = 1
        mx = -2
        for i in range(mn, u.shape[0]+mx):
            indxs.append(find_nearest_indx(y_anal, u[i,1]*Ly))

        if plotevolve:
            if rec in recds:
                ax.plot(u_anal[indxs], y_anal[indxs], 'k-x', label="vx Analytical")

        error = (u_anal[indxs] - u[mn:mx,4])/uwall

        if "dynamic" in plotevolve:
            ax.plot(error, y_anal[indxs], 'r-', label="Error")
        
        assert np.max(np.abs(error[:-1])) < 0.2
        timeerror.append(error)
        #print("Error = ",  np.sum(error))

        if "dynamic" in plotevolve:

            if ft:
                plt.legend()
                ft = False
            #plt.xlim([-0.1, 1.2])
            plt.pause(.05)
            plt.cla()
        elif "summary" in plotevolve:
            pass


    timeerror = np.array(timeerror)

    if "dynamic" in plotevolve:
        plt.ioff()
        plt.pcolormesh(timeerror[:,:-1].T, cmap=plt.cm.RdYlBu_r)
        plt.colorbar()
        plt.show()
    else:
        ax.set_xlabel("$u$")
        ax.set_ylabel("$y$")
        plt.savefig("LAMMPS_Validation.png", bbox_inches="tight")

    #Check average error less than 4%
    #print(np.mean(np.abs(timeerror),0))
    assert np.all(np.mean(np.abs(timeerror[:,:-1]),0) < 0.04)


if __name__ == "__main__":
    for uwall in [0.6,0.7,0.8,0.9,1.0]:
        fdir = "./run" + str(uwall)
        check_LAMMPS_vs_Analytical(fdir, uwall=uwall, plotevolve="dynamic")

