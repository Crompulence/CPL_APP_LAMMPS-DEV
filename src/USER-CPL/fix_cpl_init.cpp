#include "fix_cpl_init.h"
#include<iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>



fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

	if (narg < 8) error->all(FLERR,"Illegal fix cplinit command");
	int bndry_avg_mode;
	if (strcmp(arg[3],"below") == 0)
		bndry_avg_mode = AVG_MODE_BELOW;
	else if (strcmp(arg[3],"above") == 0)
		bndry_avg_mode = AVG_MODE_ABOVE;
	else if (strcmp(arg[3],"midplane") == 0)
		bndry_avg_mode = AVG_MODE_MIDPLANE;
	else
		error->all(FLERR,"Illegal fix cplinit command - averaging method should \
				   be 'below', 'above' or 'midplane' ");
    int units;
	if (strcmp(arg[5],"real") == 0)
        cplsocket.units = REAL_UNITS;
    else if (strcmp(arg[5],"lj") == 0)
        cplsocket.units = LJ_UNITS;
	else
		error->all(FLERR,"Illegal fix cplinit command - units should be 'real' or 'lj'");
   	class LAMMPS_NS::LAMMPS *lmp=lammps;
   	cplsocket.initMD(lammps);
	cplsocket.setBndryAvgMode(bndry_avg_mode);
   	nevery = cplsocket.timestep_ratio;
   	std::cout << "LAMMPS NEVERY: " << nevery <<std::endl;
}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init() {
   cplsocket.setupFixMDtoCFD(lmp);
   cplsocket.setupFixCFDtoMD(lmp);
}


void fixCPLInit::setup(int vflag) {
  	end_of_step();
}


void fixCPLInit::post_force(int vflag) {
	cplsocket.updateStress();
}

void fixCPLInit::end_of_step() {
	static int c = 0;
    // Communications
	//std::cout << "SEND: " <<  c << std::endl;
	c++;
    cplsocket.recvStress();
	cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();
}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}

