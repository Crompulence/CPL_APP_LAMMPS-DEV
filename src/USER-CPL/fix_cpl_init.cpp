#include "fix_cpl_init.h"
#include<iostream>


fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {
   class LAMMPS_NS::LAMMPS *lmp=lammps;
   cplsocket.initMD(lammps);
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
	std::cout << "SEND: " <<  c << std::endl;
	c++;
    cplsocket.recvStress();
	cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();
}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}

