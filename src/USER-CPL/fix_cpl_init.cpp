#include <iostream>
#include "fix_cpl_init.h"
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSOutgoingField.h"
#include "LAMMPSIncomingField.h"



// FixCPLInit::~FixCPLInit() { 
//
// }

FixCPLInit::FixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg), cnstFixDefined(false), bcFixDefined(false) {

    cplsocket.loadParamFile();
    CPL::get_file_param("constrain-enabled", "", cplsocket.recvEnabled);
    CPL::get_file_param("bc-enabled", "", cplsocket.sendEnabled);

   	lmp = lammps;
    cplsocket.setLammps(lmp);
   	cplsocket.init();
   	nevery = cplsocket.timestepRatio;
}

void FixCPLInit::post_constructor() {
    setas_last_fix();
}


void FixCPLInit::setas_last_fix() {
   int ifix = lmp->modify->find_fix("cplfix");
   LAMMPS_NS::Fix* fix_aux = lmp->modify->fix[ifix];
   int nfix = lmp->modify->nfix; 
   lmp->modify->fix[ifix] = lmp->modify->fix[nfix-1];
   lmp->modify->fix[nfix-1] = fix_aux;
   lmp->modify->fmask[ifix] = lmp->modify->fix[ifix]->setmask();
   lmp->modify->fmask[nfix-1] = lmp->modify->fix[nfix-1]->setmask();
}

int FixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  return mask;
}

void FixCPLInit::init() {
    if (bcFixDefined && !cplsocket.sendBuffAllocated)
        cplsocket.allocateSendBuffer(bcPool);
    if (cnstFixDefined && !cplsocket.recvBuffAllocated)
        cplsocket.allocateRecvBuffer(cnstPool);
}

void FixCPLInit::setup(int vflag) {
    end_of_step();
}

void FixCPLInit::end_of_step() {
    // By default cpl/bc and cpl/constrain had to be both
    // defined to start pack/unpack and send/recv.
    if (cnstFixDefined && bcFixDefined ||\
        !cplsocket.sendEnabled ||\
        !cplsocket.recvEnabled) {
        cplsocket.communicate(bcPool, cnstPool);
    }
}

FixCPLInit::~FixCPLInit() {
	cplsocket.finalizeComms();
}
