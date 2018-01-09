#include <iostream>
#include "fix_cpl_init.h"
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSOutgoingField.h"
#include "LAMMPSIncomingField.h"



FixCPLInit::FixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg), cnstFixDefined(false), bcFixDefined(false) {

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
    if (cnstFixDefined && bcFixDefined)
        cplsocket.allocateBuffers(bcPool, cnstPool);
}

void FixCPLInit::setup(int vflag) {
    if (cnstFixDefined && bcFixDefined)
        end_of_step();
}

void FixCPLInit::end_of_step() {
    if (cnstFixDefined && bcFixDefined)
        cplsocket.communicate(bcPool, cnstPool);
}

FixCPLInit::~FixCPLInit() {
	cplsocket.finalizeComms();
}
