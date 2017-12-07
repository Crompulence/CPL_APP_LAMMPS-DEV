#include "fix_cpl_init.h"
#include<iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSOutgoingField.h"
#include "LAMMPSIncomingField.h"



fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

    forceType = "Undefined";
    sendType = "Undefined";
    bndryAvg = "Undefined";

    //TODO: Improve line parsing
    for (int iarg=0; iarg<narg; iarg+=1){
        //std::cout << iarg << " " << arg[iarg] << std::endl;
        std::string arguments(arg[iarg]);
        if (arguments == "forcetype")
            if (iarg+1<narg)
                forceType = std::string(arg[iarg+1]);

        if (arguments == "sendtype")
            if (iarg+1<narg)
                sendType = std::string(arg[iarg+1]);

        if (arguments == "bndryavg")
            if (iarg+1<narg)
                bndryAvg = std::string(arg[iarg+1]);
    }

    if (sendType == "Undefined")
        lammps->error->all(FLERR,"Missing sendtype option in cpl/init.");
    else {
        if (sendType == "velocity") {
        }
        else if (sendType == "gran") {
        }
        else if (sendType == "granfull") {
        }
        else
            lammps->error->all(FLERR,"Missing or invalid argument in sendtype option in cpl/init.");
    }

    int bndry_avg_mode = -1;
    if (bndryAvg == "Undefined")
        lammps->error->all(FLERR,"Missing bndryavg option in cpl/init.");
    else {
        if (bndryAvg == "below")
            bndry_avg_mode = AVG_MODE_BELOW;
        else if (bndryAvg  == "above")
            bndry_avg_mode = AVG_MODE_ABOVE;
        else if (bndryAvg == "midplane")
            bndry_avg_mode = AVG_MODE_MIDPLANE;
        else
            lammps->error->all(FLERR,"Missing or invalid argument in bndryavg option in cpl/init.");
    }


   	lmp = lammps;
    cplsocket.setLammps(lmp);
   	cplsocket.init();
    cplsocket.configureBc(bndry_avg_mode);
    //TODO: Put here the configuration of sendtype if-else
    cplsocket.configureCnst(-1);
   	nevery = cplsocket.timestepRatio;

    //Instantiate pools
    bcPool = CPL::OutgoingFieldPool(cplsocket.bcPortionRegion, cplsocket.bcRegion);
    cnstPool = CPL::IncomingFieldPool(cplsocket.cnstPortionRegion, cplsocket.cnstRegion);
    if (sendType == "Undefined")
        lammps->error->all(FLERR,"Missing sendtype option in cpl/init.");
    else {
        if (sendType == "velocity") {
            (new VelOutgoingField("1velbc", cplsocket.bcPortionRegion,cplsocket.bcRegion, DepListT({"cfdbcfix"}), depPool, lammps))->addToPool(bcPool);
            (new NbinOutgoingField("2nbinbc", cplsocket.bcPortionRegion, cplsocket.bcRegion, DepListT({"cfdbcfix"}), depPool, lammps))->addToPool(bcPool);
        }
        else if (sendType == "gran") {
        }
        else if (sendType == "granfull") {
        }
        else
            lammps->error->all(FLERR,"Missing or invalid argument in sendtype option in cpl/init.");
    }
    if (forceType == "Undefined")
        lammps->error->all(FLERR,"Missing forcetype option in cpl/init.");
    else {
        if (forceType !="Flekkoy")
            lammps->error->all(FLERR,"Missing or invalid argument in forcetype option in cpl/init.");
        else {
            (new StressIncomingField("stresscnst", cplsocket.cnstPortionRegion, cplsocket.cnstRegion,DepListT({"cplforcefix"}), depPool, lammps))->addToPool(cnstPool);
        }
    }
}

void fixCPLInit::post_constructor() {
    setupDeps();
    bcPool.setupAll();
    cnstPool.setupAll();
    cplsocket.allocateBuffers(bcPool, cnstPool);
    setas_last_fix();
}


void fixCPLInit::setas_last_fix() {
   int ifix = lmp->modify->find_fix("cplfix");
   LAMMPS_NS::Fix* fix_aux = lmp->modify->fix[ifix];
   int nfix = lmp->modify->nfix; 
   lmp->modify->fix[ifix] = lmp->modify->fix[nfix-1];
   lmp->modify->fix[nfix-1] = fix_aux;
   lmp->modify->fmask[ifix] = lmp->modify->fix[ifix]->setmask();
   lmp->modify->fmask[nfix-1] = lmp->modify->fix[nfix-1]->setmask();
}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  return mask;
}

void fixCPLInit::init() {
}

void fixCPLInit::setup(int vflag) {
 	end_of_step();
}

void fixCPLInit::end_of_step() {
    cplsocket.communicate(bcPool, cnstPool);
}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}

void fixCPLInit::setupDeps() {
    (new LAMMPSDepRegion("cfdbcregion", DepListT({}), cplsocket, lmp, 
                    &cfdbcregion_depfunc))->addToPool(depPool);
    (new LAMMPSDepCompute("cfdbccompute", DepListT({"cfdbcregion"}), cplsocket,
                     lmp, &cfdbccompute_depfunc))->addToPool(depPool);
    (new LAMMPSDepFix("cfdbcfix", DepListT({"cfdbccompute"}), cplsocket,
                 lmp, &cfdbcfix_depfunc))->addToPool(depPool);
    (new LAMMPSDepRegion("cplforceregion", DepListT({}), cplsocket, lmp, 
                    &cplforceregion_depfunc))->addToPool(depPool);
    (new LAMMPSDepGroup("cplforcegroup", DepListT({"cplforceregion"}), cplsocket, lmp, 
                    &cplforcegroup_depfunc))->addToPool(depPool);
    (new LAMMPSDepFix("cplforcefix", DepListT({"cplforcegroup", "cplforceregion"}), cplsocket,
                 lmp, &cplforcefix_depfunc))->addToPool(depPool);
 
}


DEPFUNC_IMP(cfdbcregion_depfunc) {
    std::valarray<double> bounds = cplsocket.bcRegion.bounds;
    std::stringstream str_out;
    str_out << "region " << dep_name << " block "\
            << bounds[0] << " " << bounds[1] << " "\
            << bounds[2] << " " << bounds[3] << " "\
            << bounds[4] << " " << bounds[5] << " "\
            << "units box";
    return str_out.str();
}

DEPFUNC_IMP(cfdbccompute_depfunc) {
    std::valarray<double> bounds = cplsocket.bcRegion.bounds;
    std::stringstream str_out;
    str_out << "compute "  << "cfdbccompute "\
            << "all " << "chunk/atom bin/3d "\
            << "x lower " << cplsocket.dx << " "\
            << "y lower " << cplsocket.dy << " "\
            << "z lower " << cplsocket.dz << " "\
            << "ids every region " << "cfdbcregion "\
            << "units box bound y " << bounds[2] << " "\
            << bounds[3];
    return str_out.str();
}

DEPFUNC_IMP(cfdbcfix_depfunc) {
    std::stringstream str_out;
    str_out << "fix "  << "cfdbcfix "\
            << "all " << "ave/chunk "\
            << "1 1 " << cplsocket.timestepRatio << " "\
            << "cfdbccompute vx vy vz norm all";
    return str_out.str();
}

DEPFUNC_IMP(cplforceregion_depfunc) {
    std::valarray<double> bounds = cplsocket.cnstRegion.bounds;
    std::stringstream str_out;
    str_out << "region " << dep_name << " block "\
            << bounds[0] << " " << bounds[1] << " "\
            << bounds[2] << " " << bounds[3] << " "\
            << bounds[4] << " " << bounds[5] << " "\
            << "units box";
    return str_out.str();
}
DEPFUNC_IMP(cplforcegroup_depfunc) {
    return "group " + dep_name + " dynamic all region cplforceregion every 1";
}
DEPFUNC_IMP(cplforcefix_depfunc) {
    std::stringstream str_out;
    str_out << "fix "  << dep_name << " "\
            << "cplforcegroup " << "cpl/force "\
            << "region cplforceregion";
    return str_out.str();
}
