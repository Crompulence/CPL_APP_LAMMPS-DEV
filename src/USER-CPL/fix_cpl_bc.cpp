#include "fix_cpl_init.h"
#include "fix_cpl_bc.h"
#include <iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSOutgoingField.h"



FixCPLBc::FixCPLBc(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

    sendType = "Undefined";
    bndryAvg = "Undefined";

    //TODO: Improve line parsing
    for (int iarg=0; iarg<narg; iarg+=1){
        //std::cout << iarg << " " << arg[iarg] << std::endl;
        std::string arguments(arg[iarg]);

        if (arguments == "sendtype")
            if (iarg+1<narg)
                sendType = std::string(arg[iarg+1]);

        if (arguments == "bndryavg")
            if (iarg+1<narg)
                bndryAvg = std::string(arg[iarg+1]);
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


    cplsocket.configureBc(bndry_avg_mode);
    //Instantiate pools
    //TODO: Find fix_cpl_init
    int ifix = lammps->modify->find_fix("cplfix");
    fixCPLInit = static_cast<FixCPLInit*>(lmp->modify->fix[ifix]);
    fixCPLInit->bcPool = CPL::OutgoingFieldPool(cplsocket.bcPortionRegion, cplsocket.bcRegion);
    bcPool = &(fixCPLInit->bcPool);
    depPool = &(fixCPLInit->depPool);
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
}

void FixCPLBc::post_constructor() {
    // Setup dependencies
    (new LAMMPSDepRegion("cfdbcregion", DepListT({}), cplsocket, lmp, 
                     &cfdbcregion_depfunc))->addToPool(depPool);
    (new LAMMPSDepCompute("cfdbccompute", DepListT({"cfdbcregion"}), cplsocket,
                     lmp, &cfdbccompute_depfunc))->addToPool(depPool);
    (new LAMMPSDepFix("cfdbcfix", DepListT({"cfdbccompute"}), cplsocket,
                     lmp, &cfdbcfix_depfunc))->addToPool(depPool);
 
    fixCPLInit->bcPool.setupAll();
    fixCPLInit->bcFixDefined = true;
    // Call this to set fix_cpl_init always as the last fix to be call
    // in the run() loop
    fixCPLInit->setas_last_fix();
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
