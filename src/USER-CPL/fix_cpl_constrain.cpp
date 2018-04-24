#include "fix_cpl_init.h"
#include "fix_cpl_constrain.h"
#include <iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSIncomingField.h"



FixCPLConstrain::FixCPLConstrain(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

    int ifix = lammps->modify->find_fix("cplfix");
    if (ifix == -1)
        error->all(FLERR, "Fix cpl/constrain has been called before cpl/init.");
    fixCPLInit = static_cast<FixCPLInit*>(lmp->modify->fix[ifix]);
    fixCPLInit->cnstPool = CPL::IncomingFieldPool(cplsocket.cnstPortionRegion, cplsocket.cnstRegion);
    cnstPool = &(fixCPLInit->cnstPool);
    depPool = &(fixCPLInit->depPool);
    
    //Instantiate pools
    std::string force_type;
    CPL::get_file_param("constrain.momentum", "type", force_type);
    if (force_type== "Undefined")
        lammps->error->all(FLERR,"Missing forcetype option in cpl/init.");
    else {
        if (force_type !="Flekkoy")
            lammps->error->all(FLERR,"Missing or invalid argument in forcetype option in cpl/init.");
        else {
            (new StressIncomingField("stresscnst", DepListT({"cplforcefix"}), depPool, lammps))->addToPool(cnstPool);
        }
    }
}

void FixCPLConstrain::post_constructor() {
    (new LAMMPSDepRegion("cplforceregion", DepListT({}), this, lmp, 
                    &cplforceregion_depfunc))->addToPool(depPool);
    (new LAMMPSDepGroup("cplforcegroup", DepListT({"cplforceregion"}), this, lmp, 
                    &cplforcegroup_depfunc))->addToPool(depPool);
    (new LAMMPSDepFix("cplforcefix", DepListT({"cplforcegroup", "cplforceregion"}), this,
                 lmp, &cplforcefix_depfunc))->addToPool(depPool);
 
    fixCPLInit->cnstPool.setupAll();
    fixCPLInit->cnstFixDefined = true;
    // Call this to set fix_cpl_init always as the last fix to be call
    // in the run() loop
    fixCPLInit->setas_last_fix();
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
