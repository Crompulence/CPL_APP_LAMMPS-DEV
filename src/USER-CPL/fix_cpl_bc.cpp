#include "fix_cpl_init.h"
#include "fix_cpl_bc.h"
#include <iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include "LAMMPSOutgoingField.h"
#include <iomanip>

CPL::PortionField shiftBc(const CPL::PortionField& region_in) {
    double shift = 0.0;
    std::string avg_mode;
    CPL::get_file_param("bc.velocity", "compute-mode", avg_mode);
	if (avg_mode == "midplane") {
		shift = -(cplsocket.bcRegion.ly)/2.0;
	}
	else if (avg_mode == "above") {
        // Do nothing
    }
	else if (avg_mode == "below") {
		shift = -cplsocket.bcRegion.ly;
	}
    // Create a new field with the corrected BC domain 
    std::valarray<double> domain_bounds = region_in.bounds;
    domain_bounds[2] += shift;
    domain_bounds[3] += shift;
    return CPL::PortionField(CPL::Domain(domain_bounds), 
                             region_in.nCells,
                             region_in.cellBounds);
}




FixCPLBc::FixCPLBc(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

    //Instantiate pools
    int ifix = lammps->modify->find_fix("cplfix");
    if (ifix == -1)
        error->all(FLERR, "Fix cpl/bc has been called before cpl/init.");
    fixCPLInit = static_cast<FixCPLInit*>(lmp->modify->fix[ifix]);
    fixCPLInit->bcPool = CPL::OutgoingFieldPool(cplsocket.bcPortionRegion, cplsocket.bcRegion);
    bcPool = &(fixCPLInit->bcPool);
    depPool = &(fixCPLInit->depPool);
    std::string send_type;
    CPL::get_file_param("bc", "sendtype", send_type);
    if (send_type == "Undefined")
        lammps->error->all(FLERR,"Missing send_type option in cpl/init.");
    else {
        if (send_type == "velocity") {
            (new VelOutgoingField("1velocity", DepListT({"cfdbcfix"}), depPool, lammps))->addToPool(bcPool);
            (new NbinOutgoingField("2nbinbc", DepListT({"cfdbcfix"}), depPool, lammps))->addToPool(bcPool);
        }
        else if (send_type == "gran") {
        }
        else if (send_type == "granfull") {
        }
        else
            lammps->error->all(FLERR,"Missing or invalid argument in send_type option in cpl/init.");
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
    std::valarray<double> bounds = shiftBc(cplsocket.bcRegion).bounds;
    std::stringstream str_out;
    str_out << "region " << dep_name << " block "\
            << "EDGE" << " " << "EDGE" << " "\
            << bounds[2] << " " << bounds[3] << " "\
            << "EDGE" << " " << "EDGE" << " "\
            << "units box";
    return str_out.str();
}

//NOTE: Seems that 'upper' is not working for bound parameter. Use 'bounds' to get the appropriate
//      number of bins. The upper bound  is  ~1e-7 units lower due to dx*ncells is not exactly Lx.
DEPFUNC_IMP(cfdbccompute_depfunc) {
    std::valarray<double> bounds = shiftBc(cplsocket.bcRegion).bounds;
    std::stringstream str_out;
    str_out << std::setprecision(15) << "compute "  << "cfdbccompute "\
            << "all " << "chunk/atom bin/3d "\
            << "x lower " << cplsocket.dx << " "\
            << "y lower " << cplsocket.dy << " "\
            << "z lower " << cplsocket.dz << " "\
            << "ids every region " << "cfdbcregion "\
            << "units box " << " "\
            << "bound x " << bounds[0] << " " << bounds[1] << " "\
            << "bound y " << bounds[2] << " " << bounds[3] << " "\
            << "bound z " << bounds[4] << " " << bounds[5];
    std::cout << "CPL: " << str_out.str() << std::endl;
    return str_out.str();
}

DEPFUNC_IMP(cfdbcfix_depfunc) {
    std::stringstream str_out;
    str_out << "fix "  << "cfdbcfix "\
            << "all " << "ave/chunk "\
            << "1 1 " << cplsocket.timestepRatio << " "\
            << "cfdbccompute vx vy vz norm all file vels.debug";
    return str_out.str();
}
