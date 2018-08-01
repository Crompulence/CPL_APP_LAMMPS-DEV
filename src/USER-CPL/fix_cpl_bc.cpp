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


FixCPLBc:: ~FixCPLBc() {
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
            (new VelOutgoingField("1velocity", DepListT({"cfdbc_fix"}), depPool, lammps))->addToPool(bcPool);
            (new NbinOutgoingField("2nbinbc", DepListT({"cfdbc_fix"}), depPool, lammps))->addToPool(bcPool);
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
    (new LAMMPSDepRegion("cfdbcregion", DepListT({}), this, lmp, 
                     &cfdbcregion_depfunc))->addToPool(depPool);
    (new LAMMPSDepCompute("cfdbc_chunks", DepListT({"cfdbcregion"}), this,
                     lmp, &cfdbc_chunks_depfunc))->addToPool(depPool);
    (new LAMMPSDepCompute("cfdbc_property", DepListT({"cfdbc_chunks"}), this,
                     lmp, &cfdbc_property_depfunc))->addToPool(depPool);
    (new LAMMPSDepCompute("cfdbc_vcom", DepListT({"cfdbc_chunks"}), this,
                     lmp, &cfdbc_vcom_depfunc))->addToPool(depPool);
    (new LAMMPSDepFix("cfdbc_fix", DepListT({"cfdbc_chunks", "cfdbc_property", "cfdbc_vcom"}), this,
                     lmp, &cfdbc_fix_depfunc))->addToPool(depPool);

    // (new LAMMPSDepFix("cfdbc_fix", DepListT({"cfdbc_chunks"}), this,
    //                  lmp, &cfdbc_fix_depfunc))->addToPool(depPool);
 
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
DEPFUNC_IMP(cfdbc_chunks_depfunc) {
    std::valarray<double> bounds = shiftBc(cplsocket.bcRegion).bounds;
    std::stringstream str_out;
    str_out << std::setprecision(16) << "compute "  << "cfdbc_chunks "\
            << "all " << "chunk/atom bin/3d "\
            << "x lower " << cplsocket.dx << " "\
            << "y lower " << cplsocket.dy << " "\
            << "z lower " << cplsocket.dz << " "\
            << "ids every region " << "cfdbcregion "\
            << "units box " << " "\
            << "bound x " << bounds[0] << " " << bounds[1] << " "\
            << "bound y " << bounds[2] << " " << bounds[3] << " "\
            << "bound z " << bounds[4] << " " << bounds[5];
    return str_out.str();
}

DEPFUNC_IMP(cfdbc_property_depfunc) {
    std::stringstream str_out;
    str_out << "compute cfdbc_property all property/chunk "\
            << "cfdbc_chunks coord1 coord2 coord3 count";
    std::cout << "1 :" << str_out.str();
    return str_out.str();
}

DEPFUNC_IMP(cfdbc_vcom_depfunc) {
    std::stringstream str_out;
    str_out << "compute cfdbc_vcom all vcm/chunk cfdbc_chunks";
    std::cout << "2:" << str_out.str();
    return str_out.str();
}

DEPFUNC_IMP(cfdbc_fix_depfunc) {
    int samples;
    CPL::get_file_param("bc.velocity", "samples", samples);
    if (samples <= 0 || samples > cplsocket.timestepRatio)
        lmp->error->all(FLERR,"Number of samples has to be >= 0 and <= timestep ratio.");
    if (cplsocket.timestepRatio % samples != 0)
        lmp->error->all(FLERR,"Number of samples has to be a divisor of timestep ratio.");
    int sample_every = cplsocket.timestepRatio / samples;
    std::stringstream str_out;
    str_out << "fix "  << "cfdbc_fix "\
            << "all " << "ave/time "\
            << sample_every << " " << samples << " "\
            << cplsocket.timestepRatio << " "\
            << "c_cfdbc_property[*] " << "c_cfdbc_vcom[*][1] "\
            << "mode vector "\
            << "file velocity.debug";
    std::cout << "3:" << str_out.str();
    return str_out.str();
}

// DEPFUNC_IMP(cfdbc_fix_depfunc) {
//     std::stringstream str_out;
//     str_out << "fix "  << "cfdbc_fix "\
//             << "all " << "ave/chunk "\
//             << "1 " << cplsocket.timestepRatio << " "\
//             << cplsocket.timestepRatio << " "\
//             << "cfdbc_chunks vx vy vz norm all file vels.debug";
//     std::cout << "CPL: " << str_out.str() << std::endl;
//     return str_out.str();
// }
