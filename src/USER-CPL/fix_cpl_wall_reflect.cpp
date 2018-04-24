#include "fix_cpl_wall_reflect.h"
#include <iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>
#include <unistd.h>
#include "LAMMPSDep.h"
#include <iomanip>
#include "force.h"



FixCPLWallReflect::FixCPLWallReflect(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

  double *global_domain = lammps->domain->prd;
  double region_width;
  if (narg != 6)
      error->all(FLERR,"Illegal fix wall/reflect command");
  // std::cout << narg << ": " << arg[0] << "," << arg[1] << "," << arg[2] << std::endl;
  alpha = force->numeric(FLERR, arg[3]);
  region_width = force->numeric(FLERR, arg[4]);
  nbins = (int) force->numeric(FLERR, arg[5]);
  hi = global_domain[1];
  low = hi - region_width;
  dy = (hi - low) / nbins;
  // std::cout << narg << ": " << alpha << " " << low << " " << hi << " " << nbins << std::endl;
}

void FixCPLWallReflect::post_constructor() {
    // Setup dependencies
    (new LAMMPSDepRegion("wreflect_region", DepListT({}), this, lmp, 
                     &wr_region_depfunc))->addToPool(&depPool);

    (new LAMMPSDepCompute("wreflect_chunks", DepListT({"wreflect_region"}), this,
                     lmp, &wr_chunks_depfunc))->addToPool(&depPool);
    (new LAMMPSDepFix("wreflect_forces", DepListT({"wreflect_chunks"}), this,
                     lmp, &wr_forces_depfunc))->addToPool(&depPool);
    // depPool["wreflect_region"]->load();
    // depPool["wreflect_chunks"]->load();
    depPool["wreflect_forces"]->load();
 
}


DEPFUNC_IMP(wr_region_depfunc) {
    FixCPLWallReflect* wall_reflect = static_cast<FixCPLWallReflect*>(obj);
    std::stringstream str_out;
    str_out << std::setprecision(15) << "region " << dep_name << " block "\
            << "EDGE" << " " << "EDGE" << " "\
            << wall_reflect->low << " " << wall_reflect->hi << " "\
            << "EDGE" << " " << "EDGE" << " "\
            << "units box";
    std::cout << str_out.str() << std::endl;
    return str_out.str();
}

DEPFUNC_IMP(wr_chunks_depfunc) {
    FixCPLWallReflect* wall_reflect = static_cast<FixCPLWallReflect*>(obj);
    std::stringstream str_out;
    str_out << std::setprecision(15) << "compute "  << dep_name << " "\
            << "all " << "chunk/atom bin/1d "\
            << "y " << wall_reflect->low << " " << wall_reflect->dy << " "\
            << "ids every region " << "wreflect_region "\
            << "units box "\
            << "bound y " << wall_reflect->low << " " << wall_reflect->hi;
    std::cout << "CPL: " << str_out.str() << std::endl;
    return str_out.str();
}

DEPFUNC_IMP(wr_forces_depfunc) {
    std::stringstream str_out;
    str_out << "fix "  << dep_name << " "\
            << "all " << "ave/chunk "\
            << "1 1 1 "\
            << "wreflect_chunks fy norm all file wreflect_forces.debug";
    std::cout << "CPL: " << str_out.str() << std::endl;
    return str_out.str();
}
