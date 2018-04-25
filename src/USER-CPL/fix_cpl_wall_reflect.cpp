#include "fix_cpl_wall_reflect.h"
#include <iostream>
#include <string.h>
#include "update.h"
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
  Fb = std::valarray<double>(nbins);
  // std::cout << narg << ": " << alpha << " " << low << " " << hi << " " << nbins << std::endl;
  fd.open ("fb.debug");
  // std::cout << "open: " << fd.good() << std::endl;
  fd << "# Chunk-averaged data for fix wreflect_forces and group norm" << std::endl;
  fd << "# Timestep Number-of-chunks" << std::endl;
  fd << "# Chunk Coord1 fb" << std::endl;
  begin_tstep = update->ntimestep;
  // std::cout << "TSTEP:" << begin_tstep << std::endl;
}



void FixCPLWallReflect::post_constructor() {
    // Setup dependencies
    (new LAMMPSDepRegion("wreflect_region", DepListT({}), this, lmp, 
                     &wr_region_depfunc))->addToPool(&depPool);

    (new LAMMPSDepCompute("wreflect_chunks", DepListT({"wreflect_region"}), this,
                     lmp, &wr_chunks_depfunc))->addToPool(&depPool);
    (new LAMMPSDepFix("wreflect_forces", DepListT({"wreflect_chunks"}), this,
                     lmp, &wr_forces_depfunc))->addToPool(&depPool);
    depPool["wreflect_forces"]->load();
    chunk_compute = static_cast<LAMMPS_NS::ComputeChunkAtom*>(\
                    static_cast<LAMMPSDepCompute*>(\
                    depPool["wreflect_chunks"])->compute);
    forces_fix = static_cast<LAMMPSDepFix*>(\
                 depPool["wreflect_forces"])->fix;
 

 
}

void FixCPLWallReflect::setup(int vflag) {
    post_force(vflag);
    end_of_step();
}

int FixCPLWallReflect::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  return mask;
}

void FixCPLWallReflect::end_of_step() {
    if ((begin_tstep - update->ntimestep) % 1000 == 0) {
        fd << update->ntimestep << " " << nbins << std::endl;
        for (int i; i < nbins; i++)
            fd << "  " << i+1 << " " << dy*i+dy/2.0+low << " " << Fb[i] << std::endl;
    }
}


void FixCPLWallReflect::post_force(int vflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  // Needed to get the appropriate chunking in the current step
  // ave/chunk updates in end_of_step()
  int nchunk = chunk_compute->setup_chunks();
  chunk_compute->compute_ichunk();
  int *ichunk = chunk_compute->ichunk;
  double **f = atom->f;
  double **x = atom->x;
  // std::cout << "VALS: " << forces_fix->compute_array(0, 0) << " " << forces_fix->compute_array(0,1) << " " << forces_fix->compute_array(0,2) << std::endl;
  // std::cout << "Timestep:" << update->ntimestep << std::endl;
  // std::cout << "VAL1: " << forces_fix->compute_array(0, 2) << " " << forces_fix->compute_array(0, 1) << std::endl;
  forces_fix->end_of_step();
  // std::cout << "VAL2: " << forces_fix->compute_array(0, 2) << " " << forces_fix->compute_array(0, 1) << std::endl;

  // Update Force
  for (int i = 0; i < nchunk; i++)
    if (forces_fix->compute_array(i, 1) >= 1){
        Fb[i] = alpha * forces_fix->compute_array(i, 2) + (1 - alpha) * Fb[i];
    }

  double val = 0.0;
  int c = 0;

  // Force correction
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && ichunk[i] > 0)
        // if (ichunk[i] == 1) {
        //     val += f[i][1];
        //     c+=1;
        // }
        f[i][1] -= Fb[ichunk[i]-1];
  }
  // std::cout << "VALS3: " << val/c << " " << c <<std::endl << std::endl;
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
            // << "wreflect_chunks fy norm all file wreflect_forces.debug";
            << "wreflect_chunks fy norm all";
    std::cout << "CPL: " << str_out.str() << std::endl;
    return str_out.str();
}



