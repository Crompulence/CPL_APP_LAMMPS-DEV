#ifdef FIX_CLASS

FixStyle(cpl/wallreflect, FixCPLWallReflect)

#else

#ifndef LMP_FIX_CPL_WALL_REFLECT_H
#define LMP_FIX_CPL_WALL_REFLECT_H

#include "fix.h"
#include <memory>
#include "LAMMPSDep.h"
#include <valarray>
#include "compute_chunk_atom.h"
#include <fstream>

class FixCPLWallReflect: public LAMMPS_NS::Fix {

public:

    FixCPLWallReflect(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	virtual ~FixCPLWallReflect(){fd.close();};
    void post_constructor();
    int setmask();
    void setup(int vflag);
    void post_force(int vflag);
    void end_of_step();

    // Public attributes
    double alpha;
    int nbins;
    double hi, low;
    double dy;
    DepPoolT depPool;
    LAMMPS_NS::ComputeChunkAtom* chunk_compute;
    LAMMPS_NS::Fix* forces_fix;
    std::ofstream fd;
    int begin_tstep;

protected:
    std::valarray<double> Fb;

};

DEPFUNC_DEF(wr_region_depfunc);
DEPFUNC_DEF(wr_chunks_depfunc);
DEPFUNC_DEF(wr_forces_depfunc);

#endif
#endif
