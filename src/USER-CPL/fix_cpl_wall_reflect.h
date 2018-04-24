#ifdef FIX_CLASS

FixStyle(cpl/wallreflect, FixCPLWallReflect)

#else

#ifndef LMP_FIX_CPL_WALL_REFLECT_H
#define LMP_FIX_CPL_WALL_REFLECT_H

#include "fix.h"
#include <memory>
#include "LAMMPSDep.h"
#include <valarray>

class FixCPLWallReflect: public LAMMPS_NS::Fix {

public:

    FixCPLWallReflect(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	virtual ~FixCPLWallReflect(){std::cout << "DESTRUCTROR WALLREFLECT"<< std::endl;};
    void post_constructor();
    int setmask();
    void setup(int vflag);
    void post_force(int vflag);

    // Public attributes
    double alpha;
    int nbins;
    double hi, low;
    double dy;
    DepPoolT depPool;

// protected:
    // std::valarray<double> Fb;

};

DEPFUNC_DEF(wr_region_depfunc);
DEPFUNC_DEF(wr_chunks_depfunc);
DEPFUNC_DEF(wr_forces_depfunc);

#endif
#endif
