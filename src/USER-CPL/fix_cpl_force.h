#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include "fix.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
//#include "cpl/CPL_force.h"
#include <memory>

class FixCPLForce : public LAMMPS_NS::Fix {

public:

    FixCPLForce ( class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
    int setmask();
    void apply(); //todo add override <=== es205 17/01/17 WTF does this mean?
//	void post_force(int vflag);
	void setup (const CPL::DoubNdArray& stress, const CPL::IntVector& portion, double units);
    void updateProcPortion (const CPL::IntVector& portion);

private:

	const CPL::DoubNdArray* cfdStress;
    CPL::IntVector procPortion;
    double flekkoyGWeight (double y, double ymin, double ymax);
    double units_factor;

};

#endif
#endif
