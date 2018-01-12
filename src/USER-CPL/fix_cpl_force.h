#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include "fix.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
//#include "cpl/CPL_force.h"
#include "cpl/TransmittingField.h"
#include <memory>

class FixCPLForce : public LAMMPS_NS::Fix {

public:

    FixCPLForce (class LAMMPS_NS::LAMMPS *lammps,
                 int narg, char **arg);
    int setmask();
    void setup(int vflag);
    //TODO: Make this const
    void setup(CPL::TransmittingField& field) {cplField=&field;};
    CPL::TransmittingField* cplField;
    void post_force(int vflag);
    bool conversionDisabled;

private:
    double flekkoyGWeight (double y, double ymin, double ymax);

};

#endif
#endif
