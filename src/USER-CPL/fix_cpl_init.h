#ifdef FIX_CLASS

FixStyle(cpl/init, FixCPLInit)

#else

#ifndef LMP_FIX_CPL_INIT_H
#define LMP_FIX_CPL_INIT_H

#include "fix.h"
#include "cpl/cpl.h"
#include "CPLSocketLAMMPS.h"
#include <memory>
#include "LAMMPSDep.h"

class FixCPLInit : public LAMMPS_NS::Fix {

public:

    FixCPLInit(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~FixCPLInit();
    int setmask();
    void setup (int vflag); 
    void end_of_step();
    void post_constructor();
    void setas_last_fix();
	void init();

    // Pool of boundary conditions for CFD
    CPL::OutgoingFieldPool bcPool;
    CPL::IncomingFieldPool cnstPool;
    
    std::string forceType;
    std::string sendType;
    std::string bndryAvg;
    DepPoolT depPool;
    bool cnstFixDefined, bcFixDefined;

private:
    LAMMPS_NS::LAMMPS* lmp;

};

#endif
#endif
