#include "LAMMPSIncomingField.h"
#include "fix_cpl_force.h"

void StressIncomingField::setup() {
    data_size = 9;
    loadDeps(depList, *depPool);
    FixCPLForce* cplforcefix = static_cast<FixCPLForce*>(\
                               static_cast<LAMMPSDepFix*>(\
                               (*depPool)["cplforcefix"])->fix);
    cplforcefix->setup(*this);
}

// NOTE: NO unpack is necessary as post_force() is being called
// in fix_cpl_force at the each time-step.
