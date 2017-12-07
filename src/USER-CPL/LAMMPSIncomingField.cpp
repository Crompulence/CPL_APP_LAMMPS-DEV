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


