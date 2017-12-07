#include "LAMMPSOutgoingField.h"

void VelOutgoingField::pack_(const std::vector<int>& glob_cell,
                             const std::vector<int>& loc_cell,
                             const std::valarray<double>& coord) {
    LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbcfix"])->fix;
	int row = glob_cell[0]*region.nCells[1]\
                          *region.nCells[2] + glob_cell[1]\
                          *region.nCells[2] + glob_cell[2];
    double vx = cfdbcfix->compute_array(row, 4);  
    double vy = cfdbcfix->compute_array(row, 5);  
    double vz = cfdbcfix->compute_array(row, 6);  
    buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]) = vx;
    buffer(1, loc_cell[0], loc_cell[1], loc_cell[2]) = vy;
    buffer(2, loc_cell[0], loc_cell[1], loc_cell[2]) = vz;
}

void VelOutgoingField::setup() {
    loadDeps(depList, *depPool);
    data_size = 3;
}
void NbinOutgoingField::pack_(const std::vector<int>& glob_cell,
                              const std::vector<int>& loc_cell,
                              const std::valarray<double>& coord) {
    LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbcfix"])->fix;
	int row = glob_cell[0]*region.nCells[1]\
                          *region.nCells[2] + glob_cell[1]\
                          *region.nCells[2] + glob_cell[2];
	double ncount = cfdbcfix->compute_array(row, 3);  
    buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]) = ncount;
}

void NbinOutgoingField::setup() {
    loadDeps(depList, *depPool);
    data_size = 1;
}
