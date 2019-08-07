#include "LAMMPSOutgoingField.h"

void TemperatureOutgoingField::setup() {
    data_size = 1;
    CPL::get_file_param("bc", "dimension", bc_dimension);
    loadDeps(depList, *depPool);
}

//TODO: Refactor into a function to iterate ave/time vector output array
void TemperatureOutgoingField::pack_(const std::vector<int>& glob_cell,
                             const std::vector<int>& loc_cell,
                             const std::valarray<double>& coord) {

    LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbc_fix"])->fix;
    double T;
    if (bc_dimension == "3d" || bc_dimension == "3D") {
            int row = glob_cell[0]*region.nCells[1]\
                                  *region.nCells[2] + glob_cell[1]\
                                  *region.nCells[2] + glob_cell[2];
        T = cfdbcfix->compute_array(row, 7);  
    }
    // For 1D case
    else {
       T = cfdbcfix->compute_array(0, 7);  
    }
    buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]) = T;
}

void VelOutgoingField::update() {
    // Iterate over the BC region to average velocities
    //
    if (bc_dimension == "3d" || bc_dimension == "3D") {
        if (averageVels) {
            std::vector<int> bc_region = region.cellBounds;
            if (CPL::is_proc_inside(bc_region.data())) {
                LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbc_fix"])->fix;
                vAvg = 0.0;
                int row;
                int N = region.nCells[0] * region.nCells[1] * region.nCells[2];
                for (int i = bc_region[0]; i <= bc_region[1]; i++)
                for (int j = bc_region[2]; j <= bc_region[3]; j++)
                for (int k = bc_region[4]; k <= bc_region[5]; k++) {
                    row = i*region.nCells[1]\
                           *region.nCells[2] + j\
                           *region.nCells[2] + k; 
                    vAvg[0] += cfdbcfix->compute_array(row, 4);
                    vAvg[1] += cfdbcfix->compute_array(row, 5);
                    vAvg[2] += cfdbcfix->compute_array(row, 6);
                }
                vAvg = vAvg / double(N);
                
            }
        }
    }
}

void VelOutgoingField::pack_(const std::vector<int>& glob_cell,
                             const std::vector<int>& loc_cell,
                             const std::valarray<double>& coord) {
    LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbc_fix"])->fix;
    double vx, vy, vz;
    if (bc_dimension == "3d" || bc_dimension == "3D") {
        int row = glob_cell[0]*region.nCells[1]\
                              *region.nCells[2] + glob_cell[1]\
                              *region.nCells[2] + glob_cell[2];
        if (averageVels) {
            vx = vAvg[0];
            vy = vAvg[1];
            vz = vAvg[2];
        }
        else {
            vx = cfdbcfix->compute_array(row, 4);  
            vy = cfdbcfix->compute_array(row, 5);  
            vz = cfdbcfix->compute_array(row, 6);  
        }
    }
    // For 1D case
    else {
        vx = cfdbcfix->compute_array(0, 4);  
        vy = cfdbcfix->compute_array(0, 5);  
        vz = cfdbcfix->compute_array(0, 6);  
    }
    buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]) = vx;
    buffer(1, loc_cell[0], loc_cell[1], loc_cell[2]) = vy;
    buffer(2, loc_cell[0], loc_cell[1], loc_cell[2]) = vz;
}

void VelOutgoingField::setup() {
    data_size = 3;
    CPL::get_file_param("bc.velocity", "average-vels", averageVels);
    CPL::get_file_param("bc", "dimension", bc_dimension);
    vAvg = std::valarray<double>(3);
    loadDeps(depList, *depPool);
}
void NbinOutgoingField::pack_(const std::vector<int>& glob_cell,
                              const std::vector<int>& loc_cell,
                              const std::valarray<double>& coord) {
    LAMMPS_NS::Fix* cfdbcfix = static_cast<LAMMPSDepFix*>((*depPool)["cfdbc_fix"])->fix;
	double ncount;
    if (bc_dimension == "3d" || bc_dimension == "3D") {
        int row = glob_cell[0]*region.nCells[1]\
                              *region.nCells[2] + glob_cell[1]\
                              *region.nCells[2] + glob_cell[2];
        ncount = cfdbcfix->compute_array(row, 3);  
    }
    else {
	    ncount = cfdbcfix->compute_array(0, 3);  
    }
    buffer(0, loc_cell[0], loc_cell[1], loc_cell[2]) = ncount;
}

void NbinOutgoingField::setup() {
    data_size = 1;
    CPL::get_file_param("bc", "dimension", bc_dimension);
    loadDeps(depList, *depPool);
}
