#include "local_regrid.h"

namespace sledgehamr {

LocalRegrid::LocalRegrid(Sledgehamr* owner) {
    sim = owner;
    CreateCommMatrix();
}

bool LocalRegrid::AttemptRegrid(const int lev) {
    // Veto if volume threshold such that it disables local regrid.
    if (volume_threshold_strong <= 1. ) {
        veto_level = lev - 1;
        return false;
    }

    // Reset veto level.
    veto_level = -1;

    // New level doesn't exist yet so we need to force global regrid.
    if (lev==sim->finest_level) {
        amrex::Print() << "Skip local regrid as the level to be regridded does "
                       << "not yet exist." << std::endl;
        return false;
    }

    // Skip regrid after a restart.  TODO: Check logic here !
    if (numPts[lev+1] == 0 && (lev==1 || force_global_regrid_at_restart)) {
        amrex::Print() << "Skipping local regrid after a restart." << std::endl;
        return false;
    }

    // Now that we are sure we really want to attempt a local regrid initialize
    // data structures.
    for (int l=sim->shadow_hierarchy+1; l<=sim->finest_level; ++l) {
    }

    return true;
}

void LocalRegrid::CreateCommMatrix() {
    int N = amrex::ParallelDescriptor::NProcs();
    comm_matrix.resize(N, std::vector<int>(N, 0));

    int max = std::log2(N);
    for (int i = 0; i < max; ++i) {
        int s = pow(2, i);
        for(int j = 0; j < s; ++j) {
            for(int k = 0; k < s; ++k) {
                comm_matrix[j+s][k  ] = comm_matrix[j][k] + s;
                comm_matrix[j  ][k+s] = comm_matrix[j][k] + s;
                comm_matrix[j+s][k+s] = comm_matrix[j][k];
            }
        }
    }
}

}; // namespace sledgehamr
