#include "local_regrid.h"

namespace sledgehamr {

LocalRegrid::LocalRegrid(Sledgehamr* owner) {
    sim = owner;
    ParseInput();
    CreateCommMatrix();
}

bool LocalRegrid::AttemptRegrid(const int lev) {
    bool res = DoAttemptRegrid(lev);
    // TODO Check this isn't causing a memory leak.
    layouts.clear();
    return res;
}

bool LocalRegrid::DoAttemptRegrid(const int lev) {
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

    // Check if we have new levels. At the (re)start all levels will be
    // considered new.
    while( numPts.size() < sim->finest_level ) {
        int l = numPts.size();
        numPts.push_back(sim->grid_new[l].boxArray().numPts());
        no_local_regrid.push_back(false);
    }

    // Skip regrid after a restart. TODO: Read force_global_regrid_at_restart.
    if (force_global_regrid_at_restart) {
        amrex::Print() << "Skipping local regrid after a restart." << std::endl;

        if (lev==sim->shadow_hierarchy)
            force_global_regrid_at_restart = false;
        return false;
    }

    // Now that we are sure we really want to attempt a local regrid initialize
    // data structures.
    layouts.resize(sim->finest_level - sim->shadow_hierarchy);
    for (int l=sim->shadow_hierarchy+1; l<=sim->finest_level; ++l) {
        int Np = sim->dimN[l] / sim->blocking_factor[l][0];
        for (int f=0; f<omp_get_max_threads(); ++f) {
            std::unique_ptr<UniqueLayout> ptr =
                    std::make_unique<UniqueLayout>(this, Np);
            layouts[l].push_back( std::move(ptr) );
        }
    }

    std::vector<double> min_distance(sim->finest_level, -1.);
    for (int l=lev; l<sim->finest_level && !no_local_regrid[l]; ++l) {
        min_distance[l+1] = DetermineNewBoxArray(l);
    }

    return true;
}

void LocalRegrid::ParseInput() {
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("force_global_regrid_at_restart",
                  force_global_regrid_at_restart);
    pp_amr.query("n_error_buf", n_error_buf);
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

int LocalRegrid::DetermineNewBoxArray(const int lev) {
    // TODO: Create custom function.
    double threshold = n_error_buf * n_error_buf;


    return 0;
}

}; // namespace sledgehamr
