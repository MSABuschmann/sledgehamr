#include "boost/multi_array.hpp"

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
        WrapIndices(l);
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
            layouts[l].emplace_back( std::make_unique<UniqueLayout>(this, Np) );
        }
    }

    // Get the new required box array for each level. Might still violate
    // nesting.
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

void LocalRegrid::WrapIndices(const int lev) {
    const int dimN = sim->dimN[lev];
    const int bf = sim->blocking_factor[lev][0];
    const int N = dimN/bf;

    std::vector<int> indices(N+2);
    indices[0] = dimN - bf/2;
    indices[N+1] = bf/2;
    for (int i=1; i<=N; ++i)
        indices[i] = (static_cast<double>(i)-0.5) * static_cast<double>(bf);

    wrapped_index.push_back( indices );
}

int LocalRegrid::DetermineNewBoxArray(const int lev) {
    // TODO: Create custom function.
    const double threshold = n_error_buf * n_error_buf;

    const double dimNf = sim->dimN[lev+1];
    const int ibff = sim->blocking_factor[lev+1][0];
    const double bff = static_cast<double>(ibff);

    const amrex::BoxArray& ba = sim->grid_new[lev+1].boxArray();
    const LevelData& state = sim->grid_new[lev];

    // Tag cells.
    amrex::TagBoxArray tags(state.boxArray(), state.DistributionMap());
    sim->ErrorEst(lev, tags, state.t, 0);

    // Manual OpenMP reduction since it's all correlated (ugh).
    std::vector<int> min_distance2(omp_get_max_threads(), INT_MAX);
    std::vector<int> min_i(omp_get_max_threads(), -1);
    std::vector<int> min_j(omp_get_max_threads(), -1);
    std::vector<int> min_k(omp_get_max_threads(), -1);

#pragma omp parallel
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Array4<double const>& state_fab = state.array(mfi);
        const amrex::Array4<char>& tag_arr = tags.array(mfi);

        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Dim3 lo = amrex::lbound(tilebox);
        const amrex::Dim3 hi = amrex::ubound(tilebox);

        // Determine if tilebox is near a C/F border.
        const int cx0 = static_cast<int>(static_cast<double>(lo.x*2)/bff) + 1;
        const int cy0 = static_cast<int>(static_cast<double>(lo.y*2)/bff) + 1;
        const int cz0 = static_cast<int>(static_cast<double>(lo.z*2)/bff) + 1;
        const int cx1 = static_cast<int>(static_cast<double>(hi.x*2)/bff) + 1;
        const int cy1 = static_cast<int>(static_cast<double>(hi.y*2)/bff) + 1;
        const int cz1 = static_cast<int>(static_cast<double>(hi.z*2)/bff) + 1;

        const int xs = cx1 - cx0 + 1;
        const int ys = cy1 - cy0 + 1;
        const int zs = cz1 - cz0 + 1;

        const int xsb = xs + 2;
        const int ysb = ys + 2;
        const int zsb = zs + 2;

        int remaining = 0;
        boost::multi_array<bool, 3> border;
        border.resize(boost::extents[xsb][ysb][zsb]);

        for (int k=-1; k<=zs; ++k) {
            for (int j=-1; j<=ys; ++j) {
                for (int i=-1; i<=xs; ++i) {
                    border[i+1][j+1][k+1] = !ba.contains( amrex::IntVect(
                            wrapped_index[lev+1][cx0+i],
                            wrapped_index[lev+1][cy0+j],
                            wrapped_index[lev+1][cz0+k]));

                    remaining += border[i+1][j+1][k+1];
                }
            }
        }

        // Now find tags and determine distances.
        for (int k = lo.z; k <= hi.z && remaining > 0; ++k) {
            for (int j = lo.y; j <= hi.y && remaining > 0; ++j) {
                // Probably no SIMD here but hey, we can try.
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (remaining == 0)
                        break;

                    if (tag_arr(i,j,k) != amrex::TagBox::SET)
                        continue;

                    const int i0 = i*2;
                    const int j0 = j*2;
                    const int k0 = k*2;

                    const int cx = static_cast<int>(
                            static_cast<double>(i0)/bff) + 1;
                    const int cy = static_cast<int>(
                            static_cast<double>(j0)/bff) + 1;
                    const int cz = static_cast<int>(
                            static_cast<double>(k0)/bff) + 1;

                    // Check each possible border.
                    double smallest2 = 1e99;
                    for (int kk=-1; kk<=1; ++kk) {
                        for (int jj=-1; jj<=1; ++jj) {
                            for (int ii=-1; ii<=1; ++ii) {
                                // You want to know what exactly we are
                                // calculating here? So do I! I wrote this
                                // almost a year ago and now I can't remember.
                                // It works though ...
                                const int refx = cx - cx0 + ii + 1;
                                const int refy = cy - cy0 + jj + 1;
                                const int refz = cz - cz0 + kk + 1;

                                if (!border[refx][refy][refz])
                                    continue;

                                amrex::IntVect smt(ibff*(cx+ii-1),
                                                   ibff*(cy+jj-1),
                                                   ibff*(cz+kk-1));
                                amrex::IntVect bgt(ibff*(cx+ii)-1,
                                                   ibff*(cy+jj)-1,
                                                   ibff*(cz+kk)-1);

                                int dx2 = smt[0]<i0 && i0<bgt[0] ? 0 :
                                        std::min((i0-smt[0])*(i0-smt[0]),
                                                 (i0-bgt[0])*(i0-bgt[0]));
                                int dy2 = smt[1]<j0 && j0<bgt[1] ? 0 :
                                        std::min((j0-smt[1])*(j0-smt[1]),
                                                 (j0-bgt[1])*(j0-bgt[1]));
                                int dz2 = smt[2]<k0 && k0<bgt[2] ? 0 :
                                        std::min((k0-smt[2])*(k0-smt[2]),
                                                 (k0-bgt[2])*(k0-bgt[2]));

                                int d2 =  dx2 + dy2 + dx2;
                                if (d2<min_distance2[omp_get_thread_num()]) {
                                    min_distance2[omp_get_thread_num()] = d2;
                                    min_i[omp_get_thread_num()] = i;
                                    min_j[omp_get_thread_num()] = j;
                                    min_k[omp_get_thread_num()] = k;
                                }

                                // Add new box if below threshold.
                                if (d2<threshold) {
                                    border[refx][refy][refz] = false;
                                    remaining--;
                                    layouts[lev+1][omp_get_thread_num()]->Add(
                                        static_cast<int>(static_cast<double>(
                                        wrapped_index[lev+1][cx+ii])/bff-0.5),
                                        static_cast<int>(static_cast<double>(
                                        wrapped_index[lev+1][cy+jj])/bff-0.5),
                                        static_cast<int>(static_cast<double>(
                                        wrapped_index[lev+1][cz+kk])/bff-0.5));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Combine box layouts.
    layouts[lev+1][0]->Merge(layouts[lev+1]);
    layouts[lev+1][0]->Distribute();

    // Collect global stats.
    int global_min_distance2 = INT_MAX;
    int global_min_i = -1, global_min_j = -1, global_min_k = -1;
    for (int f=0; f<omp_get_max_threads(); f++) {
        if (global_min_distance2 > min_distance2[f]) {
            global_min_distance2 = min_distance2[f];
            global_min_i = min_i[f];
            global_min_j = min_j[f];
            global_min_k = min_k[f];
        }
    }

    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    std::vector<int> all_min_distance2(nprocs);
    std::vector<int> all_min_i(nprocs);
    std::vector<int> all_min_j(nprocs);
    std::vector<int> all_min_k(nprocs);
    amrex::ParallelDescriptor::Gather(&global_min_distance2, 1,
                                      &all_min_distance2[0], 1, ioproc);
    amrex::ParallelDescriptor::Gather(&global_min_i, 1, &all_min_i[0], 1,
                                      ioproc);
    amrex::ParallelDescriptor::Gather(&global_min_j, 1, &all_min_j[0], 1,
                                      ioproc);
    amrex::ParallelDescriptor::Gather(&global_min_k, 1, &all_min_k[0], 1,
                                      ioproc);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (int n=0; n<nprocs; ++n) {
            if (global_min_distance2 > all_min_distance2[n]) {
                global_min_distance2 = all_min_distance2[n];
                global_min_i = all_min_i[n];
                global_min_j = all_min_j[n];
                global_min_k = all_min_k[n];
            }
        }

        amrex::Print() << "Shortest distance to C/F boundary: "
                       << sqrt(global_min_distance2) << " grid sites @ ("
                       << global_min_i << "," << global_min_j << ","
                       << global_min_k << ")" << std::endl;
    }

    return 0;
}

}; // namespace sledgehamr
