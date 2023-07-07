#include "boost/multi_array.hpp"

#include "local_regrid.h"

namespace sledgehamr {

LocalRegrid::LocalRegrid(Sledgehamr* owner) {
    sim = owner;
    ParseInput();
    CreateCommMatrix();

    no_local_regrid.resize(sim->max_level+1,false);
    do_global_regrid.resize(sim->max_level+1, false);
    nregrids = 0;
}

bool LocalRegrid::AttemptRegrid(const int lev) {
    amrex::Print() << std::endl;
    bool res = DoAttemptRegrid(lev);
    ClearLayout();
    return res;
}

void LocalRegrid::DidGlobalRegrid(const int lev) {
    for (int l = 0; l < sim->finest_level; ++l) {
        no_local_regrid[l] = false;
        do_global_regrid[l] = false;
        nregrids = 0;
    }
}

void LocalRegrid::InitializeLayout(const int max_lev) {
    layouts.resize(sim->finest_level + 1);
    for (int l = 1; l <= sim->finest_level; ++l) {
        int Np = sim->dimN[l] / sim->blocking_factor[l][0];
        for (int f = 0; f < omp_get_max_threads(); ++f) {
            layouts[l].emplace_back( std::make_unique<UniqueLayout>(this, Np) );
        }
    }
}

void LocalRegrid::ClearLayout() {
    layouts.clear();  
}

void LocalRegrid::JoinBoxArrays(const int lev, amrex::BoxArray& ba) {
    amrex::BoxList bl = layouts[lev][0]->BoxList(sim->blocking_factor[lev][0]);
    bl.simplify(true);
    amrex::Vector<amrex::Box> bv = std::move(bl.data());
    amrex::AllGatherBoxes(bv);
    bl = amrex::BoxList(std::move(bv));
    bl.simplify(false);
    ba = amrex::BoxArray(std::move(bl));
}

void LocalRegrid::AddToLayout(const int lev, const int thread, const int i,
                              const int j, const int k) {
    layouts[lev][thread]->Add(i, j, k);
}

void LocalRegrid::FinalizeLayout(const int lev) {
    layouts[lev][0]->Merge(layouts[lev]);
    layouts[lev][0]->Distribute();
}

bool LocalRegrid::DoAttemptRegrid(const int lev) {
    if (nregrids++ >= max_local_regrids) {
        if (max_local_regrids > 0) {
            amrex::Print() << "Maximum number of local regrids reached: "
                           << max_local_regrids << std::endl;
        }
        return false;
    }

    // Veto if volume threshold such that it disables local regrid.
    if (volume_threshold_strong <= 1. ) {
        amrex::Print() << "Local regrid disabled." << std::endl;
        veto_level = lev - 1;
        return false;
    }

    // Reset veto level.
    veto_level = -1;

    if (do_global_regrid[lev]) {
        amrex::Print() << "Skip local regrid in favour of global regrid."
                       << std::endl;
        return false;
    }

    // New level doesn't exist yet so we need to force global regrid.
    if (lev == sim->finest_level) {
        amrex::Print() << "Skip local regrid as the level to be regridded does "
                       << "not yet exist." << std::endl;
        return false;
    }

    // Skip regrid after a restart.
    if (force_global_regrid_at_restart) {
        amrex::Print() << "Skipping local regrid after a restart." << std::endl;

        if (lev == 0)
            force_global_regrid_at_restart = false;
        return false;
    }

    amrex::Print() << std::endl << "Attempting local regrid at level " << lev+1
                   << " and higher." << std::endl;

    // Check if we have new levels. At the (re)start all levels will be
    // considered new.
    while (last_numPts.size() <= sim->finest_level) {
        int l = last_numPts.size();
        last_numPts.push_back(sim->grid_new[l].boxArray().numPts());
        WrapIndices(l);
    }

    // Now that we are sure we really want to attempt a local regrid initialize
    // data structures.
    InitializeLayout(sim->finest_level);

    // Get the new required box array for each level. Might still violate
    // nesting.
    std::vector<double> min_distance(sim->finest_level + 1, -1.);
    for (int l = lev; l < sim->finest_level && !no_local_regrid[l]; ++l) {
        min_distance[l+1] = DetermineNewBoxArray(l);
    }

    // Make sure we do not violate nesting requirements.
    for (int l = sim->finest_level; l > 1; --l) {
        FixNesting(l);
    }

    // Join all boxes across MPI ranks.
    std::vector<amrex::BoxArray> box_arrays(sim->finest_level+1);
    for (int l = 1; l <= sim->finest_level; ++l) {
        JoinBoxArrays(l, box_arrays[l]);
    }

    // Check if we want to go ahead and add those boxes.
    bool veto = false;
    std::vector<double> latest_possible_regrid_time(sim->finest_level+1, -1.);
    for (int l = 1; l <= sim->finest_level; ++l) {
        double Nb = box_arrays[l].numPts();
        double Nc = sim->grid_new[l].boxArray().numPts();
        double Nr = last_numPts[l];

        if (Nr == 0)
            Nr = Nc;

        double dV = Nb / Nc;
        double fV = (Nb+Nc) / Nr;

        // This level exceeds strong threshold so we veto.
        if (fV > volume_threshold_strong)
            veto = true;

        // In case any level vetos the local regrid, this is the level on which
        // we want to perform the global regrid.
        if( fV > volume_threshold_weak && veto_level == -1 )
            veto_level = l - 1;

        if (l > lev) {
            double dx_c = n_error_buf;
            double regrid_dt = sim->time_stepper->regrid_dt[l];
            double dt_delay = DBL_MAX;
            if (min_distance[l] >= 0)
                dt_delay = min_distance[l] / dx_c * regrid_dt;

            latest_possible_regrid_time[l] = sim->grid_new[l].t + dt_delay;
        }

        amrex::Print() << "  Additional boxes on level " << l << " required: "
                       << box_arrays[l].size() << std::endl
                       << "    Instantanous volume increase: " << dV
                       << std::endl
                       << "    Volume increase since last global regrid: "
                       << fV << ". Threshold: " << volume_threshold_strong
                       << std::endl;

        if( l>lev ) {
            if (latest_possible_regrid_time[l] > sim->grid_new[l].t ) {
                if (min_distance[l] >= 0)
                    amrex::Print() << "    Could delay regridding this level "
                                   << "until: t = "
                                   << latest_possible_regrid_time[l]
                                   << std::endl;
            } else {
                amrex::Print() << "    Cannot delay this regrid."
                               << std::endl;
            }
        }
    }

    // Check what to do if we veto local regrid.
    if (veto) {
        amrex::Print() << "Local regrid has been vetoed. "
                       << "Global regrid on level " << veto_level
                       << " (adjusting level " << veto_level+1<<") "
                       << "deemed optimal." << std::endl;

        // We can do a globl regrid right away.
        if (veto_level >= lev)
            return false;

        // Request a global regrid. This flag will be checked by the time
        // stepper module.
        do_global_regrid[veto_level] = true;

        // Figure out when is the next possible time we can do a global regrid
        // at the requested level. For a sim with shadow hierarchy this is
        // either in one or two time steps from now to allow for computation of
        // truncation errors.
        double Nsteps = 2;
        if (sim->grid_new[veto_level].istep%2 == 0)
            Nsteps = 1;

        // For sim without shadow level we can do it once we have sync'ed with
        // that level.
        if (sim->shadow_hierarchy)
            Nsteps = 0;

        double regrid_target_time = sim->grid_new[veto_level].t
                                    + Nsteps * sim->dt[veto_level];

        // Check if we can satisfy target time.
        bool possible_to_wait = true;
        for (int l = lev+1; l <= sim->finest_level; ++l) {
            if (latest_possible_regrid_time[l] < regrid_target_time) {
                amrex::Print() << "Regrid cannot wait until "
                               << regrid_target_time
                               << " so will perform local regrid followed by "
                               << "global." << std::endl;
                possible_to_wait = false;
                break;
            }
        }

        if (possible_to_wait) {
            amrex::Print() << "Possible to delay regrid until "
                           << regrid_target_time << std::endl;

            for (int l = lev; l <= sim->finest_level; ++l)
                no_local_regrid[l] = true;

            return true;
        }
    }

    // Finally add new boxes to each grid.
    for (int l = 1; l <= sim->finest_level; ++l) {
        if (box_arrays[l].size() > 0) {
            AddBoxes(l, box_arrays[l]);
        }
    }

    return true;
}

void LocalRegrid::ParseInput() {
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("force_global_regrid_at_restart",
                  force_global_regrid_at_restart);
    pp_amr.query("n_error_buf", n_error_buf);
    pp_amr.query("max_local_regrids", max_local_regrids);
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
#pragma omp parallel for
    for (int i=1; i<=N; ++i)
        indices[i] = (static_cast<double>(i)-0.5) * static_cast<double>(bf);

    wrapped_index.push_back( indices );
}

double LocalRegrid::DetermineNewBoxArray(const int lev) {
    const double threshold = (n_error_buf+1) * (n_error_buf+1);

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

        if (remaining == 0)
            continue;

        // Now find tags and determine distances.
        for (int k = lo.z; k <= hi.z && remaining > 0; ++k) {
            for (int j = lo.y; j <= hi.y && remaining > 0; ++j) {
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
                    for (int kk=-1; kk<=1; ++kk) {
                        for (int jj=-1; jj<=1; ++jj) {
                            for (int ii=-1; ii<=1; ++ii) {
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

                                int d2 =  dx2 + dy2 + dz2;

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
    FinalizeLayout(lev+1);

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

        if( global_min_i >= 0 ) {
            amrex::Print() << "  Shortest distance to C/F boundary: "
                           << sqrt(static_cast<double>(global_min_distance2))
                           << " grid sites @ (" << global_min_i << ","
                           << global_min_j << "," << global_min_k << ")"
                           << std::endl;
        } else {
            amrex::Print() << "  Shortest distance to C/F boundary: "
                           << "> blocking_factor (" << ibff << " cells)"
                           << std::endl;
        }
    }

    if (global_min_i >= 0 )
        return sqrt(static_cast<double>(global_min_distance2));
    else
        return -1;
}

void LocalRegrid::FixNesting(const int lev) {
    amrex::BoxArray nest_ba = layouts[lev][0]->BoxArray(
                                    sim->blocking_factor[lev][0]);
    nest_ba.growcoarsen(sim->nghost+4, amrex::IntVect(2,2,2));
    nest_ba = WrapBoxArray(nest_ba, sim->dimN[lev-1]);
    amrex::BoxArray bak = sim->grid_new[lev-1].boxArray();
    const double bfc = static_cast<double>(sim->blocking_factor[lev-1][0]);

#pragma omp parallel for
    for (int b = 0; b < nest_ba.size(); ++b) {
        amrex::Box box = nest_ba.boxList().data()[b];
        const int cxl = static_cast<double>(box.smallEnd(0))/bfc + 1.;
        const int cyl = static_cast<double>(box.smallEnd(1))/bfc + 1.;
        const int czl = static_cast<double>(box.smallEnd(2))/bfc + 1.;
        const int cxh = static_cast<double>(box.bigEnd(0))/bfc + 2.;
        const int cyh = static_cast<double>(box.bigEnd(1))/bfc + 2.;
        const int czh = static_cast<double>(box.bigEnd(2))/bfc + 2.;

        for (int cxi = cxl; cxi <= cxh; ++cxi) {
            for (int cyi = cyl; cyi <= cyh; ++cyi) {
                for (int czi = czl; czi <= czh; ++czi) {
                    amrex::IntVect ct(wrapped_index[lev-1][cxi],
                                      wrapped_index[lev-1][cyi],
                                      wrapped_index[lev-1][czi]);

                    if (!bak.contains(ct)) {
                        layouts[lev-1][omp_get_thread_num()]->Add(
                            static_cast<double>(wrapped_index[lev-1][cxi])/bfc
                                    - 0.5,
                            static_cast<double>(wrapped_index[lev-1][cyi])/bfc
                                    - 0.5,
                            static_cast<double>(wrapped_index[lev-1][czi])/bfc
                                    - 0.5);
                    }
                }
            }
        }
    }

    FinalizeLayout(lev-1);
}

amrex::BoxArray LocalRegrid::WrapBoxArray(amrex::BoxArray& ba, int N) {
    amrex::BoxList new_bl;

    for (int i = -1; i <= 1; ++i) {
        for (int j =-1 ; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                amrex::Box bO(amrex::IntVect(i*N, j*N, k*N),
                              amrex::IntVect((i+1)*N-1, (j+1)*N-1, (k+1)*N-1));
                amrex::BoxList bl2w = ba.boxList().intersect(bO);
                if (bl2w.isNotEmpty()) {
                    bl2w.shift(0, -i*N);
                    bl2w.shift(1, -j*N);
                    bl2w.shift(2, -k*N);
                    new_bl.join(bl2w);
                }
            }
        }
    }

    return amrex::BoxArray(new_bl);
}

void LocalRegrid::AddBoxes(const int lev, amrex::BoxArray& ba) {
    // Create temporary distribution mapping, box array and multifab with only
    // the new boxes.
    amrex::DistributionMapping dm(ba, amrex::ParallelDescriptor::NProcs());
    amrex::MultiFab mf_new_tmp(ba, dm, sim->scalar_fields.size(), sim->nghost);
    amrex::MultiFab mf_old_tmp(ba, dm, sim->scalar_fields.size(), sim->nghost);

    // Fill temporary mf with data.
    sim->level_synchronizer->FillPatch(lev, sim->grid_new[lev].t, mf_new_tmp);
    sim->level_synchronizer->FillPatch(lev, sim->grid_old[lev].t, mf_old_tmp);

    // Create new joint box array.
    amrex::BoxList new_bl = sim->grid_new[lev].boxArray().boxList();
    new_bl.join(ba.boxList());
    amrex::BoxArray new_ba(std::move(new_bl));

    // Create new joint distribution mapping.
    amrex::Vector<int> new_pmap = sim->dmap[lev].ProcessorMap();
    amrex::Vector<int> pmap = dm.ProcessorMap();
    std::move(pmap.begin(), pmap.end(), std::back_inserter(new_pmap));
    amrex::DistributionMapping new_dm(new_pmap);

    // Create new MultiFab and fill it with data.
    LevelData new_mf(new_ba, new_dm, sim->scalar_fields.size(), sim->nghost,
                     amrex::MFInfo().SetAlloc(false));
    LevelData old_mf(new_ba, new_dm, sim->scalar_fields.size(), sim->nghost,
                     amrex::MFInfo().SetAlloc(false));

    const int offset = new_ba.size() - ba.size();
    for (amrex::MFIter mfi(new_mf); mfi.isValid(); ++mfi) {
        if (mfi.index() < offset) {
            amrex::FArrayBox& new_old_fab = sim->grid_new[lev][mfi.index()] ;
            new_mf.setFab(mfi, std::move(new_old_fab));

            amrex::FArrayBox& old_old_fab = sim->grid_old[lev][mfi.index()] ;
            old_mf.setFab(mfi, std::move(old_old_fab));
        } else {
            amrex::FArrayBox& new_old_fab = mf_new_tmp[mfi.index() - offset] ;
            new_mf.setFab(mfi, std::move(new_old_fab));

            amrex::FArrayBox& old_old_fab = mf_old_tmp[mfi.index() - offset] ;
            old_mf.setFab(mfi, std::move(old_old_fab));
        }
    }

    // Swap old MulftiFab with new one
    new_mf.t = sim->grid_new[lev].t;
    old_mf.t = sim->grid_old[lev].t;
    std::swap(sim->grid_new[lev], new_mf);
    std::swap(sim->grid_old[lev], old_mf);
    sim->SetBoxArray(lev, new_ba);
    sim->SetDistributionMap(lev, new_dm);
    sim->grid_old[lev].contains_truncation_errors = false;
}

}; // namespace sledgehamr
