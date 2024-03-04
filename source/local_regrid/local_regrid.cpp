#include "local_regrid.h"

namespace sledgehamr {

/** @brief Initalize all needed structures.
 * @param   owner   Pointer to simulation.
 */
LocalRegrid::LocalRegrid(Sledgehamr* owner) : sim(owner) {
    ParseInput();
    CreateCommMatrix();

    no_local_regrid.resize(sim->max_level+1, false);
    do_global_regrid.resize(sim->max_level+1, false);
    nregrids = 0;
}

/** @brief Perform a local regrid on level lev after checking various
 *         criteria.
 * @param   lev Coarsest level for tagging.
 * @return Whether we regridded or not.
 */
bool LocalRegrid::AttemptRegrid(const int lev) {
    amrex::Print() << std::endl;
    bool res = DoAttemptRegrid(lev);
    ClearLayout();
    return res;
}

/** @brief To be called if a global regrid on level lev has been performed.
 *         Needed in case this LocalRegrid module invoked the global regrid
 *         such that we can reset the relevant flags.
 * @param   lev Coarsest level on which we tagged for the global regrid.
 */
void LocalRegrid::DidGlobalRegrid(const int lev) {
    for (int l = 0; l <= sim->max_level; ++l) {
        no_local_regrid[l] = false;
        do_global_regrid[l] = false;
        nregrids = 0;
    }

    for (int l = lev+1; l < last_numPts.size(); ++l) {
        last_numPts[l] = sim->grid_new[l].boxArray().numPts();
        sim->grid_old[l].contains_truncation_errors = false;
    }
}

/** @brief Creates all layout structures needed for the local regrid.
 * @param   max_lev Up to what level we need the layout structures [deprecated].
 */
void LocalRegrid::InitializeLayout(const int max_lev) {
    layouts.resize(sim->finest_level + 1);
    for (int l = 1; l <= sim->finest_level; ++l) {
        int Np = sim->dimN[l] / sim->blocking_factor[l][0];
        for (int f = 0; f < omp_get_max_threads(); ++f) {
            layouts[l].emplace_back( std::make_unique<UniqueLayout>(this, Np) );
        }
    }
}

/** @brief Cleans up all layout structures after the local regrid.
 */
void LocalRegrid::ClearLayout() {
    layouts.clear();
}

/** @brief Joins BoxArray's from multiple nodes. Simplifies box lists only 
 *         locally for efficiency.
 * @param   lev Corresponding level.
 * @param   ba  Local BoxArray to be merged globally.
 */
void LocalRegrid::JoinBoxArrays(const int lev, amrex::BoxArray& ba) {
    amrex::BoxList bl = layouts[lev][0]->BoxList(sim->blocking_factor[lev][0]);
    bl.simplify(true);
    amrex::Vector<amrex::Box> bv = std::move(bl.data());
    amrex::AllGatherBoxes(bv);
    bl = amrex::BoxList(std::move(bv));
    bl.simplify(false);
    ba = amrex::BoxArray(std::move(bl));
}

/** @brief Adds a particular box to the layout structure.
 * @param   lev     Corresponding level.
 * @param   thread  Local OpenMP thread id.
 * @param   i       i-th box index.
 * @param   j       j-th box index.
 * @param   k       k-th box index.
 */
void LocalRegrid::AddToLayout(const int lev, const int thread, const int i,
                              const int j, const int k) {
    layouts[lev][thread]->Add(i, j, k);
}

/** @brief Synchronizes all layout structures across OpenMP threads and MPI
 *         ranks for a particular level.
 * @param   lev Corresponding level.
 */
void LocalRegrid::FinalizeLayout(const int lev) {
    layouts[lev][0]->Merge(layouts[lev]);
    layouts[lev][0]->Distribute();
}

/** @brief Checks pre-conditions whether we even need to attempt a local
 *         regrid.
 * @param   lev Coarsest level for tagging.
 * @return Whether we want to attempt a local regrid.
 */
bool LocalRegrid::Prechecks(const int lev) {
    if (nregrids++ >= max_local_regrids) {
        if (max_local_regrids > 0) {
            amrex::Print() << "Maximum number of local regrids reached: "
                           << max_local_regrids << std::endl;
        }
        return false;
    }

    // Veto if volume threshold such that it disables local regrid.
    if (volume_threshold_accumulated <= 1. ) {
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

    return true;
}

/** @brief Sets up everything needed prior to a local regrid.
 */
void LocalRegrid::InitializeLocalRegrid() {
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
}

/** @brief Determines new layout on all levels based on taggin criteria.
 * @param   lev Coarsest level for tagging.
 */
void LocalRegrid::DetermineAllBoxArrays(const int lev) {
    // Get the new required box array for each level. Still allowed to violate
    // nesting at this point.
    min_distance.clear();
    min_distance.resize(sim->finest_level + 1, -1.);
    for (int l = lev; l < sim->finest_level && !no_local_regrid[l]; ++l) {
        min_distance[l+1] = DetermineNewBoxArray(l);
    }
}

/** @brief Ensures that all levels are properly nested within the lower levels.
 */
void LocalRegrid::FixAllNesting() {
    // Make sure we do not violate nesting requirements.
    for (int l = sim->finest_level; l > 1; --l) {
        FixNesting(l);
    }
}

/** @brief Combines all BoxArray's on each level.
 * @param   box_arrays  Array of BoxArray's for each level.
 */
void LocalRegrid::JoinAllBoxArrays(std::vector<amrex::BoxArray>& box_arrays) {
    // Join all boxes across MPI ranks.
    for (int l = 1; l <= sim->finest_level; ++l) {
        JoinBoxArrays(l, box_arrays[l]);
    }
}

/** @brief Adds all boxes on all levels to the grid.
 * @param   box_arrays  Array of BoxArray's to be added to each level.
 */
void LocalRegrid::AddAllBoxes(std::vector<amrex::BoxArray>& box_arrays) {
    // Finally add new boxes to each grid.
    for (int l = 1; l <= sim->finest_level; ++l) {
        if (box_arrays[l].size() > 0) {
            AddBoxes(l, box_arrays[l]);
        }
    }
}

/** @brief Checks all thresholds and determines whether we should continue
 *         with the local regrid.
 * @param   lev Level to be checked.
 * @param   ba  BoxArray to be added to the given level.
 * @return Whether we should veto the local regrid.
 */
bool LocalRegrid::CheckThresholds(const int lev, amrex::BoxArray& ba) {
    bool veto = false;

    double Nb = ba.numPts();
    double Nc = sim->grid_new[lev].boxArray().numPts();
    double Nr = last_numPts[lev];

    if (Nr == 0)
        Nr = Nc;

    double dV = Nb / Nc;
    double fV = (Nb+Nc) / Nr;

    // This level exceeds strong threshold so we veto.
    if (fV > volume_threshold_accumulated)
        veto = true;

    // In case any level vetos the local regrid, this is the level on which
    // we want to perform the global regrid.
    if( fV > volume_threshold_single && veto_level == -1 )
        veto_level = lev - 1;

    amrex::Print() << "  Additional boxes on level " << lev << " required: "
                   << ba.size() << std::endl
                   << "    Instantanous volume increase: " << dV
                   << std::endl
                   << "    Volume increase since last global regrid: "
                   << fV << ". Threshold: " << volume_threshold_accumulated
                   << std::endl;

    return veto;
}

/** @brief Computes the most pessimistiv time by which we should have performed
 *         any regrid or otherwise we run risk of having features propagate
 *         outside a refined region.
 * @param   l   Level to be used for computation.
 * @param   lev Coarsest level on which we tagged.
 */
void LocalRegrid::ComputeLatestPossibleRegridTime(const int l, const int lev) {
    if (l <= lev)
        return;

    double dx_c = n_error_buf;
    double regrid_dt = sim->time_stepper->regrid_dt[l];
    double dt_delay = DBL_MAX;
    if (min_distance[l] >= 0)
        dt_delay = min_distance[l] / dx_c * regrid_dt;

    latest_possible_regrid_time[l] = sim->grid_new[l].t + dt_delay;

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

/** @brief Indentifies whether any part of the code wants to veto the local
 *         regrid.
 * @param   lev         Level to check for vetos.
 * @param   box_arrays  Array of BoxArray's containing the local additions.
 * @return Whether we found a veto anywhere.
 */
bool LocalRegrid::CheckForVeto(const int lev,
                               std::vector<amrex::BoxArray>& box_arrays) {
    bool veto = false;
    latest_possible_regrid_time.clear();
    latest_possible_regrid_time.resize(sim->finest_level+1, -1.);

    for (int l = 1; l <= sim->finest_level; ++l) {
        veto = veto || CheckThresholds(l, box_arrays[l]);
        ComputeLatestPossibleRegridTime(l, lev);
    }

    return veto;
}

/** @brief Determines the best course of action in case the local regrid was
 *         vetoed. This could be either to do nothing, do a global regrid, or
 *         or to do a local regrid anyway followed by a global regrid at the
 *         next opportunity.
 * @param   lev Coarsest level for tagging.
 * @return Suggested course of action (no/local/global regrid).
 */
LocalRegrid::VetoResult LocalRegrid::DealWithVeto(const int lev) {
    amrex::Print() << "Local regrid has been vetoed. "
                   << "Global regrid on level " << veto_level
                   << " (adjusting level " << veto_level+1<<") "
                   << "deemed optimal." << std::endl;

    // We can do a globl regrid right away.
    if (veto_level >= lev)
        return VetoResult::DoGlobalRegrid;

    // Since we only create a shadow level when needed triggering a global
    // regrid at this level is non-trivial. For now, just give up and do
    // global at higher levels and wait for the coarse level when it is
    // scheduled normally.
    if (veto_level == 0)
        return VetoResult::DoGlobalRegrid;

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

        for (int l = lev; l <= sim->finest_level; ++l) {
            no_local_regrid[l] = true;
        }

        return VetoResult::DoNoRegrid;
    }

    return VetoResult::DoLocalRegrid;
}

/** @brief Attempts a local regrid.
 * @param   lev Coarsest level for tagging.
 * @return Whether the local regrid was succesfull. A local regrid is considered
 *         successfull even if it has not been performed but no further action
  *        has been deemed necessary.
 */
bool LocalRegrid::DoAttemptRegrid(const int lev) {
    if (!Prechecks(lev)) return false;

    amrex::Print() << std::endl << "Attempting local regrid at level " << lev+1
                   << " and higher." << std::endl;

    InitializeLocalRegrid();
    DetermineAllBoxArrays(lev);
    FixAllNesting();
    std::vector<amrex::BoxArray> box_arrays(sim->finest_level+1);
    JoinAllBoxArrays(box_arrays);

    if( CheckForVeto(lev, box_arrays) ) {
        switch (DealWithVeto(lev)) {
            case VetoResult::DoGlobalRegrid:
                return false;
            case VetoResult::DoNoRegrid:
                return true;
            case VetoResult::DoLocalRegrid:
                //[[fallthrough]];
            default:
                break;
        }
    }

    AddAllBoxes(box_arrays);

    return true;
}

/** @brief Parses all input parameters related to a local regrid.
 */
void LocalRegrid::ParseInput() {
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("force_global_regrid_at_restart",
                  force_global_regrid_at_restart);
    pp_amr.query("n_error_buf", n_error_buf);
    pp_amr.query("max_local_regrids", max_local_regrids);
    pp_amr.query("volume_threshold_strong", volume_threshold_single);
    pp_amr.query("volume_threshold_weak", volume_threshold_accumulated);
}

/** @brief Creates a N x N communication matrix M_{ij} between N MPI ranks. The
 *         i-th MPI rank talks to the M_{ij}-th MPI rank during the j-th
 *         communication cycle. Each rank talks to each other rank exactly once.
 */
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

/** @brief Wraps box indicies across the periodic boundary conditions. 
 * @param   lev Level on which to perform wrapping.
 */
void LocalRegrid::WrapIndices(const int lev) {
    if (wrapped_index.size() != lev)
        return;

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

/** @brief Determines whether a box has neighbouring coarse/fine boundaries.
 * @param   tilebox Current box to check.
 * @param   c0      Tilebox lower box index.
 * @param   c1      Tilebox upper box index.
 * @param   lev     Current level.
 * @param   border  Boolean map indicating whether a coarse/fine boundary is
 *                  present.
 * @return Total number of neighbouring coarse fine boundaries.
 */
inline int LocalRegrid::GetBoxCoarseFineBorders(
        const amrex::Box& tilebox, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int lev,
        boost::multi_array<bool, 3>& border) {
    const amrex::BoxArray& ba = sim->grid_new[lev+1].boxArray();

    const amrex::IntVect xs = c1 - c0 + 1;
    const amrex::IntVect xsb = xs + 2;

    int remaining = 0;
    border.resize(boost::extents[xsb[0]][xsb[1]][xsb[2]]);
    for (int k=-1; k<=xs[2]; ++k) {
        for (int j=-1; j<=xs[1]; ++j) {
            for (int i=-1; i<=xs[0]; ++i) {
                border[i+1][j+1][k+1] = !ba.contains( amrex::IntVect(
                        wrapped_index[lev+1][c0[0]+i],
                        wrapped_index[lev+1][c0[1]+j],
                        wrapped_index[lev+1][c0[2]+k]));

                remaining += border[i+1][j+1][k+1];
            }
        }
    }

    return remaining;
}

/** @brief Finds tagged cells within a box and measures its distance of any
 *         coarse/fine boundary.
 * @param   lo          Lower coordinate of current box.
 * @param   hi          Higher coordinate of current box.
 * @param   remaining   Number of adjecent coarse/fine boundaries that could
 *                      still be pushed further.
 * @param   tag_arr     TagArray containing all tags.
 * @param   c0          Lower box index of current box.
 * @param   c1          Upper box index of current box.
 * @param   lev         Current coarse level.
 * @param   ibff        Fine coarse level blocking factor of type int.
 * @param   bff         Fine coarse level blocking factor but cast to double.
 * @param   border      Map that indicates locations of coarse/fine boundaries.
 * @param   closest_locations   Array containing the locations closest to a
 *                              coarse/fine boundary on each OpenMP rank.                             
 * @param   threshold           Distance threshold beyond which we don't really
 *                              care about the exact distance anymore.
 * @param   omp_thread_num      Current OpenMP thread.
 */
inline void LocalRegrid::TagAndMeasure(
        const amrex::Dim3& lo, const amrex::Dim3& hi, int remaining,
        const amrex::Array4<char>& tag_arr, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int lev, const int ibff,
        const double bff, boost::multi_array<bool, 3>& border,
        std::vector<Location>& closest_locations, const double threshold,
        const int omp_thread_num) {
    // Now find tags and determine distances.
    for (int k = lo.z; k <= hi.z && remaining > 0; ++k) {
        for (int j = lo.y; j <= hi.y && remaining > 0; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (remaining == 0)
                    break;

                if (tag_arr(i,j,k) != amrex::TagBox::SET)
                    continue;

                amrex::IntVect ci(i,j,k);
                CheckBorders(ci, c0, c1, ibff, bff, remaining, lev, border,
                             closest_locations, threshold, omp_thread_num);
            }
        }
    }
}

/** @brief Given a location we measure the distance to all adjacent coarse/fine
 *         boundaries. We add to our layout structure if the location is too 
 *         close.
 * @param   ci  Index to be checked.
 * @param   c0          Lower box index of current box.
 * @param   c1          Upper box index of current box.
 * @param   ibff        Fine coarse level blocking factor of type int.
 * @param   bff         Fine coarse level blocking factor but cast to double.
 * @param   remaining   Number of adjecent coarse/fine boundaries that could
 *                      still be pushed further.
 * @param   lev         Current coarse level.
 * @param   border      Map that indicates locations of coarse/fine boundaries.
 * @param   closest_locations   Array containing the locations closest to a
 *                              coarse/fine boundary on each OpenMP rank.                             
 * @param   threshold           Distance threshold beyond which we don't really
 *                              care about the exact distance anymore.
 * @param   omp_thread_num      Current OpenMP thread.
 */
inline void LocalRegrid::CheckBorders(
        const amrex::IntVect& ci, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int ibff, const double bff,
        int remaining, const int lev, boost::multi_array<bool, 3>& border,
        std::vector<Location>& closest_locations, const double threshold,
        const int omp_thread_num) {
    amrex::IntVect fi = ci*2;
    amrex::IntVect cfi(static_cast<int>(static_cast<double>(fi[0])/bff)+1,
                       static_cast<int>(static_cast<double>(fi[1])/bff)+1,
                       static_cast<int>(static_cast<double>(fi[2])/bff)+1);

    // Check each possible border.
    for (int kk=-1; kk<=1; ++kk) {
        for (int jj=-1; jj<=1; ++jj) {
            for (int ii=-1; ii<=1; ++ii) {
                amrex::IntVect cii(ii, jj, kk);
                amrex::IntVect ref = cfi - c0 + cii + 1;

                if (!border[ref[0]][ref[1]][ref[2]])
                    continue;

                amrex::IntVect smt = ibff*(cfi + cii - 1);
                amrex::IntVect bgt = ibff*(cfi + cii) - 1;

                int d_sq = 0;
                for(int d=0; d < 3; d++) {
                    d_sq += smt[d] < fi[d] && fi[d] < bgt[d] ? 0 :
                                std::min((fi[d] - smt[d])*(fi[d] - smt[d]),
                                         (fi[d] - bgt[d])*(fi[d] - bgt[d]));
                }

                closest_locations[omp_thread_num].SelectClosest(
                        ci[0], ci[1], ci[2], d_sq);

                // Add new box if below threshold.
                if (d_sq < threshold) {
                    border[ref[0]][ref[1]][ref[2]] = false;
                    remaining--;
                    layouts[lev+1][omp_thread_num]->Add(
                        static_cast<int>(static_cast<double>(
                        wrapped_index[lev+1][cfi[0]+ii])/bff-0.5),
                        static_cast<int>(static_cast<double>(
                        wrapped_index[lev+1][cfi[1]+jj])/bff-0.5),
                        static_cast<int>(static_cast<double>(
                        wrapped_index[lev+1][cfi[2]+kk])/bff-0.5));
                }
            }
        }
    }
}

/** @brief Tags cells on a level and measure distance to coarse/fine boundary.
 *         If the tagged cell is too close we expand our layout structure at 
 *         that location.
 * @param   lev Level on which to tag.
 * @return Globally shortest distance between a tagged cell and coarse/fine
 *         boundary. Will be -1 if nothing has been tagged or distance is too
 *         large to reliably determine distance (> blocking_factor).
 */
double LocalRegrid::DetermineNewBoxArray(const int lev) {
    const double threshold = (n_error_buf+1) * (n_error_buf+1);
    const double dimNf = sim->dimN[lev+1];
    const amrex::BoxArray& ba = sim->grid_new[lev+1].boxArray();
    const LevelData& state = sim->grid_new[lev];
    const int ibff = sim->blocking_factor[lev+1][0];
    const double bff = static_cast<double>(ibff);

    // Tag cells.
    amrex::TagBoxArray tags(state.boxArray(), state.DistributionMap());
    sim->ErrorEst(lev, tags, state.t, 0);

    std::vector<Location> closest_locations(omp_get_max_threads());

#pragma omp parallel
    for (amrex::MFIter mfi(state, false); mfi.isValid(); ++mfi) {
        const amrex::Array4<double const>& state_fab = state.array(mfi);
        const amrex::Array4<char>& tag_arr = tags.array(mfi);

        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Dim3 lo = amrex::lbound(tilebox);
        const amrex::Dim3 hi = amrex::ubound(tilebox);
        const amrex::IntVect c0(
                static_cast<int>(static_cast<double>(lo.x*2)/bff) + 1,
                static_cast<int>(static_cast<double>(lo.y*2)/bff) + 1,
                static_cast<int>(static_cast<double>(lo.z*2)/bff) + 1);
        const amrex::IntVect c1(
                static_cast<int>(static_cast<double>(hi.x*2)/bff) + 1,
                static_cast<int>(static_cast<double>(hi.y*2)/bff) + 1,
                static_cast<int>(static_cast<double>(hi.z*2)/bff) + 1);

        // Determine if tilebox is near a C/F border.
        boost::multi_array<bool, 3> border;
        int remaining = GetBoxCoarseFineBorders(tilebox, c0, c1, lev, border);
        if (remaining == 0)
            continue;

        TagAndMeasure(lo, hi, remaining, tag_arr, c0, c1, lev, ibff, bff,
                      border, closest_locations, threshold,
                      omp_get_thread_num());
    }


    // Combine box layouts.
    FinalizeLayout(lev+1);
    Location closest = Location::FindClosestGlobally(closest_locations);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        if( closest.distance_sq < 46000) {
            amrex::Print() << "  Shortest distance to C/F boundary: "
                           << sqrt(static_cast<double>(closest.distance_sq))
                           << " grid sites @ (" << closest.i << ","
                           << closest.j << "," << closest.k << ")"
                           << std::endl;
        } else {
            amrex::Print() << "  Shortest distance to C/F boundary: "
                           << "> blocking_factor (" << ibff << " cells)"
                           << std::endl;
        }
    }

    if (closest.i >= 0 )
        return sqrt(static_cast<double>(closest.distance_sq));
    else
        return -1;
}

/** @brief Will expand lower level locally to ensure the fine level is nested
 *         properly within.
 * @param   lev Fine level.
 */
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

/** @brief Wraps a box array accross periodic boundary conditions.
 * @param   ba  BoxArray.
 * @param   N   Total number of potential boxes along an axis.
 */
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

/** @brief Explicitly add boxes to a given level and fill it with data.
 * @param   lev Current level.
 * @param   ba  BoxArray to be added to level.
 */
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
    new_mf.istep = sim->grid_new[lev].istep;
    old_mf.istep = sim->grid_old[lev].istep;
    std::swap(sim->grid_new[lev], new_mf);
    std::swap(sim->grid_old[lev], old_mf);
    sim->SetBoxArray(lev, new_ba);
    sim->SetDistributionMap(lev, new_dm);
    sim->grid_old[lev].contains_truncation_errors = false;
}

}; // namespace sledgehamr
