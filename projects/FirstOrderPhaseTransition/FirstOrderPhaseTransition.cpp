#include "FirstOrderPhaseTransition.h"

namespace FirstOrderPhaseTransition{

void FirstOrderPhaseTransition::Init() {
    ParseVariables();
    ParseBubbles();
    ComputeParameters();
    SetProjections();

    idx_perfmon_add_bubbles = performance_monitor->timer.size();
    performance_monitor->timer.push_back(sledgehamr::Timer("InjectBubbles"));
}

void FirstOrderPhaseTransition::ComputeParameters() {
    double numer = 3. + sqrt(9. - 8*lambda_bar);
    quadratic = -1.;
    cubic = 3.*numer/(4.*lambda_bar);
    quartic = -numer*numer/(8.*lambda_bar);
}

void FirstOrderPhaseTransition::SetProjections() {
    sledgehamr::Projection proj(dPhi2, "dPhi2");
    io_module->projections.push_back(proj);
}

void FirstOrderPhaseTransition::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("lambda_bar", lambda_bar);

    pp_prj.queryarr("bubbles_to_inject", bubbles_to_inject);
}


void FirstOrderPhaseTransition::SetParamsRhs(
        std::vector<double>& params, const double time, const int lev) {
    params.resize(3);
    params[0] = quadratic;
    params[1] = cubic;
    params[2] = quartic;
}

void FirstOrderPhaseTransition::BeforeTimestep(const double time) {
    InjectBubbles(time);
}

void FirstOrderPhaseTransition::ParseBubbles() {
    std::string file = "";
    amrex::ParmParse pp("input");
    pp.get("bubbles", file);

    if( file == "")
        return;

    amrex::Print() << "Read bubble information: " << file << std::endl;

    double header[5];
    sledgehamr::IOModule::ReadFromHDF5(file, {"Header"}, header);

    int NB = header[2];
    if( NB == 0 )
        return;

    std::vector<double> xlocs(NB);
    std::vector<double> ylocs(NB);
    std::vector<double> zlocs(NB);
    std::vector<double> ts(NB);
    std::vector<int> use_profile(NB);

    sledgehamr::IOModule::ReadFromHDF5(file, {"xlocs"}, &xlocs[0]);
    sledgehamr::IOModule::ReadFromHDF5(file, {"ylocs"}, &ylocs[0]);
    sledgehamr::IOModule::ReadFromHDF5(file, {"zlocs"}, &zlocs[0]);
    sledgehamr::IOModule::ReadFromHDF5(file, {"t"}, &ts[0]);
    sledgehamr::IOModule::ReadFromHDF5(file, {"use_profile"}, &use_profile[0]);

    const int ncomp = grid_new[0].nComp();
    for (int b = 0; b < NB; ++b) {
        Bubble B;
        B.x = xlocs[b];
        B.y = ylocs[b];
        B.z = zlocs[b];
        B.t = ts[b];

        if (use_profile[b] == b) {
            double profile_header[4];
            std::string str_profile_header = "profile_header_"
                                           + std::to_string(b);
            sledgehamr::IOModule::ReadFromHDF5(file, {str_profile_header},
                                               profile_header);

            B.inv_dx = profile_header[1];
            B.L = profile_header[2];
            B.finest_level = profile_header[3];

            B.level.resize((int)profile_header[0]);
            std::string str_data_level = "profile_level_" + std::to_string(b);
            sledgehamr::IOModule::ReadFromHDF5(file, {str_data_level}, &B.level[0]);

            for (int n = 0; n < ncomp; ++n) {
                std::string sname = scalar_fields[n]->name;
                // TODO Remove old convention.
                if (sname == "Phi")  {
                    sname = "Psi1";
                } else if (sname == "dPhi") {
                    sname = "Pi1";
                } else {
                    continue;
                }
                std::string str_data = "profile_" + sname + "_"
                                     + std::to_string(b);

                std::vector<double> profile((int)profile_header[0]);
                sledgehamr::IOModule::ReadFromHDF5(file, {str_data},
                                                   &profile[0]);
                B.profile.push_back( profile );
            }
        }

        bubbles.push_back(B);
    }

    for (int b = 0; b < NB; ++b){
        bubbles[b].p_bubble = &(bubbles[use_profile[b]]);
    }

    //std::sort( bubbles.begin(), bubbles.end() );
    amrex::Print() << bubbles.size() << " bubble(s) found to be injected."
                   << std::endl;
}

void FirstOrderPhaseTransition::InjectBubbles(const double time) {
    performance_monitor->Start(idx_perfmon_add_bubbles);

    // Check if bubbles need to be added.
    std::vector<int> ab;
    int skip = 0;
    for (int b = next_bubble; b < bubbles.size(); ++b) {
        if (bubbles[b].t > time)
            continue;

        if (bubbles_to_inject.size() > 0) {
            amrex::Print() << "Test: " << b << std::endl;
            for(int x=0;x<bubbles_to_inject.size();++x)
                amrex::Print() << bubbles_to_inject[x] << std::endl;

            if(std::find(bubbles_to_inject.begin(), bubbles_to_inject.end(), b)
                == bubbles_to_inject.end()) {
                continue;
            }
        }

        if (time - bubbles[b].t < dt[0]) {
            ab.push_back(b);
        } else {
            skip++;
        }

        next_bubble = b + 1;
    }

    if (skip > 0) {
        amrex::Print() << "Skipping " << skip
                       << " bubble(s) that have been injected earlier already."
                       << std::endl;
    }

    performance_monitor->Stop(idx_perfmon_add_bubbles);

    if (ab.size() == 0)
        return;

    amrex::Print() << "Injecting " << ab.size() << " bubble(s) ... "
                   << std::endl;

    performance_monitor->Start(idx_perfmon_add_bubbles);

    InjectBubbleLevels(ab);
    AddBubbleValues(ab);

    // Make levels consistent
    for (int lev = finest_level - 1; lev >= 0; lev--)
        level_synchronizer->AverageDownTo(lev);

    performance_monitor->Stop(idx_perfmon_add_bubbles);
}

void FirstOrderPhaseTransition::InjectBubbleLevels(std::vector<int> ab)
{
    // Check what finest injection level is.
    int finest_bubble_level = 0;
    for (int b : ab) {
        finest_bubble_level = std::max(finest_bubble_level,
                                       bubbles[b].GetFinestLevel());
    }

    // Stop if no level needs to be injected.
    if (finest_bubble_level < 1)
        return;

    // Initialize entirely new levels if needed to zero for local regrid.
    for (int lev = finest_level + 1; lev <= finest_bubble_level; ++lev) {
        const int ncomp = grid_new[0].nComp();
        const int nghost = grid_new[0].nGrow();
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;

        grid_new[lev].define(ba, dm, ncomp, nghost, grid_new[0].t);
        grid_old[lev].define(ba, dm, ncomp, nghost);

        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);
    }

    int old_finest_level = finest_level;
    finest_level = finest_bubble_level;

    // Initialize layout structure from local regrid.
    time_stepper->local_regrid->InitializeLayout(finest_bubble_level);

    for (int lev = 0; lev <= finest_bubble_level; ++lev) {
        time_stepper->local_regrid->WrapIndices(lev);
    }

    // Fill layout structure for each level.
    for (int lev = 1; lev <= finest_bubble_level; ++lev) {
        FillBubbleLayout(lev, ab);
    }

    // Ensure nesting.
    for (int lev = finest_level; lev > 1; lev--) {
        time_stepper->local_regrid->FixNesting(lev);
    }

    // Join all boxes across MPI ranks.
    std::vector<amrex::BoxArray> box_arrays(finest_level+1);
    for (int lev = 1; lev <= finest_level; ++lev) {
        time_stepper->local_regrid->JoinBoxArrays(lev, box_arrays[lev]);
    }

    // Add boxes to level for existing levels.
    for (int lev = 1; lev <= finest_level; ++lev) {
        if (box_arrays[lev].size() > 0) {
            time_stepper->local_regrid->AddBoxes(lev, box_arrays[lev]);
            geom[lev] = geom[lev-1];
            geom[lev].refine(amrex::IntVect(2,2,2));
        }
    }

    time_stepper->local_regrid->ClearLayout();
}

void FirstOrderPhaseTransition::FillBubbleLayout(const int lev,
                                                    std::vector<int> ab) {
    const amrex::BoxArray ba = grid_new[lev].boxArray();
    int Nbs = dimN[lev] / blocking_factor[lev][0];
    double dxb = L / (double)Nbs;

    // For now all nodes do the same work. Could be distributed, but workload
    // probably not too bad either way.
#pragma omp parallel for
    for (int i = 0; i < Nbs; ++i) {
        for (int j = 0; j < Nbs; ++j) {
            for (int k = 0; k < Nbs; ++k) {
                // Check if already exists.
                amrex::IntVect ce(i,j,k);
                if (ba.contains(ce*blocking_factor[lev]))
                    continue;

                // Check all boxes.
                for (int b = 0; b < ab.size(); ++b) {
                    Bubble *bub = &(bubbles[ab[b]]);
                    double minD = DBL_MAX;
                    double maxD = 0;
                    const double Dx[2] = {Distance(i*dxb,     bub->x, L),
                                          Distance((i+1)*dxb, bub->x, L)};
                    const double Dy[2] = {Distance(j*dxb,     bub->y, L),
                                          Distance((j+1)*dxb, bub->y, L)};
                    const double Dz[2] = {Distance(k*dxb,     bub->z, L),
                                          Distance((k+1)*dxb, bub->z, L)};

                    for (int ii = 0; ii <= 1; ++ii) {
                        for (int jj = 0; jj <= 1; ++jj) {
                            for (int kk = 0; kk <= 1; ++kk) {
                                double D = std::sqrt(Dx[ii]*Dx[ii]
                                                   + Dy[jj]*Dy[jj]
                                                   + Dz[kk]*Dz[kk]);
                                minD = std::min(D,minD);
                                maxD = std::max(D,maxD);
                            }
                        }
                    }

                    int ind0 = bub->GetPos(minD);
                    if (ind0 == -1)
                        continue;

                    int ind1 = bub->GetPos(maxD);
                    if (ind1 == -1 )
                        ind1 = bub->GetNBins() - 1;

                    for (int ind = ind0; ind <= ind1; ++ind) {
                        if (bub->GetLevel(ind) >= lev) {
                            time_stepper->local_regrid->AddToLayout(
                                    lev, omp_get_thread_num(), i, j, k);
                            b = ab.size();
                            break;
                        }
                    }
                }
            }
        }
    }

    time_stepper->local_regrid->FinalizeLayout(lev);
}

void FirstOrderPhaseTransition::AddBubbleValues(std::vector<int> ab) {
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::MultiFab& mf = grid_new[lev];

#pragma omp parallel
        for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            const auto& fab = mf.array(mfi);

            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        for(int b=0;b<ab.size();++b) {
                            AddBubble(i, j, k, dx[lev], L, fab, bubbles[ab[b]]);
                        }
                    }
                }
            }
        }
    }
}

}; // namespace FirstOrderPhaseTransition
