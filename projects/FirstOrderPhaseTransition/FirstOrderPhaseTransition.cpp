#include <hdf5_utils.h>

#include "FirstOrderPhaseTransition.h"
#include "spectrum_modifier.h"

namespace FirstOrderPhaseTransition {

void FirstOrderPhaseTransition::Init() {
    ParseVariables();
    ParseBubbles();
    ComputeParameters();
    SetProjections();
    AddSpectrumModification();

    idx_perfmon_add_bubbles = performance_monitor->timer.size();
    performance_monitor->timer.emplace_back("InjectBubbles");
}

void FirstOrderPhaseTransition::ComputeParameters() {
    double numer = 3. + sqrt(9. - 8 * lambda_bar);
    quadratic = -1.;
    cubic = 3. * numer / (4. * lambda_bar);
    quartic = -numer * numer / (8. * lambda_bar);
}

void FirstOrderPhaseTransition::SetProjections() {
    io_module->projections.emplace_back(dPhi2, "dPhi2");
}

void FirstOrderPhaseTransition::ParseVariables() {
    amrex::ParmParse pp_prj("project");

    pp_prj.query("potential_type", potential_type);
    switch (potential_type) {
    case PotentialType::PureLambdaBar:
        pp_prj.get("lambda_bar", lambda_bar);
        break;
    case PotentialType::Piecewise:
        pp_prj.get("lambda_bar", lambda_bar);
        pp_prj.get("vbar", vbar);
        pp_prj.get("vareps", vareps);
        pp_prj.get("phiesc", phiesc);
        break;
    default:
        amrex::Abort("Unkown potential type!");
    }

    pp_prj.queryarr("bubbles_to_inject", bubbles_to_inject);
    pp_prj.query("tc", tc);
    pp_prj.query("t0", t0);
}

void FirstOrderPhaseTransition::AddSpectrumModification() {
    io_module->output.emplace_back(
        "gw_spec_u_times_k",
        OUTPUT_FCT(FirstOrderPhaseTransition::GwSpectrum_UtimesK));

    io_module->output.emplace_back(
        "gw_spec_two_bubbles_from_one",
        OUTPUT_FCT(FirstOrderPhaseTransition::GwSpectrum_2BubblesFrom1));
}

bool FirstOrderPhaseTransition::GwSpectrum_UtimesK(double time,
                                                   std::string prefix) {
    if (!with_gravitational_waves)
        return false;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    std::unique_ptr<sledgehamr::GravitationalWavesSpectrumModifier> modifier =
        std::make_unique<SpectrumModifier_UtimesK>();
    gravitational_waves->ComputeSpectrum(file_id, modifier.get());

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

bool FirstOrderPhaseTransition::GwSpectrum_2BubblesFrom1(double time,
                                                         std::string prefix) {
    if (!with_gravitational_waves || bubbles.size() < 2)
        return false;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    const double d[3] = {bubbles[1].x - bubbles[0].x,
                         bubbles[1].y - bubbles[0].y,
                         bubbles[1].z - bubbles[0].z};

    std::unique_ptr<sledgehamr::GravitationalWavesSpectrumModifier> modifier =
        std::make_unique<SpectrumModifier_2BubblesFrom1>(d);
    gravitational_waves->ComputeSpectrum(file_id, modifier.get());

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

void FirstOrderPhaseTransition::SetParamsRhs(std::vector<double> &params,
                                             const double time, const int lev) {
    params.resize(7);
    params[0] = static_cast<double>(potential_type);
    params[1] = quadratic;
    params[2] = cubic;
    params[3] = quartic;
    params[4] = vbar;
    params[5] = vareps;
    params[6] = phiesc;
}

void FirstOrderPhaseTransition::SetParamsGravitationalWaveRhs(
    std::vector<double> &params, const double time, const int lev) {
    params.resize(2);
    params[0] = tc;
    params[1] = t0;
}

void FirstOrderPhaseTransition::SetParamsTruncationModifier(
    std::vector<double> &params, const double time, const int lev) {
    if (maxima_time != time) {
        sledgehamr::LevelData &ld = grid_new[0];
        const int ncomp = ld.nComp();

        if (field_maxima.size() != ncomp) {
            field_maxima.resize(ncomp);
            comp_vector.resize(ncomp);
            for (int c = 0; c < ncomp; ++c) {
                comp_vector[c] = c;
            }
        }

        field_maxima = ld.norm0(comp_vector);
        maxima_time = time;

        for (int c = 0; c < ncomp; ++c) {
            amrex::Print() << "Maximum field value of " << GetScalarFieldName(c)
                           << ": " << field_maxima[c] << std::endl;
        }
    }

    params = field_maxima;
}

void FirstOrderPhaseTransition::BeforeTimestep(const double time) {
    InjectBubbles(time);
}

void FirstOrderPhaseTransition::ParseBubbles() {
    std::string file = "";
    amrex::ParmParse pp("input");
    pp.get("bubbles", file);

    if (file == "")
        return;

    amrex::Print() << "Read bubble information: " << file << std::endl;

    double header[5];
    sledgehamr::utils::hdf5::Read(file, {"Header"}, header);

    int NB = header[2];
    if (NB == 0)
        return;

    std::vector<double> xlocs(NB);
    std::vector<double> ylocs(NB);
    std::vector<double> zlocs(NB);
    std::vector<double> ts(NB);
    std::vector<int> use_profile(NB);

    sledgehamr::utils::hdf5::Read(file, {"xlocs"}, &xlocs[0]);
    sledgehamr::utils::hdf5::Read(file, {"ylocs"}, &ylocs[0]);
    sledgehamr::utils::hdf5::Read(file, {"zlocs"}, &zlocs[0]);
    sledgehamr::utils::hdf5::Read(file, {"t"}, &ts[0]);
    sledgehamr::utils::hdf5::Read(file, {"use_profile"}, &use_profile[0]);

    const int ncomp = grid_new[0].nComp();
    for (int b = 0; b < NB; ++b) {
        Bubble B;
        // Intentional flip (temporary).
        B.x = zlocs[b] - zlocs[0];
        B.y = ylocs[b] - ylocs[0];
        B.z = xlocs[b] - xlocs[0];
        B.t = ts[b];

        if (use_profile[b] == b) {
            double profile_header[4];
            std::string str_profile_header =
                "profile_header_" + std::to_string(b);
            sledgehamr::utils::hdf5::Read(file, {str_profile_header},
                                          profile_header);

            B.inv_dx = profile_header[1];
            B.L = profile_header[2];
            B.finest_level = profile_header[3];

            B.level.resize((int)profile_header[0]);
            std::string str_data_level = "profile_level_" + std::to_string(b);
            sledgehamr::utils::hdf5::Read(file, {str_data_level}, &B.level[0]);

            for (int n = 0; n < ncomp; ++n) {
                std::string sname = scalar_fields[n]->name;
                // TODO Remove old convention.
                if (sname == "Phi") {
                    sname = "Psi1";
                } else if (sname == "dPhi") {
                    sname = "Pi1";
                } else {
                    continue;
                }
                std::string str_data =
                    "profile_" + sname + "_" + std::to_string(b);

                std::vector<double> profile((int)profile_header[0]);
                sledgehamr::utils::hdf5::Read(file, {str_data}, &profile[0]);
                B.profile.push_back(profile);
            }
        }

        bubbles.push_back(B);
    }

    for (int b = 0; b < NB; ++b) {
        bubbles[b].p_bubble = &(bubbles[use_profile[b]]);
    }

    MoveBubblesToCentre();

    // std::sort( bubbles.begin(), bubbles.end() );
    amrex::Print() << bubbles.size() << " bubble(s) found to be injected."
                   << std::endl;
}

void FirstOrderPhaseTransition::MoveBubblesToCentre() {
    if (bubbles.size() > 2) {
        return;
    }

    double C = L / 2.;

    if (bubbles.size() == 1) {
        bubbles[0].x = C;
        bubbles[0].y = C;
        bubbles[0].z = C;
        return;
    }

    double cx = (bubbles[0].x + bubbles[1].x) / 2.;
    double cy = (bubbles[0].y + bubbles[1].y) / 2.;
    double cz = (bubbles[0].z + bubbles[1].z) / 2.;
    for (Bubble &bubble : bubbles) {
        bubble.x += C - cx;
        bubble.y += C - cy;
        bubble.z += C - cz;
    }
}

std::vector<int> FirstOrderPhaseTransition::FindBubbles(const double time) {
    std::vector<int> ab;
    int skip = 0;
    for (int b = next_bubble; b < bubbles.size(); ++b) {
        if (bubbles[b].t > time)
            continue;

        if (bubbles_to_inject.size() > 0) {
            for (int x = 0; x < bubbles_to_inject.size(); ++x)
                amrex::Print() << bubbles_to_inject[x] << std::endl;

            if (std::find(bubbles_to_inject.begin(), bubbles_to_inject.end(),
                          b) == bubbles_to_inject.end()) {
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

    return ab;
}

void FirstOrderPhaseTransition::InjectBubbles(const double time) {
    performance_monitor->Start(idx_perfmon_add_bubbles);

    std::vector<int> ab = FindBubbles(time);
    if (ab.size() == 0) {
        performance_monitor->Stop(idx_perfmon_add_bubbles);
        return;
    }

    amrex::Print() << "Injecting " << ab.size() << " bubble(s) ... "
                   << std::endl;

    InjectBubbleLevels(ab);
    AddBubbleValues(ab);

    // Make levels consistent
    for (int lev = finest_level - 1; lev >= 0; lev--)
        level_synchronizer->AverageDownTo(lev);

    performance_monitor->Stop(idx_perfmon_add_bubbles);
}

void FirstOrderPhaseTransition::InjectBubbleLevels(std::vector<int> ab) {
    // Check what finest injection level is.
    int finest_bubble_level = 0;
    for (int b : ab) {
        finest_bubble_level =
            std::max(finest_bubble_level, bubbles[b].GetFinestLevel());
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
    std::vector<amrex::BoxArray> box_arrays(finest_level + 1);
    for (int lev = 1; lev <= finest_level; ++lev) {
        time_stepper->local_regrid->JoinBoxArrays(lev, box_arrays[lev]);
    }

    // Add boxes to level for existing levels.
    for (int lev = 1; lev <= finest_level; ++lev) {
        if (box_arrays[lev].size() > 0) {
            time_stepper->local_regrid->AddBoxes(lev, box_arrays[lev]);
            geom[lev] = geom[lev - 1];
            geom[lev].refine(amrex::IntVect(2, 2, 2));
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
                amrex::IntVect ce(i, j, k);
                if (ba.contains(ce * blocking_factor[lev]))
                    continue;

                // Check all boxes.
                for (int b = 0; b < ab.size(); ++b) {
                    Bubble *bub = &(bubbles[ab[b]]);
                    double minD = DBL_MAX;
                    double maxD = 0;
                    const double Dx[2] = {Distance(i * dxb, bub->x, L),
                                          Distance((i + 1) * dxb, bub->x, L)};
                    const double Dy[2] = {Distance(j * dxb, bub->y, L),
                                          Distance((j + 1) * dxb, bub->y, L)};
                    const double Dz[2] = {Distance(k * dxb, bub->z, L),
                                          Distance((k + 1) * dxb, bub->z, L)};

                    for (int ii = 0; ii <= 1; ++ii) {
                        for (int jj = 0; jj <= 1; ++jj) {
                            for (int kk = 0; kk <= 1; ++kk) {
                                double D = std::sqrt(Dx[ii] * Dx[ii] +
                                                     Dy[jj] * Dy[jj] +
                                                     Dz[kk] * Dz[kk]);
                                minD = std::min(D, minD);
                                maxD = std::max(D, maxD);
                            }
                        }
                    }

                    int ind0 = bub->GetPos(minD);
                    if (ind0 == -1)
                        continue;

                    int ind1 = bub->GetPos(maxD);
                    if (ind1 == -1)
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
        amrex::MultiFab &mf = grid_new[lev];

#pragma omp parallel
        for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
            const amrex::Box &bx = mfi.tilebox();
            const auto &fab = mf.array(mfi);

            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        for (int b = 0; b < ab.size(); ++b) {
                            AddBubble(i, j, k, dx[lev], L, fab, bubbles[ab[b]]);
                        }
                    }
                }
            }
        }
    }
}

}; // namespace FirstOrderPhaseTransition
