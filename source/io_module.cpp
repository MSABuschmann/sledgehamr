#include <filesystem>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>

#include "io_module.h"
#include "sledgehamr_utils.h"
#include "output_types/slices.h"
#include "output_types/level_writer.h"

namespace sledgehamr {

IOModule::IOModule(Sledgehamr* owner) : sim(owner) {
    ParseParams();
    CheckIfOutputAlreadyExists(output_folder);
    CheckIfOutputAlreadyExists(alternative_output_folder);
    amrex::ParallelDescriptor::Barrier();
    CreateOutputFolder(output_folder);
    CreateOutputFolder(alternative_output_folder);
    AddOutputModules();
}

void IOModule::ParseParams() {
    amrex::ParmParse pp("output");
    pp.get("output_folder", output_folder);
    pp.query("alternative_output_folder", alternative_output_folder);
    pp.query("checkpoints.rolling", rolling_checkpoints);

    amrex::ParmParse pp_sim("sim");
    pp_sim.query("delete_restart_checkpoint", delete_restart_checkpoint);
}

void IOModule::AddOutputModules() {
    idx_slices = output.size();
    output.emplace_back("slices", OUTPUT_FCT(IOModule::WriteSlices));

    idx_coarse_box = output.size();
    output.emplace_back("coarse_box", OUTPUT_FCT(IOModule::WriteCoarseBox));

    idx_full_box = output.size();
    output.emplace_back("full_box", OUTPUT_FCT(IOModule::WriteFullBox));

    idx_slices_truncation_error = output.size();
    output.emplace_back("slices_truncation_error",
                        OUTPUT_FCT(IOModule::WriteSlicesTruncationError));

    idx_coarse_box_truncation_error = output.size();
    output.emplace_back("coarse_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteCoarseBoxTruncationError));

    idx_full_box_truncation_error = output.size();
    output.emplace_back("full_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteFullBoxTruncationError));

    idx_projections = output.size();
    output.emplace_back("projections", OUTPUT_FCT(IOModule::WriteProjections));

    idx_spectra = output.size();
    output.emplace_back("spectra", OUTPUT_FCT(IOModule::WriteSpectra));

    idx_gw_spectra = output.size();
    output.emplace_back("gw_spectra",
                        OUTPUT_FCT(IOModule::WriteGravitationalWaveSpectrum));

    idx_performance_monitor = output.size();
    output.emplace_back("performance_monitor",
                        OUTPUT_FCT(IOModule::WritePerformanceMonitor));

    // Checkpoint. Always add checkpoints last.
    idx_checkpoints = output.size();
    output.emplace_back("checkpoints", OUTPUT_FCT(IOModule::WriteCheckpoint));

    bool write_at_start = false;
    amrex::ParmParse pp("output");
    pp.query("write_at_start", write_at_start);
    if (!write_at_start) {
        for (OutputModule& out : output) {
            out.SetLastTimeWritten( sim->t_start );
        }
    }
}

void IOModule::CheckIfOutputAlreadyExists(std::string folder) {
    if (folder=="")
        return;

    // Determine and create output folder.
    bool rename_old = false;
    amrex::ParmParse pp("output");
    pp.query("rename_old_output", rename_old);

    if (amrex::FileExists(folder) && !sim->restart_sim &&
        amrex::ParallelDescriptor::IOProcessor() && !rename_old)  {
        std::string msg =
                "sledgehamr::IOModule: Output folder " + folder + " "
                "already exists!\n"
                "If you intended to restart the simulation from the latest "
                "checkpoint within this folder please add 'sim.restart = 1' "
                "to your input file.\nOtherwise please chose a different "
                "directory.";
        amrex::Abort(msg);
    }
}

void IOModule::CreateOutputFolder(std::string folder) {
    if (folder=="")
        return;

    if (!sim->restart_sim) {
        std::string tmp = folder;
        while(tmp.back() == '/') {
            tmp.pop_back();
        }

        amrex::Print() << "Create output folder: " << tmp << std::endl;
        amrex::UtilCreateCleanDirectory(tmp, true);
    }
}

void IOModule::Write(bool force) {
    // Make sure checkpoints are written last to ensure most up-to-date meta
    // data.
    for(int i = 0; i < output.size(); ++i) {
        if (i != idx_checkpoints) {
            sim->performance_monitor->Start(
                    sim->performance_monitor->idx_output, i);
            output[i].Write(sim->grid_new[0].t, force);
            sim->performance_monitor->Stop(
                    sim->performance_monitor->idx_output, i);
        }
    }

    sim->performance_monitor->Start(
            sim->performance_monitor->idx_output, idx_checkpoints);
    output[idx_checkpoints].Write(sim->grid_new[0].t, force);
    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_output, idx_checkpoints);
}

void IOModule::FillLevelFromFile(int lev) {
    amrex::ParmParse pp("input");
    std::string initial_state_file = "";
    pp.query("initial_state", initial_state_file);

    if (amrex::FileExists(initial_state_file + "/Meta.hdf5")) {
        FillLevelFromCheckpointFile(lev, initial_state_file);
    } else {
        FillLevelFromHdf5File(lev, initial_state_file);
    }
}

void IOModule::FillLevelFromCheckpointFile(int lev, std::string folder) {
    Checkpoint chk(sim, folder);
    chk.Read();

    if (delete_restart_checkpoint)  {
        old_checkpoint = folder;
    }
}

void IOModule::FillLevelFromHdf5File(int lev, std::string initial_state_file) {
    amrex::ParmParse pp("input");

    // Iterate over fields but introduce offset such that each node grabs a
    // different file first.
    const int ncomp = sim->scalar_fields.size();
    int lr = amrex::ParallelDescriptor::MyProc();
    int mr = amrex::ParallelDescriptor::NProcs();
    for (int f=0; f<ncomp; ++f) {
        int f2 = (f+lr) % ncomp;

        // Determine which input file to use if any.
        std::string initial_state_file_component = "";
        std::string query = "initial_state_" + sim->scalar_fields[f2]->name;
        pp.query(query.c_str(), initial_state_file_component);

        if (initial_state_file_component == "")
            initial_state_file_component = initial_state_file;

        // If no file found, fill level with 0's. Otherwise read file and fill
        // LevelData.
        if(initial_state_file_component == "") {
            FillLevelFromConst(lev, f2, 0);
        } else {
            amrex::Print() << "Reading initial state for "
                           << sim->scalar_fields[f2]->name
                           << " from " << initial_state_file_component
                           << std::endl;

            // Test if chunks exist.
            std::string chunk1 = sim->scalar_fields[f2]->name + "_"
                               + std::to_string(lr);
            std::string chunk2 = "data_" + std::to_string(lr);
            std::string existing_chunk = ExistingDataset(
                    initial_state_file_component, {chunk1, chunk2});

            if (existing_chunk == "") {
                const int dimN = sim->dimN[lev];
                const long long dsetsize = dimN*dimN*dimN;
                double* input_data = new double [dsetsize] ();

                if (!ReadFromHDF5(initial_state_file_component,
                        {sim->scalar_fields[f2]->name, "data"},
                        input_data)) {
                    std::string msg =
                               "Sledgehamr::IOModule::FillLevelFromHdf5File: "
                               "Could not find initial state data for "
                             + sim->scalar_fields[f2]->name + "!";
                    amrex::Abort(msg);
                }

                FillLevelFromArray(lev, f2, input_data, dimN);
                delete[] input_data;
            } else {
                amrex::Box bx = sim->grid_new[lev].boxArray()[lr];
                const long long dsetsize= bx.numPts();
                double* input_data = new double [dsetsize] ();

               if (!ReadFromHDF5(initial_state_file_component,
                        {chunk1, chunk2},
                        input_data)) {
                    std::string msg =
                               "Sledgehamr::IOModule::FillLevelFromHdf5File: "
                               "Could not find initial state chunk "
                             + chunk1 + "!";
                    amrex::Abort(msg);
                }

                FillChunkFromArray(lev, f2, input_data);
                delete[] input_data;
            }
        }
    }
}

void IOModule::FillLevelFromArray(int lev, const int comp, double* data,
                                  const long long dimN) {
    LevelData& state = sim->grid_new[lev];

#pragma omp parallel
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = state.array(mfi);

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    long long ind =  static_cast<long long>(i) * dimN*dimN
                                   + static_cast<long long>(j) * dimN
                                   + static_cast<long long>(k);
                    state_arr(i,j,k,comp) = data[ind];
                }
            }
        }
    }
}

void IOModule::FillChunkFromArray(int lev, const int comp, double* data) {
    LevelData& state = sim->grid_new[lev];

#pragma omp parallel
    for (amrex::MFIter mfi(state, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = state.array(mfi);

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);
        const int lx = bx.length(0);
        const int ly = bx.length(1);
        const int lz = bx.length(2);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    long long ind =  static_cast<long long>(i-lo.x) * lz*ly
                                   + static_cast<long long>(j-lo.y) * lz
                                   + static_cast<long long>(k-lo.z);

                    state_arr(i,j,k,comp) = data[ind];
                }
            }
        }
    }
}

void IOModule::FillLevelFromConst(int lev, const int comp, const double c) {
    LevelData& state = sim->grid_new[lev];

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid();
            ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = state.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            noexcept {
            state_arr(i,j,k,comp) = c;
        });
    }
}

bool IOModule::WriteSlices(double time, std::string prefix) {
    Slices slices(sim, prefix, false);
    slices.Write();
    return true;
}

bool IOModule::WriteSlicesTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    Slices slices(sim, prefix, true);
    slices.Write();
    return true;
}

bool IOModule::WriteCoarseBox(double time, std::string prefix) {
    LevelWriter writer(sim, prefix, idx_coarse_box);
    writer.Write();
    return true;
}

bool IOModule::WriteCoarseBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_coarse_box_truncation_error);
    writer.Write();
    return true;
}

bool IOModule::WriteFullBox(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_full_box);
    writer.Write();
    return true;
}

bool IOModule::WriteFullBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_full_box_truncation_error);
    writer.Write();
    return true;
}

bool IOModule::WriteProjections(double time, std::string prefix) {
    if (projections.empty())
        return false;

    amrex::Print() << "Write projections: " << prefix << std::endl;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/projections.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    for (int p = 0; p < projections.size(); ++p)
        projections[p].Compute(p, file_id, sim);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

bool IOModule::WriteSpectra(double time, std::string prefix) {
    if (spectra.empty())
        return false;

    amrex::Print() << "Compute spectra: " << prefix << std::endl;

    sim->ReadSpectrumKs();

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    for (int p = 0; p < spectra.size(); ++p)
        spectra[p].Compute(p, file_id, sim);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

bool IOModule::WriteGravitationalWaveSpectrum(double time, std::string prefix) {
    if (!sim->with_gravitational_waves)
        return false;

    amrex::Print() << "Compute gravitational wave spectrum: " << prefix
                   << std::endl;

    sim->ReadSpectrumKs();

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    sim->gravitational_waves->ComputeSpectrum(file_id);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

bool IOModule::WriteCheckpoint(double time, std::string prefix) {
    Checkpoint chk(sim, prefix);
    chk.Write();

    if (rolling_checkpoints) {
        if (old_checkpoint != "") {
            Checkpoint chk_del(sim, old_checkpoint);
            chk_del.Delete();
        }

        old_checkpoint = prefix;
    }

    return true;
}

void IOModule::RestartSim() {
    amrex::ParmParse pp("sim");
    std::string selected_chk = "None Selected";
    pp.query("select_checkpoint", selected_chk);

    if (selected_chk == "None Selected") {
        int latest     = FindLatestCheckpoint(output_folder);
        int latest_alt = FindLatestCheckpoint(alternative_output_folder);

        if (latest > latest_alt) {
            initial_chk = output_folder + "/checkpoints/"
                        + std::to_string(latest);
        } else {
            initial_chk = alternative_output_folder + "/checkpoints/"
                        + std::to_string(latest_alt);
        }

        if (std::max(latest, latest_alt) == -1 &&
            amrex::ParallelDescriptor::IOProcessor()) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "No checkpoint found!";
            amrex::Abort(msg);
        }
   } else {
        if (amrex::is_integer(selected_chk.c_str())) {

            initial_chk = output_folder + "/checkpoints/" + selected_chk;

            if (!amrex::FileExists(initial_chk)) {
                initial_chk = alternative_output_folder + "/checkpoints/"
                            + selected_chk;
            }
        } else {
            initial_chk = selected_chk;
        }

        if (!amrex::FileExists(initial_chk)) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "Selected checkpoint not found!";
            amrex::Abort(msg);
        }
    }

    amrex::ParallelDescriptor::Barrier();
    Checkpoint chk(sim, initial_chk);
    chk.Read();

    if (delete_restart_checkpoint) {
        old_checkpoint = initial_chk;
    }
}

void IOModule::UpdateOutputModules() {
    Checkpoint chk(sim, initial_chk);
    chk.UpdateOutputModules();
}

int IOModule::FindLatestCheckpoint(std::string folder) {
    if (folder == "")
        return -1;

    std::string prefix = folder + "/checkpoints/";
    std::vector<std::string> folders = GetDirectories(prefix);

    double latest_time = -DBL_MAX;
    int latest_chk = -1;
    for (std::string& str : folders) {
        str.erase(0, prefix.size());

        if (amrex::is_integer(str.c_str())) {
            std::string chk_folder = prefix + str;

            Checkpoint chk(sim, chk_folder);
            if (!chk.ReadHeader())
                continue;

            double time = chk.GetTime();
            if (time > latest_time) {
                latest_time = time;
                latest_chk = std::stoi(str);
            }
        }
    }

    return latest_chk;
}

std::vector<std::string> IOModule::GetDirectories(const std::string prefix) {
    std::vector<std::string> res;
    for (auto& p : std::filesystem::recursive_directory_iterator(prefix))
        if (p.is_directory())
            res.push_back(p.path().string());
    return res;
}

bool IOModule::WritePerformanceMonitor(double time, std::string prefix) {
    if (!sim->performance_monitor->IsActive())
        return false;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/log.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    sim->performance_monitor->Log(file_id);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        H5Fclose(file_id);
    }

    return true;
}

void IOModule::WriteBoxArray(amrex::BoxArray& ba) {
    if (!amrex::ParallelDescriptor::IOProcessor()) 
        return;

    const int nba = ba.size();
    std::vector<int> x0(nba), y0(nba), z0(nba), x1(nba), y1(nba), z1(nba);
    for (int b = 0; b < nba; ++b) {
        x0[b] = ba[b].smallEnd(0);
        y0[b] = ba[b].smallEnd(1);
        z0[b] = ba[b].smallEnd(2);
        x1[b] = ba[b].bigEnd(0);
        y1[b] = ba[b].bigEnd(1);
        z1[b] = ba[b].bigEnd(2);
    }
    int header[1] = {nba}; 
 
    std::string filename = output_folder + "/box_layout.h5";
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    WriteToHDF5(file_id, "header", header, 1);
    WriteToHDF5(file_id, "x0", x0.data(), nba);
    WriteToHDF5(file_id, "y0", y0.data(), nba);
    WriteToHDF5(file_id, "z0", z0.data(), nba);
    WriteToHDF5(file_id, "x1", x1.data(), nba);
    WriteToHDF5(file_id, "y1", y1.data(), nba);
    WriteToHDF5(file_id, "z1", z1.data(), nba);
    H5Fclose(file_id);
}

std::string IOModule::ExistingDataset(std::string filename,
                                            std::vector<std::string> dnames) {
    // Try and open HDF5 file.
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == H5I_INVALID_HID) {
        amrex::Abort("#error: Could not open file: " + filename);
    }

    // Try and find dataset. Iterate over vector, use first to be found.
    std::string dname_found = "";
    for (std::string dname : dnames) {
        htri_t exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);

        if (exists > 0)
            return dname;
    }

    return "";
}

}; // namespace sledgehamr
