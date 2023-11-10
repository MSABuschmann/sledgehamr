#include <filesystem>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>

#include "io_module.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

IOModule::IOModule(Sledgehamr* owner) {
    sim = owner;

    // Determine and create output folder.
    amrex::ParmParse pp("output");
    pp.get("output_folder", output_folder);
    bool rename_old = false;
    pp.query("rename_old_output", rename_old);
    pp.query("rolling_checkpoints", rolling_checkpoints);

    amrex::ParmParse pp_sim("sim");
    pp_sim.query("delete_restart_checkpoint", delete_restart_checkpoint);

    if (amrex::FileExists(output_folder) && !sim->restart_sim &&
        amrex::ParallelDescriptor::IOProcessor() && !rename_old)  {
        const char* msg =
                "sledgehamr::IOModule: Output folder already exists!\n"
                "If you intended to restart the simulation from the latest "
                "checkpoint within this folder please add 'sim.restart = 1' "
                "to your input file.\nOtherwise please chose a different "
                "directory.";
        amrex::Abort(msg);
    }

    amrex::ParallelDescriptor::Barrier();

    if (!sim->restart_sim) {
        std::string tmp = output_folder;
        while(tmp.back() == '/') {
            tmp.pop_back();
        }

        amrex::Print() << "Create output folder: " << tmp << std::endl;
        amrex::UtilCreateCleanDirectory(tmp, true);
    }

    // Add various output formats.
    idx_slices = output.size();
    output.emplace_back(output_folder, "slices",
                        OUTPUT_FCT(IOModule::WriteSlices));

    // Full coarse box.
    idx_coarse_box = output.size();
    output.emplace_back(output_folder, "coarse_box",
                        OUTPUT_FCT(IOModule::WriteCoarseBox));

    // Entire volume.
    idx_full_box = output.size();
    output.emplace_back(output_folder, "full_box",
                        OUTPUT_FCT(IOModule::WriteFullBox));

    // Slices of truncation errors.
    idx_slices_truncation_error = output.size();
    output.emplace_back(output_folder, "slices_truncation_error",
                        OUTPUT_FCT(IOModule::WriteSlicesTruncationError));

    // Full coarse box of truncation errors.
    idx_coarse_box_truncation_error = output.size();
    output.emplace_back(output_folder, "coarse_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteCoarseBoxTruncationError));

    // Entire volume of truncation errors.
    idx_full_box_truncation_error = output.size();
    output.emplace_back(output_folder, "full_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteFullBoxTruncationError));

    // Projections.
    idx_projections = output.size();
    output.emplace_back(output_folder, "projections",
                        OUTPUT_FCT(IOModule::WriteProjections));

    // Spectra.
    double interval_spectra = -1;
    pp.query("interval_spectra", interval_spectra);
    idx_spectra = output.size();
    output.emplace_back(output_folder, "spectra",
                        OUTPUT_FCT(IOModule::WriteSpectra),
                        interval_spectra);

    // GW spectra.
    double interval_gw_spectra = -1;
    pp.query("interval_gw_spectra", interval_gw_spectra);
    idx_gw_spectra = output.size();
    output.emplace_back(output_folder, "gw_spectra",
                        OUTPUT_FCT(IOModule::WriteGravitationalWaveSpectrum),
                        interval_gw_spectra);

    // GW spectra.
    double interval_performance_monitor = -1;
    pp.query("interval_performance_monitor", interval_performance_monitor);
    idx_performance_monitor = output.size();
    output.emplace_back(output_folder, "performance_log",
                        OUTPUT_FCT(IOModule::WritePerformanceMonitor),
                        interval_performance_monitor);

    // Checkpoint. Always add checkpoints last.
    double interval_checkpoints = -1;
    pp.query("interval_checkpoints", interval_checkpoints);
    idx_checkpoints = output.size();
    output.emplace_back(output_folder, "checkpoints",
                        OUTPUT_FCT(IOModule::WriteCheckpoint),
                        interval_checkpoints);

    bool write_at_start = false;
    pp.query("write_at_start", write_at_start);
    if (!write_at_start) {
        for (OutputModule& out : output) {
            out.SetLastTimeWritten( sim->t_start );
        }
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
    amrex::Print() << "Write slices: " << prefix << std::endl;
    DoWriteSlices(time, prefix, false);
    return true;
}

bool IOModule::WriteSlicesTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    amrex::Print() << "Write truncation error slices: " << prefix << std::endl;
    DoWriteSlices(time, prefix, true);
    return true;
}

void IOModule::DoWriteSlices(double time, std::string prefix,
                             bool with_truncation_errors) {
    for (int lev = 0; lev <= sim->finest_level; ++lev) {
        // Create folder and file.
        std::string folder = prefix + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(folder.c_str(), 0755);

        std::string filename = folder + "/"
                        + std::to_string(amrex::ParallelDescriptor::MyProc())
                        + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        // Write field data.
        const LevelData* state = &sim->grid_new[lev];
        WriteSingleSlice(time, state, lev, file_id, "x", 0, 1, 2, false);
        WriteSingleSlice(time, state, lev, file_id, "y", 1, 0, 2, false);
        WriteSingleSlice(time, state, lev, file_id, "z", 2, 0, 1, false);

        if (with_truncation_errors) {
            // Write truncation errors.
            const LevelData* state = &sim->grid_old[lev];
            WriteSingleSlice(time, state, lev, file_id, "te_x", 0, 1, 2, true);
            WriteSingleSlice(time, state, lev, file_id, "te_y", 1, 0, 2, true);
            WriteSingleSlice(time, state, lev, file_id, "te_z", 2, 0, 1, true);
        }

        H5Fclose(file_id);
    }
}

void IOModule::WriteSingleSlice(double time, const LevelData* state, int lev,
                                hid_t file_id, std::string ident, int d1,
                                int d2, int d3, bool is_truncation_errors) {
    std::vector<int> le1, he1, le2, he2;
    const int ndist = is_truncation_errors ? 2 : 1;

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();

        // Find box that cuts through selected plane.
        if (bx.smallEnd(d1) == 0) {
            const auto& state_arr = state->array(mfi);

            // Get box dimensions
            int l1 = bx.smallEnd(d2);
            int l2 = bx.smallEnd(d3);
            int h1 = bx.bigEnd(d2)+1;
            int h2 = bx.bigEnd(d3)+1;
            le1.push_back(l1);
            le2.push_back(l2);
            he1.push_back(h1);
            he2.push_back(h2);

            int dim1 = (h1-l1)/ndist;
            int dim2 = (h2-l2)/ndist;

            // Copy data into flattened array for each scalar field.
            long len = dim1 * dim2;
            // TODO Adjust output type.
            float *output_arr = new float[len]();
            for (int f=0; f<sim->scalar_fields.size(); ++f) {
                for (int j=l2; j<h2; ++j) {
                    for (int i=l1; i<h1; ++i) {
                        if (is_truncation_errors && (i%2 != 0 || j%2 != 0))
                            continue;

                        int ind = (i-l1)/ndist*dim2 + (j-l2)/ndist;

                        if (d1 == 0) {
                            output_arr[ind] = state_arr(0,i,j,f);
                        } else if (d1 == 1) {
                            output_arr[ind] = state_arr(i,0,j,f);
                        } else if (d1 == 2) {
                            output_arr[ind] = state_arr(i,j,0,f);
                        }
                    }
                }

                std::string dset_name = sim->scalar_fields[f]->name
                                      + "_" + ident
                                      + "_" + std::to_string(le1.size());
                WriteToHDF5(file_id, dset_name, output_arr, len);
            }

            delete[] output_arr;
        }
    }

    // Write header information for this slice.
    const int nparams = 5;
    double header_data[nparams] = {time,
                (double)amrex::ParallelDescriptor::NProcs(),
                (double)(sim->finest_level),
                (double)sim->dimN[lev],
                (double)le1.size()};
    WriteToHDF5(file_id, "Header_"+ident, header_data, nparams);

    // Write box dimensions so we can reassemble slice.
    if (le1.size() > 0) {
        WriteToHDF5(file_id, "le1_"+ident, (int*)&(le1[0]), le1.size());
        WriteToHDF5(file_id, "le2_"+ident, (int*)&(le2[0]), le2.size());
        WriteToHDF5(file_id, "he1_"+ident, (int*)&(he1[0]), he1.size());
        WriteToHDF5(file_id, "he2_"+ident, (int*)&(he2[0]), he2.size());
    }
}

bool IOModule::WriteCoarseBox(double time, std::string prefix) {
    amrex::Print() << "Write coarse level box: " << prefix << std::endl;
    DoWriteCoarseBox(time, prefix, coarse_box_downsample_factor, false);
    return true;
}

bool IOModule::WriteCoarseBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    amrex::Print() << "Write truncation errors on coarse level box: "
                   << prefix << std::endl;
    DoWriteCoarseBox(time, prefix,
                     coarse_box_truncation_error_downsample_factor, true);
    return true;
}

void IOModule::DoWriteCoarseBox(double time, std::string prefix,
                                int downsample_factor,
                                bool with_truncation_errors) {
    const int lev = 0;

    // Create folder and file.
    std::string filename = prefix + "/"
                    + std::to_string(amrex::ParallelDescriptor::MyProc())
                    + ".hdf5";
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                              H5P_DEFAULT);

    // Write fields.
    const LevelData* state = &sim->grid_new[lev];
    WriteLevel(time, state, lev, file_id, "data", downsample_factor, false);

    if (with_truncation_errors) {
        const LevelData* state = &sim->grid_old[lev];
        WriteLevel(time, state, lev, file_id, "te", downsample_factor, true);
    }

    H5Fclose(file_id);
}

void IOModule::WriteLevel(double time, const LevelData* state, int lev,
                          hid_t file_id, std::string ident, 
                          int downsample_factor, bool is_truncation_errors) {
    const int ndist = is_truncation_errors ? 2 : 1;
    const int grid_density = downsample_factor * ndist;

    std::vector<int> lex, hex, ley, hey, lez, hez;

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = state->array(mfi);

        // Get box dimensions
        int lx = bx.smallEnd(0);
        int ly = bx.smallEnd(1);
        int lz = bx.smallEnd(2);
        int hx = bx.bigEnd(0)+1;
        int hy = bx.bigEnd(1)+1;
        int hz = bx.bigEnd(2)+1;
        lex.push_back(lx);
        ley.push_back(ly);
        lez.push_back(lz);
        hex.push_back(hx);
        hey.push_back(hy);
        hez.push_back(hz);


        int dimx = (hx-lx) / grid_density;
        int dimy = (hy-ly) / grid_density;
        int dimz = (hz-lz) / grid_density;
        double volfac = 1./std::pow(downsample_factor, 3);

        // Copy data into flattened array for each scalar field.
        long len = dimx * dimy * dimz;

        // TODO Adjust output type.
        for (int f=0; f<sim->scalar_fields.size(); ++f) {
            float *output_arr = new float[len]();
            for (int k=lz; k<hz; ++k) {
                for (int j=ly; j<hy; ++j) {
                    for (int i=lx; i<hx; ++i) {
                        if (is_truncation_errors
                            && (i%2 != 0 || j%2 != 0 || k%2 != 0))
                            continue;

                        int ind = (i-lx)/grid_density*dimy*dimz
                                + (j-ly)/grid_density*dimz
                                + (k-lz)/grid_density;

                        if (is_truncation_errors)
                            output_arr[ind] = std::max(
                                    static_cast<float>(state_arr(i,j,k,f)),
                                    output_arr[ind]);
                        else
                            output_arr[ind] += state_arr(i,j,k,f) * volfac;
                    }
                }
            }

            std::string dset_name = sim->scalar_fields[f]->name + "_" + ident
                                  + "_" + std::to_string(lex.size());
            WriteToHDF5(file_id, dset_name, output_arr, len);
            delete[] output_arr;
        }

    }

    // Write header information for this slice.
    const int nparams = 6;
    double header_data[nparams] = {time,
                (double)amrex::ParallelDescriptor::NProcs(),
                (double)sim->finest_level,
                (double)sim->dimN[lev],
                (double)downsample_factor,
                (double)lex.size()};
    WriteToHDF5(file_id, "Header_" + ident, header_data, nparams);

    // Write box dimensions so we can reassemble slice.
    if (lex.size() > 0) {
        WriteToHDF5(file_id, "lex_"+ident, (int*)&(lex[0]), lex.size());
        WriteToHDF5(file_id, "ley_"+ident, (int*)&(ley[0]), ley.size());
        WriteToHDF5(file_id, "lez_"+ident, (int*)&(lez[0]), lez.size());
        WriteToHDF5(file_id, "hex_"+ident, (int*)&(hex[0]), hex.size());
        WriteToHDF5(file_id, "hey_"+ident, (int*)&(hey[0]), hey.size());
        WriteToHDF5(file_id, "hez_"+ident, (int*)&(hez[0]), hez.size());
    }
}

bool IOModule::WriteFullBox(double time, std::string prefix) {
    amrex::Print() << "Write full box at all levels: " << prefix << std::endl;
    DoWriteFullBox(time, prefix, full_box_downsample_factor, false);
    return true;
}

bool IOModule::WriteFullBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    amrex::Print() << "Write truncation errors at all levels: "
                   << prefix << std::endl;
    DoWriteFullBox(time, prefix, full_box_truncation_error_downsample_factor,
                   true);
    return true;
}

void IOModule::DoWriteFullBox(double time, std::string prefix,
                              int downsample_factor,
                              bool is_truncation_errors) {
    for (int lev = 0; lev <= sim->finest_level; ++lev) {
        // Create folder and file.
        std::string folder = prefix + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(folder.c_str(), 0755);

        std::string filename = folder + "/"
                        + std::to_string(amrex::ParallelDescriptor::MyProc())
                        + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        // Write fields.
        const LevelData* state = &sim->grid_new[lev];
        WriteLevel(time, state, lev, file_id, "data", 
                   downsample_factor, false);

        if (is_truncation_errors) {
            // Write truncation errors.
            const LevelData* state = &sim->grid_old[lev];

            WriteLevel(time, state, lev, file_id, "te", 
                       downsample_factor, true);
        }

        H5Fclose(file_id);
    }
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

void IOModule::CheckDownsampleFactor(int factor, std::string name,
                                     int max_level) {
    if (!amrex::ParallelDescriptor::IOProcessor())
        return;

    if (!utils::IsPowerOfTwo(factor)) {
        std::string msg = "sledgehamr::IOModule: Downsample factor output."
                        + name + " is not a power of 2";
        amrex::Abort(msg);
    }

    for (int lev = 0; lev <= max_level; ++lev) {
        if (factor > sim->blocking_factor[lev][0]) {
            std::string msg = "sledgehamr::IOModule: Downsample factor output."
                            + name + " exceeds blocking factor";
            amrex::Abort(msg);
        }
    }
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
    int selected_chk = -1;
    pp.query("select_checkpoint", selected_chk);

    chk_id = 0;
    if (selected_chk < 0) {
        int latest = FindLatestCheckpoint();

        if (latest == -1 && amrex::ParallelDescriptor::IOProcessor()) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "No checkpoint found!";
            amrex::Abort(msg);
        }

        chk_id = latest;
   } else {
        std::string folder = output_folder + "/checkpoints/"
                           + std::to_string(selected_chk);
        if (!amrex::FileExists(folder)) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "Selected checkpoint not found!";
            amrex::Abort(msg);
        }

        chk_id = selected_chk;
    }

    amrex::ParallelDescriptor::Barrier();
    Checkpoint chk(sim, output_folder, chk_id);
    chk.Read();

    if (delete_restart_checkpoint) {
        old_checkpoint = output_folder + "/checkpoints/"
                       + std::to_string(chk_id);
    }
}

void IOModule::UpdateOutputModules() {
    Checkpoint chk(sim, output_folder, chk_id);
    chk.UpdateOutputModules();
}

int IOModule::FindLatestCheckpoint() {
    std::string prefix = output_folder + "/checkpoints/";
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
