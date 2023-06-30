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

    if (amrex::FileExists(output_folder) && !sim->restart_sim &&
        amrex::ParallelDescriptor::IOProcessor())  {
        const char* msg = 
                "sledgehamr::IOModule: Output folder already exists!\n"
                "If you intended to restart the simulation from the latest "
                "checkpoint within this folder please add 'sim.restart = 1' "
                "to your input file.\nOtherwise please chose a different "
                "directory.";
        amrex::Abort(msg);
    }

    amrex::ParallelDescriptor::Barrier();
    amrex::UtilCreateDirectory(output_folder.c_str(), 0755);

    // Add various output formats.
    // Slices.
    double interval_slices = -1;
    pp.query("interval_slices", interval_slices);
    idx_slices = output.size();
    output.emplace_back(output_folder + "/slices",
                        OUTPUT_FCT(IOModule::WriteSlices),
                        interval_slices);

    // Full coarse box.
    double interval_coarse_box = -1;
    pp.query("interval_coarse_box", interval_coarse_box);
    pp.query("coarse_box_downsample_factor", coarse_box_downsample_factor);
    CheckDownsampleFactor(coarse_box_downsample_factor,
                         "coarse_box_downsample_factor", 0);
    idx_coarse_box = output.size();
    output.emplace_back(output_folder + "/coarse_box",
                        OUTPUT_FCT(IOModule::WriteCoarseBox),
                        interval_coarse_box);

    // Entire volume.
    double interval_full_box = -1;
    pp.query("interval_full_box", interval_full_box);
    pp.query("full_box_downsample_factor", full_box_downsample_factor);
    CheckDownsampleFactor(full_box_downsample_factor,
                         "full_box_downsample_factor", sim->max_level);
    idx_full_box = output.size();
    output.emplace_back(output_folder + "/full_box",
                        OUTPUT_FCT(IOModule::WriteFullBox),
                        interval_full_box);

    // Slices of truncation errors.
    double interval_slices_truncation_error = -1;
    pp.query("interval_slices_truncation_error",
             interval_slices_truncation_error);
    idx_slices_truncation_error = output.size();
    output.emplace_back(output_folder + "/slices_truncation_error",
                        OUTPUT_FCT(IOModule::WriteSlicesTruncationError),
                        interval_slices_truncation_error);

    // Full coarse box of truncation errors.
    double interval_coarse_box_truncation_error = -1;
    pp.query("interval_coarse_box_truncation_error",
             interval_coarse_box_truncation_error);
    pp.query("coarse_box_truncation_error_downsample_factor",
             coarse_box_truncation_error_downsample_factor);
    CheckDownsampleFactor(coarse_box_truncation_error_downsample_factor,
                         "coarse_box_truncation_error_downsample_factor", 0);
    idx_coarse_box_truncation_error = output.size();
    output.emplace_back(output_folder + "/coarse_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteCoarseBoxTruncationError),
                        interval_coarse_box_truncation_error);

    // Entire volume of truncation errors.
    double interval_full_box_truncation_error = -1;
    pp.query("interval_full_box_truncation_error",
             interval_full_box_truncation_error);
    pp.query("full_box_truncation_error_downsample_factor",
             full_box_truncation_error_downsample_factor);
    CheckDownsampleFactor(full_box_truncation_error_downsample_factor,
                         "full_box_truncation_error_downsample_factor",
                          sim->max_level);
    idx_full_box_truncation_error = output.size();
    output.emplace_back(output_folder + "/full_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteFullBoxTruncationError),
                        interval_full_box_truncation_error);

/*
   // yt output
    double interval_yt = -1;
    pp.query("interval_yt", interval_yt);
    idx_yt = output.size();
    output.emplace_back(output_folder + "/yt",
                        OUTPUT_FCT(IOModule::WriteYt),
                        interval_yt);
*/

    // Projections.
    double interval_projections = -1;
    pp.query("interval_projections", interval_projections);
    idx_projections = output.size();
    output.emplace_back(output_folder + "/projections",
                        OUTPUT_FCT(IOModule::WriteProjections),
                        interval_projections);

    // Spectra.
    double interval_spectra = -1;
    pp.query("interval_spectra", interval_spectra);
    idx_spectra = output.size();
    output.emplace_back(output_folder + "/spectra",
                        OUTPUT_FCT(IOModule::WriteSpectra),
                        interval_spectra);

    // GW spectra.
    double interval_gw_spectra = -1;
    pp.query("interval_gw_spectra", interval_gw_spectra);
    idx_gw_spectra = output.size();
    output.emplace_back(output_folder + "/gw_spectra",
                        OUTPUT_FCT(IOModule::WriteGravitationalWaveSpectrum),
                        interval_gw_spectra);

    // Checkpoint. Always add checkpoints last.
    double interval_checkpoints = -1;
    pp.query("interval_checkpoints", interval_checkpoints);
    idx_checkpoints = output.size();
    output.emplace_back(output_folder + "/checkpoints",
                        OUTPUT_FCT(IOModule::WriteCheckpoint),
                        interval_checkpoints);
}

void IOModule::Write(bool force) {
    // Make sure checkpoints are written last to ensure most up-to-date meta
    // data.
    for(int i = 0; i < output.size(); ++i) {
        if (i != idx_checkpoints)
            output[i].Write(sim->grid_new[0].t, force);
    }
    output[idx_checkpoints].Write(sim->grid_new[0].t, force);
}

void IOModule::FillLevelFromFile (int lev)
{
    // Figure out how files are organized.
    amrex::ParmParse pp("input");

    int N_chunks = 1;
    pp.query("N_chunks", N_chunks);

    if (N_chunks == 1) {
        FillLevelFromFile_NoChunks(lev);
    } else if (N_chunks > 1) {
        FillLevelFromFile_Chunks(lev);
    }
}

void IOModule::FillLevelFromFile_Chunks(int lev) {
    // TODO: Include version that supports files split in chunks.
}

void IOModule::FillLevelFromFile_NoChunks(int lev) {
    amrex::ParmParse pp("input");
    std::string initial_state_file = "";
    pp.query("initial_state", initial_state_file);

    const int ncomp = sim->scalar_fields.size();
    const int dimN = sim->dimN[lev];
    const long long dsetsize = dimN*dimN*dimN;
    double* input_data = new double [dsetsize] ();

    // Iterate over fields but introduce offset such that each node grabs a
    // different file first.
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
            ReadFromHDF5(initial_state_file_component,
                    {sim->scalar_fields[f2]->name, "data"},
                    input_data);
            FillLevelFromArray(lev, f2, input_data, dimN);
        }
    }

    delete[] input_data;
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
        float *output_arr = new float[len]();
        for (int f=0; f<sim->scalar_fields.size(); ++f) {
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
        }

        delete[] output_arr;
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
    amrex::Print() << "Writing checkpoint " << prefix << std::endl;

    const int nlevels = sim->finest_level + 1;
    amrex::PreBuildDirectorHierarchy(prefix, "Level_", nlevels, true);

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(prefix + "/BoxArrays");
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                           std::ofstream::trunc |
                                           std::ofstream::binary);
        if( !HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        for (int lev = 0; lev < nlevels; ++lev) {
            sim->boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }

        std::string filename = prefix + "/Meta.hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                              H5P_DEFAULT);

        const int nparams = 8;
        double header_data[nparams] = {
                time,
                (double)amrex::ParallelDescriptor::NProcs(),
                (double)sim->finest_level,
                (double)sim->dimN[0],
                (double)sim->nghost,
                (double)sim->scalar_fields.size(),
                (double)output.size(),
                (double)idx_checkpoints
        };
        WriteToHDF5(file_id, "Header", header_data, nparams);

        // levels: blocking_factor, istep
        std::vector<int> bf(nlevels), istep(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            bf[lev] = sim->blocking_factor[lev][0];
            istep[lev] = sim->grid_new[lev].istep;
        }
        WriteToHDF5(file_id, "blocking_factors", &(bf[0]), nlevels);
        WriteToHDF5(file_id, "isteps", &(istep[0]), nlevels);

        // outputs: last id, last time written
        const int noutput = output.size();
        std::vector<int> next_id(noutput);
        std::vector<double> last_time_written(noutput);
        for (int i = 0; i < noutput; ++i) {
            next_id[i] = output[i].GetNextId();
            last_time_written[i] = output[i].GetLastTimeWritten();
        }
        WriteToHDF5(file_id, "next_id", &(next_id[0]), noutput);
        WriteToHDF5(file_id, "last_time_written", &(last_time_written[0]),
                    noutput);

        H5Fclose(file_id);
    }

    // Write the MultiFab data.
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::VisMF::Write(
                sim->grid_new[lev],
                amrex::MultiFabFileFullPrefix(
                        lev, prefix, "Level_", "Cell"));
    }

    return true;
}

void IOModule::RestartSim() {
    int latest = FindLatestCheckpoint();

    if (latest == -1 && amrex::ParallelDescriptor::IOProcessor())
        amrex::Abort("Sledgehamr::IOModule::RestartSim: No checkpoint found!");
    amrex::ParallelDescriptor::Barrier(); 

    ReadCheckpoint(latest);
}

int IOModule::FindLatestCheckpoint() {
    int latest = 0;
    std::string folder = output_folder + "/checkpoints/"
                       + std::to_string(latest);

    while (amrex::FileExists(folder)) {
        ++latest;
        folder = output_folder + "/checkpoints/" + std::to_string(latest);
    }
    return latest - 1;
}

void IOModule::ReadCheckpoint(int id) {
    std::string folder = output_folder + "/checkpoints/" + std::to_string(id);
    amrex::Print() << "Restarting from checkpoint: " << folder << std::endl;

    const int nparams = 8;
    double header[nparams];
    std::string filename = folder + "/Meta.hdf5";
    ReadFromHDF5(filename, {"Header"}, header);

    double time = header[0];
    int MPIranks = static_cast<int>(header[1]);
    int finest_level = static_cast<int>(header[2]);
    int dim0 = static_cast<int>(header[3]);
    int nghost = static_cast<int>(header[4]);
    int nscalars = static_cast<int>(header[5]);
    int noutput = static_cast<int>(header[6]);
    int npredefoutput = static_cast<int>(header[7]);

    if (nscalars != sim->scalar_fields.size()) {
        const char* msg = "Sledgehamr::IOModule::ReadCheckpoint: "
                          "Number of scalar fields has changed!";
        amrex::Abort(msg);
    }

    if (noutput != output.size() || npredefoutput != idx_checkpoints) {
        const char* msg = "Sledgehamr::IOModule::ReadCheckpoint: "
                          "Number of output types changed!";
        amrex::Abort(msg);
    }

    std::string File(folder + "/BoxArrays");
    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());
    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    sim->finest_level = finest_level;
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        amrex::DistributionMapping dm {ba, amrex::ParallelDescriptor::NProcs()};
        sim->SetBoxArray(lev, ba);
        sim->SetDistributionMap(lev, dm);

        // In case nghost changed, we can create grid_old already with the new
        // value. grid_new has to be set to the old value first since ghost
        // cells are saved in the checkpoint. We change nghost for this MultiFab
        // later below.
        sim->grid_old[lev].define(ba, dm, nscalars, sim->nghost, time);
        sim->grid_new[lev].define(ba, dm, nscalars, nghost, time);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
                sim->grid_new[lev],
                amrex::MultiFabFileFullPrefix(lev, folder, "Level_", "Cell"));
    }

    if (nghost != sim->nghost) {
        amrex::Print() << "#warning: Number of ghost cells has changed!\n"
                       << "checkpoint: " << nghost
                       << " vs input file: " << sim->nghost << std::endl;
        // TODO
    }

    if (MPIranks != amrex::ParallelDescriptor::NProcs()) {
        amrex::Print() << "#warning: Number of MPI ranks has changed. Will "
                       << "regrid coarse level to satisfy new constraint."
                       << std::endl;
        // TODO
    }
}

void IOModule::GotoNextLine(std::istream& is) {
    constexpr std::streamsize bl_ignore_max {100000};
    is.ignore(bl_ignore_max, '\n');
}

}; // namespace sledgehamr
