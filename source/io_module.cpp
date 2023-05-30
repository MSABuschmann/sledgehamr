#include "io_module.h"

namespace sledgehamr {

IOModule::IOModule(Sledgehamr* owner) {
    sim = owner;

    // Determine and create output folder.
    amrex::ParmParse pp("output");
    std::string output_folder;
    pp.get("output_folder", output_folder);
    amrex::UtilCreateDirectory(output_folder.c_str(), 0755);

    // Add various output formats.
    // Slices.
    double interval_slices = -1;
    pp.query("interval_slices", interval_slices);

    if (interval_slices >= 0) {
        OutputModule out(output_folder + "/slices",
                         OUTPUT_FCT(IOModule::WriteSlices),
                         interval_slices);
        output.push_back(out);
    }

    // Full coarse box.
    double interval_coarse_box = -1;
    pp.query("interval_coarse_box", interval_coarse_box);
    // TODO Check that power of 2 and <= blocking_factor.
    pp.query("coarse_box_downsample_factor", coarse_box_downsample_factor);

    if (interval_coarse_box >= 0) {
        OutputModule out(output_folder + "/coarse_box",
                         OUTPUT_FCT(IOModule::WriteCoarseBox),
                         interval_coarse_box);
        output.push_back(out);
    }

    // Entire volume.
    double interval_full_box = -1;
    pp.query("interval_full_box", interval_full_box);
    // TODO Check that power of 2 and <= blocking_factor.
    pp.query("full_box_downsample_factor", full_box_downsample_factor);

    if (interval_full_box >= 0) {
        OutputModule out(output_folder + "/full_box",
                         OUTPUT_FCT(IOModule::WriteFullBox),
                         interval_full_box);
        output.push_back(out);
    }

    // Slices of truncation errors.
    double interval_slices_truncation_error = -1;
    pp.query("interval_slices_truncation_error",
             interval_slices_truncation_error);

    if (interval_slices_truncation_error >= 0) {
        OutputModule out(output_folder + "/slices_truncation_error",
                         OUTPUT_FCT(IOModule::WriteSlicesTruncationError),
                         interval_slices_truncation_error);
        output.push_back(out);
    }

    // Full coarse box of truncation errors.
    double interval_coarse_box_truncation_error = -1;
    pp.query("interval_coarse_box_truncation_error",
             interval_coarse_box_truncation_error);
    // TODO Check that power of 2 and <= blocking_factor.
    pp.query("truncation_error_coarse_box_downsample_factor",
             truncation_error_coarse_box_downsample_factor);

    if (interval_coarse_box_truncation_error >= 0) {
        OutputModule out(output_folder + "/coarse_box_truncation_error",
                         OUTPUT_FCT(IOModule::WriteCoarseBoxTruncationError),
                         interval_coarse_box_truncation_error);
        output.push_back(out);
    }

/*
    // Entire volume of truncation errors.
    double interval_full_box_truncation_error = -1;
    pp.query("interval_full_box_truncation_error",
             interval_full_box_truncation_error);

    if (interval_full_box_truncation_error >= 0) {
        OutputModule out(output_folder + "/full_box_truncation_error",
                         OUTPUT_FCT(IOModule::WriteFullBoxTruncationError),
                         interval_full_box_truncation_error);
        output.push_back(out);
    }

   // yt output
    double interval_yt = -1;
    pp.query("interval_yt", interval_yt);

    if (interval_yt >= 0) {
        OutputModule out(output_folder + "/yt",
                         OUTPUT_FCT(IOModule::WriteYt),
                         interval_yt);
        output.push_back(out);
    }

    // Projections.
    double interval_projections = -1;
    pp.query("interval_projections", interval_projections);

    if (interval_projections >= 0) {
        OutputModule out(output_folder + "/projections",
                         OUTPUT_FCT(IOModule::WriteProjections),
                         interval_projections);
        output.push_back(out);
    }

     // Projections of truncation errors.
    double interval_projections_truncation_error = -1;
    pp.query("interval_projections_truncation_error",
             interval_projections_truncation_error);

    if (interval_projections_truncation_error >= 0) {
        OutputModule out(output_folder + "/projections_truncation_error",
                         OUTPUT_FCT(IOModule::WriteProjectionsTruncationError),
                         interval_projections_truncation_error);
        output.push_back(out);
    }

    // Spectra.
    double interval_spectra = -1;
    pp.query("interval_spectra", interval_spectra);

    if (interval_spectra >= 0) {
        OutputModule out(output_folder + "/spectra",
                         OUTPUT_FCT(IOModule::WriteSpectra),
                         interval_spectra);
        output.push_back(out);
    }

    // GW spectra.
    double interval_gw_spectra = -1;
    pp.query("interval_gw_spectra", interval_gw_spectra);

    if (interval_gw_spectra >= 0) {
        OutputModule out(output_folder + "/gw_spectra",
                         OUTPUT_FCT(IOModule::WriteGwSpectra),
                         interval_gw_spectra);
        output.push_back(out);
    }

    // Checkpoint.
    double interval_checkpoints = -1;
    pp.query("interval_checkpoints", interval_checkpoints);

    if (interval_checkpoints >= 0) {
        OutputModule out(output_folder + "/checkpoints",
                         OUTPUT_FCT(IOModule::WriteCheckpoint),
                         interval_checkpoints);
        output.push_back(out);
    }
*/
}

void IOModule::Write(bool force) {
    for (OutputModule& out : output)
        out.Write(sim->grid_new[0].t, force);
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

void IOModule::WriteSlices(double time, std::string prefix) {
    amrex::Print() << "Write slices: " << prefix << std::endl;
    DoWriteSlices(time, prefix, false);
}

void IOModule::WriteSlicesTruncationError(double time, std::string prefix) {
    amrex::Print() << "Write truncation error slices: " << prefix << std::endl;
    DoWriteSlices(time, prefix, false);
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

void IOModule::WriteCoarseBox(double time, std::string prefix) {
    amrex::Print() << "Write coarse level box: " << prefix << std::endl;
    DoWriteCoarseBox(time, prefix, coarse_box_downsample_factor, false);
}

void IOModule::WriteCoarseBoxTruncationError(double time, std::string prefix) {
    amrex::Print() << "Write truncation errors on coarse level box: "
                   << prefix << std::endl;
    DoWriteCoarseBox(time, prefix,
                     truncation_error_coarse_box_downsample_factor, true);
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

        int dimx = (hx-lx) / downsample_factor / ndist;
        int dimy = (hy-ly) / downsample_factor / ndist;
        int dimz = (hz-lz) / downsample_factor / ndist;
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

                        int ind = (i-lx)/ndist*dimy*dimz + (j-ly)/ndist*dimz
                                + (k-lz)/ndist;

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

void IOModule::WriteFullBox(double time, std::string prefix) {
    amrex::Print() << "Write full box at all levels: " << prefix << std::endl;

    for (int lev = 0; lev <= sim->finest_level; ++lev) {
        // Create folder and file.
        std::string folder = prefix + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(folder.c_str(), 0755);

        std::string filename = folder + "/"
                        + std::to_string(amrex::ParallelDescriptor::MyProc())
                        + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        const LevelData* state = &sim->grid_new[lev];

        // Write output to file.
        WriteLevel(time, state, lev, file_id, "data", 
                   full_box_downsample_factor, false);
        H5Fclose(file_id);
    }
}

}; // namespace sledgehamr
