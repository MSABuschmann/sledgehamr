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
    double interval_slices = -1;
    pp.query("interval_slices", interval_slices);

    if (interval_slices >= 0) {
        OutputModule out1(output_folder + "/slices",
                          OUTPUT_FCT(IOModule::WriteSlices),
                          interval_slices);
        output.push_back(out1);
    }
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

    for (int lev = sim->shadow_hierarchy; lev <= sim->finest_level; ++lev) {
        // Create folder and file.
        std::string folder = prefix + "/Level_"
                            + std::to_string( lev - sim->shadow_hierarchy );
        amrex::UtilCreateDirectory(folder.c_str(), 0755);

        std::string filename = folder + "/"
                        + std::to_string(amrex::ParallelDescriptor::MyProc())
                        + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        // TODO: Allow for writing of truncation errors .
        const LevelData* state = &sim->grid_new[lev];

        // Write output to file.
        WriteSingleSlice(time, state, lev, file_id, "x", 0, 1, 2);
        WriteSingleSlice(time, state, lev, file_id, "y", 1, 0, 2);
        WriteSingleSlice(time, state, lev, file_id, "z", 2, 0, 1);

        H5Fclose(file_id);
    }
}

void IOModule::WriteSingleSlice(double time, const LevelData* state, int lev,
                                hid_t file_id, std::string ident, int d1,
                                int d2, int d3) {
    std::vector<int> le1, he1, le2, he2;

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

            int dim1 = h1-l1;
            int dim2 = h2-l2;

            // Copy data into flattened array for each scalar field.
            long Len = dim1 * dim2;
            float *output_arr = new float[Len]();
            for (int f=0; f<sim->scalar_fields.size(); ++f) {
                for (int i=l1; i<h1; ++i) {
                    for (int j=l2; j<h2; ++j) {
                        int ind = (i-l1)*dim2 + (j-l2);

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
                WriteToHDF5(file_id, dset_name, output_arr, Len);
            }

            delete[] output_arr;
        }
    }

    // Write header information for this slice.
    double header_data[5] = {time,
                (double)amrex::ParallelDescriptor::NProcs(),
                (double)(sim->finest_level - sim->shadow_hierarchy),
                (double)sim->dimN[lev],
                (double)le1.size()};
    WriteToHDF5(file_id, "Header_"+ident, header_data, 5);

    // Write box dimensions so we can reassemble slice.
    if (le1.size() > 0) {
        WriteToHDF5(file_id, "le1_"+ident, (int*)&(le1[0]), le1.size());
        WriteToHDF5(file_id, "le2_"+ident, (int*)&(le2[0]), le2.size());
        WriteToHDF5(file_id, "he1_"+ident, (int*)&(he1[0]), he1.size());
        WriteToHDF5(file_id, "he2_"+ident, (int*)&(he2[0]), he2.size());
    }
}

}; // namespace sledgehamr
