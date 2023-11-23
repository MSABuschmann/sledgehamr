#include "fill_level.h"
#include "hdf5_utils.h"
#include "sledgehamr_utils.h"
#include "checkpoint.h"

namespace sledgehamr {

void FillLevel::FromInitialStateFile() {
    amrex::ParmParse pp("");
    std::string initial_state_file = "";
    std::string param_name = "input.initial_state"; 
    pp.query(param_name.c_str(), initial_state_file);
    utils::AssessParamOK(param_name, initial_state_file,
                         sim->do_thorough_checks);

    if (amrex::FileExists(initial_state_file + "/Meta.hdf5")) {
        FromCheckpointFile(initial_state_file);
    } else {
        FromHdf5File(initial_state_file);
    }
}

void FillLevel::FromCheckpointFile(std::string folder) {
    Checkpoint chk(sim, folder);
    chk.Read();

    bool delete_restart_checkpoint = false;
    amrex::ParmParse pp("");
    pp.query("input.delete_restart_checkpoint", delete_restart_checkpoint);

    if (delete_restart_checkpoint)  {
        sim->io_module->old_checkpoint = folder;
    }
}

void FillLevel::FromHdf5File(std::string initial_state_file) {
    amrex::ParmParse pp("input");

    // Iterate over fields but introduce offset such that each node grabs a
    // different file first.
    const int ncomp = sim->GetLevelData(lev).nComp();
    int lr = amrex::ParallelDescriptor::MyProc();
    int mr = amrex::ParallelDescriptor::NProcs();
    for (int f=0; f<ncomp; ++f) {
        int f2 = (f+lr) % ncomp;

        // Determine which input file to use if any.
        std::string initial_state_file_component = "";
        std::string scalar_name = sim->GetScalarFieldName(f2);
        std::string query = "initial_state_" + scalar_name;
        pp.query(query.c_str(), initial_state_file_component);

        if (initial_state_file_component == "")
            initial_state_file_component = initial_state_file;

        // If no file found, fill level with 0's. Otherwise read file and fill
        // LevelData.
        if(initial_state_file_component == "") {
            FromConst(f2, 0);
        } else {
            amrex::Print() << "Reading initial state for " << scalar_name
                           << " from " << initial_state_file_component
                           << std::endl;

            // Test if chunks exist.
            std::string chunk1 = scalar_name + "_" + std::to_string(lr);
            std::string chunk2 = "data_" + std::to_string(lr);
            std::string existing_chunk = utils::hdf5::FindDataset(
                    initial_state_file_component, {chunk1, chunk2});

            if (existing_chunk == "") {
                const int dimN = sim->GetDimN(lev);
                const long long dsetsize = dimN*dimN*dimN;
                double* input_data = new double [dsetsize] ();

                if (!utils::hdf5::Read(initial_state_file_component,
                        {scalar_name, "data"}, input_data)) {
                    //std::string msg =
                    //           "Sledgehamr::IOModule::FillLevelFromHdf5File: "
                    //           "Could not find initial state data for "
                    //         + scalar_name + "!";
                    //amrex::Abort(msg);
                    const int constant = 0;
                    amrex::Print() << "Dataset not found for " << scalar_name
                                   << ". Will initialize to " << constant
                                   << "." << std::endl;
                    FromConst(f2, constant);
                } else {
                    FromArray(f2, input_data, dimN);
                }

                delete[] input_data;
            } else {
                amrex::Box bx = sim->GetLevelData(lev).boxArray()[lr];
                const long long dsetsize= bx.numPts();
                double* input_data = new double [dsetsize] ();

               if (!utils::hdf5::Read(initial_state_file_component,
                        {chunk1, chunk2}, input_data)) {
                    std::string msg =
                               "Sledgehamr::IOModule::FillLevelFromHdf5File: "
                               "Could not find initial state chunk "
                             + chunk1 + "!";
                    amrex::Abort(msg);
                }

                FromArrayChunks(f2, input_data);
                delete[] input_data;
            }
        }
    }
}

void FillLevel::FromArray(const int comp, double* data, const long long dimN) {
    LevelData& state = sim->GetLevelData(lev);

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

void FillLevel::FromArrayChunks(const int comp, double* data) {
    LevelData& state = sim->GetLevelData(lev);

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

void FillLevel::FromConst(const int comp, const double c) {
    LevelData& state = sim->GetLevelData(lev);

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

}; // namespace sledgehamr
