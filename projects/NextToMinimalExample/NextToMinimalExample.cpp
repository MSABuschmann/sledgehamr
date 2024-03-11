#include "NextToMinimalExample.h"
#include <hdf5_utils.h>

namespace NextToMinimalExample {

bool NextToMinimalExample::WriteAvg(double time, std::string prefix) {
    // Grab current coarse level.
    const int lev = 0;
    const sledgehamr::LevelData& state = GetLevelData(lev);

    // Create a new state r and compute r = sqrt(Psi1**2 + Psi2**2).
    sledgehamr::LevelData r(state.boxArray(), state.DistributionMap(), 1, 0, 0);

    // Loop over all boxes.
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(r, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.tilebox();
        const auto& state_arr = state.array(mfi);
        const auto& r_arr     = r.array(mfi);

        // Loop over all cells within box.
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            noexcept {
            double Psi1 = state_arr(i, j, k, Scalar::Psi1);
            double Psi2 = state_arr(i, j, k, Scalar::Psi2);
            r_arr(i, j, k, 0) = std::sqrt(Psi1*Psi1 + Psi2*Psi2);
        });
    }

    // Compute average value:
    double avg = r.sum() / r.boxArray().d_numPts();

    // Write information to file.
    constexpr int nsize = 4;
    double data[nsize] = {time, avg};
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/vev.h5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
        sledgehamr::utils::hdf5::Write(file_id, "data", data, nsize);
        H5Fclose(file_id);
    }

    return true;
}

}; // namespace NextToMinimalExample
