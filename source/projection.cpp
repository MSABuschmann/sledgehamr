#include "projection.h"

namespace sledgehamr {

void Projection::Compute(const int id, const hid_t file_id, Sledgehamr* sim) {
    // Could be set to something lower to not include finer level.
    int mlevel = sim->finest_level;

    const int dimN = sim->dimN[mlevel];
    long long N = dimN*dimN;
    double* d_projection = new double[N]();
    int* n_projection = new int[N]();

    std::vector<double> params;
    sim->SetParamsProjections(params);

    for (int lev = 0; lev <= mlevel; ++lev) {
        const int dimN_lev = sim->dimN[lev];
        const int ratio = dimN / dimN_lev;
        const double dx = sim->dx[lev];
        const double dt = sim->dt[lev];
        const double time = sim->grid_new[lev].t;

        amrex::BoxArray ba;
        if (lev != sim->finest_level)
            ba = sim->grid_new[lev+1].boxArray();

        // No OpenMP for thread-safety.
        for (amrex::MFIter mfi(sim->grid_new[lev], false); mfi.isValid();
             ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            const auto& state_fab = sim->grid_new[lev].array(mfi);

            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        bool contd = false;

                        // Only include if cell is not refined.
                        if (lev == sim->finest_level ) {
                            contd = true;
                        } else {
                            if (!ba.contains(amrex::IntVect(i*2, j*2, k*2)))
                                contd = true;
                        }

                        if (contd) {
                            double val = fct(state_fab, i, j, k, lev, time, dt,
                                             dx, params);
                            Add(i, j, k, d_projection, n_projection, ratio,
                                dimN, val);
                        }
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(d_projection, N,
            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceIntSum(n_projection, N,
            amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        if (id == 0) {
            const int nparams = 2;
            hsize_t dims[1] = {nparams};
            double header_data[nparams] = {sim->grid_new[0].t, (double)dimN};
            IOModule::WriteToHDF5(file_id, "Header", header_data, nparams);
        }

        IOModule::WriteToHDF5(file_id, ident + "_data", d_projection, N);
        IOModule::WriteToHDF5(file_id, ident + "_n", n_projection, N);
    }

    delete[] d_projection;
    delete[] n_projection;
}

AMREX_FORCE_INLINE 
void Projection::Add(const int i, const int j, const int k, double* projection,
                     int* n_projection, const int ratio, const int dimN,
                     const double val) {
    // TODO different axis.
    int ii = i*ratio;
    int jj = j*ratio;
    for (int ia = 0; ia < ratio; ++ia) {
        for (int ja = 0; ja < ratio; ++ja) {
            long long ind = (ii + ia)*dimN + (jj+ja);
            if (mode == 0) {
                projection[ind] += ratio*val;
            } else if (mode == 1) {
                projection[ind] = std::max(projection[ind], val);
            }
            n_projection[ind] += 1;
        }
    }
}

}; // namespace sledgehamr
