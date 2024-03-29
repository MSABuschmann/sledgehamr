#include "projection.h"
#include "hdf5_utils.h"

namespace sledgehamr {

/** @brief Computes and saves a line-of-sight projection.
 * @param   id      Output counter id.
 * @param   file_id HDF5 file id.
 * @param   sim     Pointer to the simulation.
 */
void Projection::Compute(const int id, const hid_t file_id, Sledgehamr* sim) {
    int mlevel = INT_MAX;
    amrex::ParmParse pp("output.projections");
    pp.query("max_level", mlevel);
    mlevel = std::min(mlevel, sim->finest_level);

    const int dimN = sim->dimN[mlevel];
    long long N = dimN*dimN;
    amrex::Vector<double> d_projection(N);
    amrex::Vector<int> n_projection(N);

    std::vector<double> params;
    sim->SetParamsProjections(params, sim->grid_new[0].t);

    for (int lev = 0; lev <= mlevel; ++lev) {
        const int dimN_lev = sim->dimN[lev];
        const int ratio = dimN / dimN_lev;
        const double dx = sim->dx[lev];
        const double dt = sim->dt[lev];
        const double time = sim->grid_new[lev].t;

        amrex::BoxArray ba;
        if (lev != mlevel)
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
                        if (lev == mlevel ) {
                            contd = true;
                        } else {
                            if (!ba.contains(amrex::IntVect(i*2, j*2, k*2)))
                                contd = true;
                        }

                        if (contd) {
                            double val = fct(state_fab, i, j, k, lev, time, dt,
                                             dx, params);
                            Add(i, j, k, &d_projection[0], &n_projection[0], ratio,
                                dimN, val);
                        }
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(d_projection.dataPtr(), N,
            amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::ParallelDescriptor::ReduceIntSum(n_projection.dataPtr(), N,
            amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        if (id == 0) {
            const int nparams = 2;
            hsize_t dims[1] = {nparams};
            double header_data[nparams] = {sim->grid_new[0].t, (double)dimN};
            utils::hdf5::Write(file_id, "Header", header_data, nparams);
        }

        utils::hdf5::Write(file_id, ident + "_data", &d_projection[0], N);
        utils::hdf5::Write(file_id, ident + "_n", &n_projection[0], N);
    }
}

/** @brief Projection operator.
 * @param   i               i-th index of cell to project.
 * @param   j               j-th index of cell to project.
 * @param   k               k-th index of cell to project.
 * @param   projection      Existing projection.
 * @param   n_projection    Number of cells projected.
 * @param   ratio           Refinement ratio to coarse level.
 * @param   dimN            This levels number of theoretical cells along on
 *                          axis.
 * @param   val             Cells value.
 */
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
