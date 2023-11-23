#include "spectrum.h"
#include "hdf5_utils.h"
#include "fft.h"

namespace sledgehamr {

void Spectrum::Compute(const int id, const hid_t file_id, Sledgehamr* sim) {
    amrex::Print() << "Compute Spectrum: " << ident << std::endl;

    const int lev = 0;
    const int dimN = sim->dimN[lev];
    const double dx = sim->dx[lev];
    const double dt = sim->dt[lev];
    const double time = sim->grid_new[lev].t;
    const LevelData& state = sim->grid_new[lev];
    const amrex::BoxArray& ba = state.boxArray();

    amrex::MultiFab field, field_fft;
    field.define(ba, sim->dmap[lev], 1, 0);
    field_fft.define(ba, sim->dmap[lev], 1, 0);

    std::vector<double> params;
    sim->SetParamsSpectra(params, time);

#pragma omp parallel
    for (amrex::MFIter mfi(field, true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& field_arr = field.array(mfi);
        const auto& state_arr = state.array(mfi);

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    field_arr(i, j, k, 0) = fct(state_arr, i, j, k, lev, time,
                                                dt, dx, params);
                }
            }
        }
    }

    utils::Fft(field, 0, field_fft, field_fft, sim->geom[lev], true);

    double fac = pow(1./dimN, 6);
    double dk = 2.*M_PI / sim->L;
    double pre = fac*state.t/dk;

    std::vector<int>& ks = sim->spectrum_ks;
    const int kmax = ks.size();
    constexpr int NTHREADS = 16;
    const unsigned long SpecLen = kmax*NTHREADS;
    std::unique_ptr<double[]> spectrum(new double[SpecLen]);

#pragma omp parallel num_threads(std::min(NTHREADS, omp_get_max_threads()))
    for (amrex::MFIter mfi(field_fft, true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& field_fft_arr = field_fft.array(mfi);

        const int il = bx.smallEnd()[0];
        const int ih = bx.bigEnd()[0];
        const int jl = bx.smallEnd()[1];
        const int jh = bx.bigEnd()[1];
        const int kl = bx.smallEnd()[2];
        const int kh = bx.bigEnd()[2];

        for (int i = il; i <= ih; ++i) {
            for (int j = jl; j <= jh; ++j) {
                AMREX_PRAGMA_SIMD
                for (int k = kl; k <= kh; ++k) {
                    int li = i >= dimN/2 ? i-dimN : i;
                    int lj = j >= dimN/2 ? j-dimN : j;
                    int lk = k >= dimN/2 ? k-dimN : k;
                    unsigned int sq = li*li + lj*lj + lk*lk;
                    unsigned long index =
                            std::lower_bound(ks.begin(), ks.end(), sq)
                            - ks.begin() + omp_get_thread_num() * kmax;
                    spectrum[index] += pre * field_fft_arr(i,j,k,0)
                                           * field_fft_arr(i,j,k,0);
                }
            }
        }
     }

    for (int a = 1; a < NTHREADS; ++a) {
        for (int c = 0; c < kmax; ++c) {
            spectrum[c] += spectrum[a*kmax + c];
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(spectrum.get(), kmax,
            amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        if (id == 0) {
            const int nparams = 3;
            double header_data[nparams] = {sim->grid_new[lev].t, (double)dimN,
                                           (double)kmax};
            utils::hdf5::Write(file_id, "Header", header_data, nparams);
            utils::hdf5::Write(file_id, "k_sq", &(ks[0]), kmax);
        }

        utils::hdf5::Write(file_id, ident, spectrum.get(), kmax);
    }
}

}; // namespace sledgehamr
