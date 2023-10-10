#include "spectrum.h"

#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

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
    sim->SetParamsSpectra(params);

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

    Fft(field, 0, field_fft, field_fft, sim->geom[lev], true);

    double fac = pow(1./dimN, 6);
    double dk = 2.*M_PI / sim->L;
    double pre = fac*state.t/dk;

    std::vector<int>& ks = sim->spectrum_ks;
    const int kmax = ks.size();
    constexpr int NTHREADS = 16;
    const unsigned long SpecLen = kmax*NTHREADS;
    double* spectrum = new double [SpecLen] ();

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

    amrex::ParallelDescriptor::ReduceRealSum(spectrum, kmax,
            amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        if (id == 0) {
            const int nparams = 3;
            double header_data[nparams] = {sim->grid_new[lev].t, (double)dimN,
                                           (double)kmax};
            IOModule::WriteToHDF5(file_id, "Header", header_data, nparams);
            IOModule::WriteToHDF5(file_id, "k_sq", &(ks[0]), kmax);
        }

        IOModule::WriteToHDF5(file_id, ident, spectrum, kmax);
    }

    delete[] spectrum;
}

#define ALIGN 16
void Spectrum::Fft(const amrex::MultiFab& field, const int comp,
                   amrex::MultiFab& field_fft_real_or_abs,
                   amrex::MultiFab& field_fft_imag, const amrex::Geometry& geom,
                   bool abs) {
    const amrex::BoxArray& ba = field.boxArray();
    const amrex::DistributionMapping& dm = field.DistributionMap();
    amrex::MultiFab field_tmp(ba, dm, 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ){
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = field.array(mfi);
        const auto& field_tmp_arr = field_tmp.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            field_tmp_arr(i, j, k, 0) = state_arr(i, j, k, comp);
        });
    }

    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];

    amrex::Box domain(geom.Domain());
    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;

    int nboxes = nbx * nby * nbz;

    amrex::Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    for (int ib = 0; ib < nboxes; ++ib) {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
        int k = ba[ib].smallEnd(2) / nz;

        int local_index = k*nbx*nby + j*nbx + i;

        rank_mapping[local_index] = dm[ib];
    }

    int Ndims[3] = { nbz, nby, nbx };
    int n[3] = { domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD, n, Ndims, &rank_mapping[0]);
    hacc::Dfft dfft(d);

    for (amrex::MFIter mfi(field, false); mfi.isValid(); ++mfi) {
        int gid = mfi.index();

        size_t local_size  = dfft.local_size();

        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

        a.resize(nx*ny*nz);
        b.resize(nx*ny*nz);

        dfft.makePlans(&a[0], &b[0], &a[0], &b[0]);

        size_t local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
            for(size_t j=0; j<(size_t)ny; j++) {
                for(size_t i=0; i<(size_t)nx; i++) {
                    complex_t temp(field_tmp[mfi].dataPtr()[local_indx],0.);
                    a[local_indx] = temp;
                    local_indx++;
                }
            }
        }

        dfft.forward(&a[0]);
        d.redistribute_2_to_3(&a[0], &b[0], 2);
        size_t global_size  = dfft.global_size();

        local_indx = 0;

        if (abs) {
            for(size_t k=0; k<(size_t)nz; k++) {
                for(size_t j=0; j<(size_t)ny; j++) {
                    for(size_t i=0; i<(size_t)nx; i++) {
                        field_fft_real_or_abs[mfi].dataPtr()[local_indx] =
                                std::abs(b[local_indx]);
                        local_indx++;
                    }
                }
            }
        } else {
            for(size_t k=0; k<(size_t)nz; k++) {
                for(size_t j=0; j<(size_t)ny; j++) {
                    for(size_t i=0; i<(size_t)nx; i++) {
                        field_fft_real_or_abs[mfi].dataPtr()[local_indx] =
                                std::real(b[local_indx]);
                        field_fft_imag[mfi].dataPtr()[local_indx] =
                                std::imag(b[local_indx]);
                        local_indx++;
                    }
                }
            }
        }
    }
}

}; // namespace sledgehamr
