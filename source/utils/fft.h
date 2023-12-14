#ifndef SLEDGEHAMR_FFT_H_
#define SLEDGEHAMR_FFT_H_

#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

namespace sledgehamr {
namespace utils {

/** @brief This function computes the FFT of some quantity.
 * @param   field                   State to compute the FFT of.
 * @param   field_fft_real_or_abs   Contains the real or absolute part of the
 *                                  FFT.
 * @param   field_fft_imag          Contains the imaginary part of the FFT, if
 *                                  any.
 * @param   geom                    Geometry of the data.
 * @param   abs                     Wheter to keep the real and imaginary part
 *                                  of the FFT or to compute the absolute part.
 */
static void Fft(const amrex::MultiFab& field, const int comp,
                amrex::MultiFab& field_fft_real_or_abs,
                amrex::MultiFab& field_fft_imag, const amrex::Geometry& geom,
                bool abs) {
    const amrex::BoxArray& ba = field.boxArray();
    const amrex::DistributionMapping& dm = field.DistributionMap();
    amrex::MultiFab field_tmp(ba, dm, 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(field, amrex::TilingIfNotGPU());
          mfi.isValid();
          ++mfi ){
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

        constexpr int ALIGN = 16;
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

}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_FFT_H_
