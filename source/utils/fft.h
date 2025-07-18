#ifndef SLEDGEHAMR_FFT_H_
#define SLEDGEHAMR_FFT_H_

// #define OLD_FFT

#include <AMReX_BCUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#ifdef OLD_FFT
#include <AlignedAllocator.h>
#include <Dfft.H>
#include <Distribution.H>
#else
#include <AMReX_FFT.H>
#endif
#include <iterator>

#include "hdf5_utils.h"

namespace sledgehamr {
namespace utils {

static void ChopGrids(amrex::BoxArray &ba, int target_size) {
    const int max_grid_size = 8192;
    const int blocking_factor = 16;

    amrex::IntVect chunk(max_grid_size, max_grid_size, max_grid_size);
    chunk.min(ba.minimalBox().length());

    while (ba.size() < target_size) {
        amrex::IntVect chunk_prev = chunk;

        std::array<std::pair<int, int>, AMREX_SPACEDIM> chunk_dir{AMREX_D_DECL(
            std::make_pair(chunk[0], int(0)), std::make_pair(chunk[1], int(1)),
            std::make_pair(chunk[2], int(2)))};
        std::sort(chunk_dir.begin(), chunk_dir.end());

        for (int idx = AMREX_SPACEDIM - 1; idx >= 0; idx--) {
            int idim = chunk_dir[idx].second;
            int new_chunk_size = chunk[idim] / 2;
            if (new_chunk_size != 0 && new_chunk_size % blocking_factor == 0) {
                chunk[idim] = new_chunk_size;
                ba.maxSize(chunk);
                break;
            }
        }

        if (chunk == chunk_prev) {
            break;
        }
    }
}

/** @brief Writes a field. Debug function for FFT.
 * @param   state               State.
 * @param   file_id             HDF5 file to write to.
 */
static void WriteThis(const amrex::MultiFab *state, hid_t file_id, int comp) {
    const int ndist = 1;
    const int grid_density = 1;
    const int downsample_factor = 1;

    std::vector<int> lex, hex, ley, hey, lez, hez;

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox();
        const auto &state_arr = state->array(mfi);

        // Get box dimensions
        int lx = bx.smallEnd(0);
        int ly = bx.smallEnd(1);
        int lz = bx.smallEnd(2);
        int hx = bx.bigEnd(0) + 1;
        int hy = bx.bigEnd(1) + 1;
        int hz = bx.bigEnd(2) + 1;
        lex.push_back(lx);
        ley.push_back(ly);
        lez.push_back(lz);
        hex.push_back(hx);
        hey.push_back(hy);
        hez.push_back(hz);

        int dimx = (hx - lx) / grid_density;
        int dimy = (hy - ly) / grid_density;
        int dimz = (hz - lz) / grid_density;
        double volfac = 1. / std::pow(downsample_factor, 3);

        // Copy data into flattened array for each scalar field.
        long len = dimx * dimy * dimz;

        const int f = comp;
        // TODO Adjust output type.
        // std::unique_ptr<float[]> output_arr(new float[len]);
        // std::fill_n(output_arr.get(), len, 0.0f);
        std::vector<float> output_arr(len);

        for (int k = lz; k < hz; ++k) {
            for (int j = ly; j < hy; ++j) {
                for (int i = lx; i < hx; ++i) {
                    int ind = (i - lx) / grid_density * dimy * dimz +
                              (j - ly) / grid_density * dimz +
                              (k - lz) / grid_density;

                    output_arr[ind] += state_arr(i, j, k, f) * volfac;
                }
            }
        }

        std::string dset_name = "data_" + std::to_string(lex.size());
        utils::hdf5::Write(file_id, dset_name, &(output_arr[0]), len);
    }

    // Write header information for this slice.
    const int nparams = 3;
    double header_data[nparams] = {
        (double)amrex::ParallelDescriptor::NProcs(),
        (double)state->boxArray().minimalBox().length(0), (double)lex.size()};
    hdf5::Write(file_id, "Header", header_data, nparams);

    // Write box dimensions so we can reassemble slice.
    if (lex.size() == 0)
        return;

    hdf5::Write(file_id, "lex", (int *)&(lex[0]), lex.size());
    hdf5::Write(file_id, "ley", (int *)&(ley[0]), ley.size());
    hdf5::Write(file_id, "lez", (int *)&(lez[0]), lez.size());
    hdf5::Write(file_id, "hex", (int *)&(hex[0]), hex.size());
    hdf5::Write(file_id, "hey", (int *)&(hey[0]), hey.size());
    hdf5::Write(file_id, "hez", (int *)&(hez[0]), hez.size());
}
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
static int fft_counter = 0;
static void Fft(const amrex::MultiFab &field, const int comp,
                amrex::MultiFab &field_fft_real_or_abs,
                amrex::MultiFab &field_fft_imag, const amrex::Geometry &geom,
                bool abs, int zero_padding = 1) {
    const amrex::BoxArray &original_ba = field.boxArray();
    const amrex::Vector<int> &original_pmap =
        field.DistributionMap().ProcessorMap();
    const amrex::Box original_V = original_ba.minimalBox();
    const int N = original_V.length(0);
    const int N_padded = N * zero_padding;

    amrex::Vector<int> tmp_padded_pmap;
    amrex::BoxList tmp_padded_bl;
    for (int i = 0; i < zero_padding; ++i) {
        for (int j = 0; j < zero_padding; ++j) {
            for (int k = 0; k < zero_padding; ++k) {
                amrex::BoxList new_bl = original_ba.boxList();
                new_bl.shift(0, i * N);
                new_bl.shift(1, j * N);
                new_bl.shift(2, k * N);
                tmp_padded_bl.join(new_bl);

                std::copy(original_pmap.begin(), original_pmap.end(),
                          std::back_inserter(tmp_padded_pmap));
            }
        }
    }
    amrex::BoxArray tmp_padded_ba(tmp_padded_bl);
    amrex::DistributionMapping tmp_padded_dm(tmp_padded_pmap);

    amrex::Geometry padded_geom(geom);
    padded_geom.refine(
        amrex::IntVect(zero_padding, zero_padding, zero_padding));

    amrex::MultiFab tmp_field(original_ba, field.DistributionMap(), 1, 0);
    tmp_field.ParallelCopy(field, comp, 0, 1, 0, 0);
    amrex::MultiFab tmp_padded_field(tmp_padded_ba, tmp_padded_dm, 1, 0,
                                     amrex::MFInfo().SetAlloc(false));

    const int offset = original_ba.size();
    for (amrex::MFIter mfi(tmp_padded_field); mfi.isValid(); ++mfi) {
        if (mfi.index() < offset) {
            amrex::FArrayBox &fab = tmp_field[mfi.index()];
            tmp_padded_field.setFab(mfi, std::move(fab));
        } else {
            amrex::FArrayBox fab(mfi.tilebox(), 1);
            fab.setVal<amrex::RunOn::Host>(0.0);
            tmp_padded_field.setFab(mfi, std::move(fab));
        }
    }

    amrex::Box bx = tmp_padded_ba.minimalBox();
    amrex::BoxArray padded_ba = amrex::BoxArray(bx);
    ChopGrids(padded_ba, amrex::ParallelDescriptor::NProcs());
    amrex::DistributionMapping padded_dm(padded_ba,
                                         amrex::ParallelDescriptor::NProcs());

    amrex::MultiFab padded_field(padded_ba, padded_dm, 1, 0);

    amrex::Vector<amrex::BCRec> bcs;
    bcs.resize(1);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        bcs[0].setLo(i, amrex::BCType::int_dir);
        bcs[0].setHi(i, amrex::BCType::int_dir);
    }

    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(padded_geom, bcs,
                                                      bndry_func);

    amrex::Vector<amrex::MultiFab *> smf{
        static_cast<amrex::MultiFab *>(&tmp_padded_field)};
    amrex::Vector<double> stime{0};

    amrex::FillPatchSingleLevel(padded_field, 0, smf, stime, 0, 0, 1,
                                padded_geom, physbc, 0);

#ifndef OLD_FFT
    // Use the new amrex::FFT setup
    amrex::Box domain = padded_ba.minimalBox();
    amrex::FFT::R2C my_fft(domain);

    // create storage for the FFT
    auto const &[cba, cdm] = my_fft.getSpectralDataLayout();

    // amrex::Print() << "padded ba" << padded_ba << std::endl;
    // amrex::Print() << "cba" << cba << std::endl;

    field_fft_real_or_abs.define(cba, cdm, 1, 0);
    field_fft_imag.define(cba, cdm, 1, 0);

    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<amrex::Real>>> phi_fft(
        cba, cdm, 1, 0);
    my_fft.forward(padded_field, phi_fft);

    for (amrex::MFIter mfi(phi_fft); mfi.isValid(); ++mfi) {

        amrex::Array4<amrex::GpuComplex<amrex::Real>> const &phi_fft_ptr =
            phi_fft.array(mfi);
        amrex::Array4<amrex::Real> real_or_abs =
            field_fft_real_or_abs.array(mfi);
        amrex::Array4<amrex::Real> imag = field_fft_imag.array(mfi);

        const amrex::Box &bx = mfi.fabbox();

        if (abs) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                        int k) noexcept {
                real_or_abs(i, j, k, 0) = std::sqrt(
                    phi_fft_ptr(i, j, k).real() * phi_fft_ptr(i, j, k).real() +
                    phi_fft_ptr(i, j, k).imag() * phi_fft_ptr(i, j, k).imag());
            });
        } else {
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    real_or_abs(i, j, k, 0) = phi_fft_ptr(i, j, k).real();
                    imag(i, j, k, 0) = phi_fft_ptr(i, j, k).imag();
                });
        }
    }
#else

    // Now setup SWFFT
    int nx = padded_ba[0].size()[0];
    int ny = padded_ba[0].size()[1];
    int nz = padded_ba[0].size()[2];

    amrex::Box domain(padded_geom.Domain());
    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;

    int nboxes = nbx * nby * nbz;

    amrex::Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    for (int ib = 0; ib < nboxes; ++ib) {
        int i = padded_ba[ib].smallEnd(0) / nx;
        int j = padded_ba[ib].smallEnd(1) / ny;
        int k = padded_ba[ib].smallEnd(2) / nz;

        int local_index = k * nbx * nby + j * nbx + i;

        rank_mapping[local_index] = padded_dm[ib];
    }

    int Ndims[3] = {nbz, nby, nbx};
    int n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD, n, Ndims, &rank_mapping[0]);
    hacc::Dfft dfft(d);

    for (amrex::MFIter mfi(padded_field, false); mfi.isValid(); ++mfi) {
        int gid = mfi.index();

        size_t local_size = dfft.local_size();

        constexpr int ALIGN = 16;
        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN>> a;
        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN>> b;

        a.resize(nx * ny * nz);
        b.resize(nx * ny * nz);

        dfft.makePlans(&a[0], &b[0], &a[0], &b[0]);

        size_t local_indx = 0;
        for (size_t k = 0; k < (size_t)nz; k++) {
            for (size_t j = 0; j < (size_t)ny; j++) {
                for (size_t i = 0; i < (size_t)nx; i++) {
                    complex_t temp(padded_field[mfi].dataPtr()[local_indx], 0.);
                    a[local_indx] = temp;
                    local_indx++;
                }
            }
        }

        dfft.forward(&a[0]);
        d.redistribute_2_to_3(&a[0], &b[0], 2);
        size_t global_size = dfft.global_size();

        local_indx = 0;

        if (abs) {
            for (size_t k = 0; k < (size_t)nz; k++) {
                for (size_t j = 0; j < (size_t)ny; j++) {
                    for (size_t i = 0; i < (size_t)nx; i++) {
                        field_fft_real_or_abs[mfi].dataPtr()[local_indx] =
                            std::abs(b[local_indx]);
                        local_indx++;
                    }
                }
            }
        } else {
            for (size_t k = 0; k < (size_t)nz; k++) {
                for (size_t j = 0; j < (size_t)ny; j++) {
                    for (size_t i = 0; i < (size_t)nx; i++) {
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
#endif
}

}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_FFT_H_
