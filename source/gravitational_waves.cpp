#include "gravitational_waves.h"
#include "fft.h"
#include "hdf5_utils.h"
#include "sledgehamr_utils.h"
#include <AMReX_ParallelReduce.H>

namespace sledgehamr {

/** @brief Constructor that initialized all tensor components needed to
 *         simulate gravitional waves.
 * @param   owner   Pointer to the simulation.
 */
GravitationalWaves::GravitationalWaves(Sledgehamr *owner) {
    sim = owner;
    idx_offset = sim->scalar_fields.size();
    default_modifier = std::make_unique<GravitationalWavesSpectrumModifier>();

    // Don't assume ownership of pointer, it is transferred to
    // sim->scalar_fields.
    ScalarField *u_xx = new ScalarField("u_xx", sim->scalar_fields, false);
    ScalarField *u_yy = new ScalarField("u_yy", sim->scalar_fields, false);
    ScalarField *u_zz = new ScalarField("u_zz", sim->scalar_fields, false);
    ScalarField *u_xy = new ScalarField("u_xy", sim->scalar_fields, false);
    ScalarField *u_xz = new ScalarField("u_xz", sim->scalar_fields, false);
    ScalarField *u_yz = new ScalarField("u_yz", sim->scalar_fields, false);
    ScalarField *du_xx = new ScalarField("du_xx", sim->scalar_fields, true);
    ScalarField *du_yy = new ScalarField("du_yy", sim->scalar_fields, true);
    ScalarField *du_zz = new ScalarField("du_zz", sim->scalar_fields, true);
    ScalarField *du_xy = new ScalarField("du_xy", sim->scalar_fields, true);
    ScalarField *du_xz = new ScalarField("du_xz", sim->scalar_fields, true);
    ScalarField *du_yz = new ScalarField("du_yz", sim->scalar_fields, true);

    amrex::ParmParse pp("");
    std::string param_name = "output.gw_spectra.projection_type";
    pp.query(param_name.c_str(), projection_type);
    utils::ErrorState validity =
        (utils::ErrorState)(projection_type == 2 || projection_type == 3);
    std::string error_msg = "Currently only " + param_name +
                            " = 2 or 3 "
                            "implemented!";
    utils::AssessParam(validity, param_name, projection_type, error_msg, "",
                       sim->nerrors, sim->do_thorough_checks);

    pp.query("output.gw_spectra.zero_padding_factor", zero_padding);
    pp.query("output.gw_spectra.unbinned", unbinned);
}

/** @brief Will compute the gravitional wave spectrum and write the result in
 *         the given hdf5 file.
 * @param   file_id hdf5 file handle.
 */
void GravitationalWaves::ComputeSpectrum(
    hid_t file_id, GravitationalWavesSpectrumModifier *modifier) {
    sim->ReadSpectrumKs();

    if (modifier == nullptr) {
        modifier = default_modifier.get();
    }

    const int lev = 0;
    int dimN = sim->dimN[lev] * zero_padding;
    double L = sim->L * static_cast<double>(zero_padding);

    const LevelData &ld = sim->grid_new[lev];
    amrex::MultiFab du_real[6];
    amrex::MultiFab du_imag[6];
    const int mat[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    int comps[6];
    modifier->SelectComponents(comps);

    std::chrono::steady_clock::time_point start_time =
        std::chrono::steady_clock::now();

    for (int i = 0; i < 6; ++i) {
#ifdef OLD_FFT
        const amrex::BoxArray &ba = ld.boxArray();
        const amrex::DistributionMapping &dm = ld.DistributionMap();
        du_real[i].define(ba, dm, 1, 0);
        du_imag[i].define(ba, dm, 1, 0);
#endif
        utils::Fft(ld, comps[i] + idx_offset, du_real[i], du_imag[i],
                   sim->geom[lev], false, zero_padding);
    }

    std::chrono::steady_clock::time_point end_time =
        std::chrono::steady_clock::now();
    double duration_ms = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count());
    // amrex::Print() << "FFTs: " << duration_ms << std::endl;

    double dk = 2. * M_PI / L;
    double dimN6 = pow(dimN, 6);

    modifier->FourierSpaceModifications(du_real, du_imag, dk, dimN);

    int NTHREADS;
    std::vector<int> ks;
    if (unbinned) {
        ks = sim->gw_spectrum_ks;
        NTHREADS = 16;
    } else {
        for (int k = 0;
             k <= std::sqrt(3.) / 2. * static_cast<double>(dimN) + .5; ++k) {
            ks.push_back(k * k);
        }
        NTHREADS = 16;
    }
    int kmax = ks.size();
    amrex::Gpu::AsyncArray<double> async_index_to_k(sim->index_to_k.data(),
                                                    sim->index_to_k.size());
    double *index_to_k = async_index_to_k.data();

    start_time = std::chrono::steady_clock::now();

#ifdef AMREX_USE_GPU
    if (l_unbinned) {
        amrex::Abort(
            "Unbinned gw spectra on GPUs are currently not supported!");
    }

    bool l_unbinned = unbinned;

    unsigned long SpecLen = kmax;
    std::vector<double> gw_spectrum(SpecLen, 0.0);
    amrex::Gpu::DeviceVector<double> d_data(SpecLen, 0.0);
    double *const AMREX_RESTRICT dptr_data = d_data.dataPtr();

    for (amrex::MFIter mfi(du_real[0], false); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox();

        const amrex::Array4<double> du_real_arr[6] = {
            du_real[0].array(mfi), du_real[1].array(mfi),
            du_real[2].array(mfi), du_real[3].array(mfi),
            du_real[4].array(mfi), du_real[5].array(mfi)};

        const amrex::Array4<double> du_imag_arr[6] = {
            du_imag[0].array(mfi), du_imag[1].array(mfi),
            du_imag[2].array(mfi), du_imag[3].array(mfi),
            du_imag[4].array(mfi), du_imag[5].array(mfi)};

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int a, int b,
                                                    int c) noexcept {
            // To account of negative frequencies
            double multpl = (a == 0 || a == dimN / 2) ? 1. : 2.;
            int abc[3] = {a, b, c};
            double running_sum = 0;

            int li = a >= dimN / 2 ? a - dimN : a;
            int lj = b >= dimN / 2 ? b - dimN : b;
            int lk = c >= dimN / 2 ? c - dimN : c;
            unsigned int sq = li * li + lj * lj + lk * lk;
            unsigned long index;
            if (l_unbinned) {
                // index = std::lower_bound(ks.begin(), ks.end(), sq) -
                // ks.begin();
            } else {
                index =
                    static_cast<long>(std::sqrt(static_cast<double>(sq)) + .5);
            }

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int l = 0; l < 3; ++l) {
                        for (int m = 0; m < 3; ++m) {
                            running_sum +=
                                gw_GetLambda(i, j, l, m, abc, index_to_k) *
                                multpl *
                                (du_real_arr[mat[i][j]].operator()(a, b, c) *
                                     du_real_arr[mat[l][m]].operator()(a, b,
                                                                       c) +
                                 du_imag_arr[mat[i][j]].operator()(a, b, c) *
                                     du_imag_arr[mat[l][m]].operator()(a, b,
                                                                       c));
                        }
                    }
                }
            }
            amrex::HostDevice::Atomic::Add(&dptr_data[index],
                                           running_sum / dimN6);
        });
    }

    end_time = std::chrono::steady_clock::now();
    duration_ms = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count());
    // amrex::Print() << "Sum: " << duration_ms << std::endl;
    start_time = std::chrono::steady_clock::now();

    // blocking copy from device to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_data.begin(), d_data.end(),
                     gw_spectrum.begin());

    // reduced sum over mpi ranks
    amrex::ParallelDescriptor::ReduceRealSum(
        gw_spectrum.data(), static_cast<int>(gw_spectrum.size()),
        amrex::ParallelDescriptor::IOProcessorNumber());

#else // ifdef AMREX_USE_GPU

    unsigned long SpecLen = kmax * NTHREADS;
    std::vector<double> gw_spectrum(SpecLen, 0.0);

#pragma omp parallel num_threads(std::min(NTHREADS, omp_get_max_threads()))
    for (amrex::MFIter mfi(du_real[0], true); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox();

        amrex::Array4<double> const &dr0 = du_real[0].array(mfi);
        amrex::Array4<double> const &dr1 = du_real[1].array(mfi);
        amrex::Array4<double> const &dr2 = du_real[2].array(mfi);
        amrex::Array4<double> const &dr3 = du_real[3].array(mfi);
        amrex::Array4<double> const &dr4 = du_real[4].array(mfi);
        amrex::Array4<double> const &dr5 = du_real[5].array(mfi);
        amrex::Array4<double> const &di0 = du_imag[0].array(mfi);
        amrex::Array4<double> const &di1 = du_imag[1].array(mfi);
        amrex::Array4<double> const &di2 = du_imag[2].array(mfi);
        amrex::Array4<double> const &di3 = du_imag[3].array(mfi);
        amrex::Array4<double> const &di4 = du_imag[4].array(mfi);
        amrex::Array4<double> const &di5 = du_imag[5].array(mfi);

        const amrex::Array4<double> *du_real_arr[6] = {
            std::addressof(dr0), std::addressof(dr1), std::addressof(dr2),
            std::addressof(dr3), std::addressof(dr4), std::addressof(dr5)};
        const amrex::Array4<double> *du_imag_arr[6] = {
            std::addressof(di0), std::addressof(di1), std::addressof(di2),
            std::addressof(di3), std::addressof(di4), std::addressof(di5)};

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int c = lo.z; c <= hi.z; ++c) {
            for (int b = lo.y; b <= hi.y; ++b) {
                AMREX_PRAGMA_SIMD
                for (int a = lo.x; a <= hi.x; ++a) {
                    // To account of negative frequencies
                    double multpl = (a == 0 || a == dimN / 2) ? 1. : 2.;
                    int abc[3] = {a, b, c};
                    double running_sum = 0;

                    int li = a >= dimN / 2 ? a - dimN : a;
                    int lj = b >= dimN / 2 ? b - dimN : b;
                    int lk = c >= dimN / 2 ? c - dimN : c;
                    unsigned int sq = li * li + lj * lj + lk * lk;
                    unsigned long index;
                    if (unbinned) {
                        index = std::lower_bound(ks.begin(), ks.end(), sq) -
                                ks.begin() + omp_get_thread_num() * kmax;
                    } else {
                        index = static_cast<long>(
                            std::sqrt(static_cast<double>(sq)) + .5);
                    }

                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            for (int l = 0; l < 3; ++l) {
                                for (int m = 0; m < 3; ++m) {
                                    running_sum +=
                                        gw_GetLambda(i, j, l, m, abc,
                                                     index_to_k) *
                                        multpl *
                                        (du_real_arr[mat[i][j]]->operator()(
                                             a, b, c) *
                                             du_real_arr[mat[l][m]]->operator()(
                                                 a, b, c) +
                                         du_imag_arr[mat[i][j]]->operator()(
                                             a, b, c) *
                                             du_imag_arr[mat[l][m]]->operator()(
                                                 a, b, c));
                                }
                            }
                        }
                    }
                    gw_spectrum[index] += running_sum;
                }
            }
        }
    }

    end_time = std::chrono::steady_clock::now();
    duration_ms = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count());
    // amrex::Print() << "Sum: " << duration_ms << std::endl;
    start_time = std::chrono::steady_clock::now();

    for (int a = 1; a < NTHREADS; ++a) {
        for (int c = 0; c < kmax; ++c) {
            gw_spectrum[c] += gw_spectrum[a * kmax + c];
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(
        &(gw_spectrum[0]), kmax,
        amrex::ParallelDescriptor::IOProcessorNumber());

#pragma omp parallel for
    for (int c = 0; c < kmax; ++c) {
        gw_spectrum[c] /= dimN6;
    }

#endif // AMREX_USE_GPU
    end_time = std::chrono::steady_clock::now();
    duration_ms = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count());
    // amrex::Print() << "Reduce: " << duration_ms << std::endl;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const int nparams = 6;
        double header_data[nparams] = {
            ld.t, (double)dimN,         (double)kmax,
            L,    (double)zero_padding, (double)unbinned};
        utils::hdf5::Write(file_id, "Header", header_data, nparams);
        utils::hdf5::Write(file_id, "k", &(ks[0]), kmax);
        utils::hdf5::Write(file_id, "Spectrum", &gw_spectrum[0], kmax);
    }
}

}; // namespace sledgehamr
