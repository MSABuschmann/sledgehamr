#include "gravitational_waves.h"

namespace sledgehamr {

GravitationalWaves::GravitationalWaves(Sledgehamr* owner) {
    sim = owner;
    idx_offset = sim->scalar_fields.size();

    ScalarField* u_xx = new ScalarField("u_xx", sim->scalar_fields, false);
    ScalarField* u_yy = new ScalarField("u_yy", sim->scalar_fields, false);
    ScalarField* u_zz = new ScalarField("u_zz", sim->scalar_fields, false);
    ScalarField* u_xy = new ScalarField("u_xy", sim->scalar_fields, false);
    ScalarField* u_xz = new ScalarField("u_xz", sim->scalar_fields, false);
    ScalarField* u_yz = new ScalarField("u_yz", sim->scalar_fields, false);
    ScalarField* du_xx = new ScalarField("du_xx", sim->scalar_fields, true);
    ScalarField* du_yy = new ScalarField("du_yy", sim->scalar_fields, true);
    ScalarField* du_zz = new ScalarField("du_zz", sim->scalar_fields, true);
    ScalarField* du_xy = new ScalarField("du_xy", sim->scalar_fields, true);
    ScalarField* du_xz = new ScalarField("du_xz", sim->scalar_fields, true);
    ScalarField* du_yz = new ScalarField("du_yz", sim->scalar_fields, true);
}

void GravitationalWaves::ComputeSpectrum(hid_t file_id) {
    const int lev = 0;
    int dimN = sim->dimN[lev];

    const LevelData& ld = sim->grid_new[lev];
    const amrex::BoxArray& ba = ld.boxArray();
    const amrex::DistributionMapping& dm = ld.DistributionMap();

    amrex::MultiFab du_real[6], du_imag[6];
    const int mat[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    const int comps[6] = {Gw::du_xx, Gw::du_xy, Gw::du_xz,
                          Gw::du_yy, Gw::du_yz, Gw::du_zz};

    for (int i = 0; i < 6; ++i) {
        du_real[i].define(ba, dm, 1, 0);
        du_imag[i].define(ba, dm, 1, 0);
        Spectrum::Fft(ld, comps[i] + idx_offset, du_real[i], du_imag[i],
                      sim->geom[lev], false);
    }

    double dk = 2.*M_PI / sim->L;
    double dimN6 = pow(dimN, 6);

    std::vector<int>& ks = sim->spectrum_ks;
    const int kmax = ks.size();
    constexpr int NTHREADS = 16;
    const unsigned long SpecLen = kmax*NTHREADS;
    double* gw_spectrum = new double [SpecLen] ();

    // Non-trivial load-balancing here. Not sure what wins.
#pragma omp parallel num_threads(NTHREADS)
    for (amrex::MFIter mfi(du_real[0], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        // Ugly work around since arrays of references are not allowed.
        amrex::Array4<double> const& dr0 = du_real[0].array(mfi);
        amrex::Array4<double> const& dr1 = du_real[1].array(mfi);
        amrex::Array4<double> const& dr2 = du_real[2].array(mfi);
        amrex::Array4<double> const& dr3 = du_real[3].array(mfi);
        amrex::Array4<double> const& dr4 = du_real[4].array(mfi);
        amrex::Array4<double> const& dr5 = du_real[5].array(mfi);
        amrex::Array4<double> const& di0 = du_imag[0].array(mfi);
        amrex::Array4<double> const& di1 = du_imag[1].array(mfi);
        amrex::Array4<double> const& di2 = du_imag[2].array(mfi);
        amrex::Array4<double> const& di3 = du_imag[3].array(mfi);
        amrex::Array4<double> const& di4 = du_imag[4].array(mfi);
        amrex::Array4<double> const& di5 = du_imag[5].array(mfi);
        const amrex::Array4<double>* du_real_arr[6] = {&dr0, &dr1, &dr2,
                                                       &dr3, &dr4, &dr5};
        const amrex::Array4<double>* du_imag_arr[6] = {&di0, &di1, &di2,
                                                       &di3, &di4, &di5};

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int c = lo.z; c <= hi.z; ++c) {
            for (int b = lo.y; b <= hi.y; ++b) {
                AMREX_PRAGMA_SIMD
                for (int a = lo.x; a <= hi.x; ++a) {
                    int abc[3] = {a, b, c};
                    double running_sum = 0;

                    for(int i = 0; i < 3; ++i) {
                        for(int j = 0; j < 3; ++j) {
                            for(int l = 0; l < 3; ++l) {
                                for(int m = 0; m < 3; ++m) {
                                    running_sum +=
                                        GetLambda(i, j, l, m, abc, dimN)
                                        * ( (*du_real_arr)[mat[i][j]](a, b, c)
                                          * (*du_real_arr)[mat[l][m]](a, b, c)
                                          + (*du_imag_arr)[mat[i][j]](a ,b, c)
                                          * (*du_imag_arr)[mat[l][m]](a, b, c));
                                }
                            }
                        }
                    }

                    int li = a >= dimN/2 ? a-dimN : a;
                    int lj = b >= dimN/2 ? b-dimN : b;
                    int lk = c >= dimN/2 ? c-dimN : c;
                    unsigned int sq= li*li + lj*lj + lk*lk;
                    unsigned long index =
                            std::lower_bound(ks.begin(), ks.end(), sq) -
                            ks.begin() + omp_get_thread_num() * kmax;
                    gw_spectrum[index] += running_sum;
                }
            }
        }
    }

    for (int a = 1; a < NTHREADS; ++a) {
        for (int c = 0; c < kmax; ++c) {
            gw_spectrum[c] += gw_spectrum[a*kmax + c];
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(gw_spectrum, kmax,
            amrex::ParallelDescriptor::IOProcessorNumber());

#pragma omp parallel for
    for (int c = 0; c < kmax; ++c) {
        gw_spectrum[c] /= dimN6;
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const int nparams = 3;
        double header_data[nparams] = {ld.t, (double)dimN, (double)kmax};
        IOModule::WriteToHDF5(file_id, "Header", header_data, nparams);
        IOModule::WriteToHDF5(file_id, "k", &(ks[0]), kmax);
        IOModule::WriteToHDF5(file_id, "Spectrum", gw_spectrum, kmax);
    }

    delete[] gw_spectrum;
}

inline double GravitationalWaves::IndexToK(int a, int N) {
    double n_tilde = a-N <= -N/2-1 ? a : a-N;
    double two_pi_n_tilde = 2.*M_PI/static_cast<double>(N)*n_tilde;

    if (sim->nghost == 1) {
        return sin(two_pi_n_tilde);
    } else if (sim->nghost > 1) {
        return (8.*sin(two_pi_n_tilde) - sin(2.*two_pi_n_tilde)) / 6.;
    }

    return 0.;
}

inline double GravitationalWaves::GetProjection(int i, int j, int abc[3],
                                                int N) {
    if( abc[0] == 0 && abc[1] == 0 && abc[2] == 0)
        return 0.;

    double abc_d[3];
    abc_d[0] = IndexToK(abc[0], N);
    abc_d[1] = IndexToK(abc[1], N);
    abc_d[2] = IndexToK(abc[2], N);

    double norm = abc_d[0]*abc_d[0] + abc_d[1]*abc_d[1] + abc_d[2]*abc_d[2];
    int l = amrex::min(i, j);
    int m = amrex::max(i, j);

    double proj = abc_d[l] * abc_d[m];

    return static_cast<double>(l==m) - proj/norm;
}

inline double GravitationalWaves::GetLambda(int i, int j, int l, int m,
                                            int abc[3], int N) {
    return GetProjection(i, l, abc, N) * GetProjection(j, m, abc, N)
         - GetProjection(i, j, abc, N) * GetProjection(l, m, abc, N) / 2.;
}

}; // namespace sledgehamr
