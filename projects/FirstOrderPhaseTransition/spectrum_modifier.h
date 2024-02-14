#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SPECTRUM_MODIFIER_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SPECTRUM_MODIFIER_H_

#include <gravitational_waves.h>

namespace FirstOrderPhaseTransition {

struct SpectrumModifier_UtimesK
        : sledgehamr::GravitationalWavesSpectrumModifier {
     virtual void SelectComponents(int components[6]) {
        components[0] = sledgehamr::GravitationalWaves::Gw::u_xx;
        components[1] = sledgehamr::GravitationalWaves::Gw::u_xy;
        components[2] = sledgehamr::GravitationalWaves::Gw::u_xz;
        components[3] = sledgehamr::GravitationalWaves::Gw::u_yy;
        components[4] = sledgehamr::GravitationalWaves::Gw::u_yz;
        components[5] = sledgehamr::GravitationalWaves::Gw::u_zz;
    };

    virtual void FourierSpaceModifications(
            amrex::MultiFab du_real[6], amrex::MultiFab du_imag[6],
            const double dk, const int dimN) {

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        for (amrex::MFIter mfi(du_real[0], true); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
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
                        int kx = a >= dimN/2 ? a-dimN : a;
                        int ky = b >= dimN/2 ? b-dimN : b;
                        int kz = c >= dimN/2 ? c-dimN : c;
                        double k = std::sqrt(kx*kx + ky*ky + kz*kz) * dk;

                        for (int i = 0; i < 6; ++i) {
                            (*du_real_arr)[i](a, b, c) *= k;
                            (*du_imag_arr)[i](a, b, c) *= k;
                        }
                    }
                }
            }
        }
    };
};

struct SpectrumModifier_2BubblesFrom1
        : sledgehamr::GravitationalWavesSpectrumModifier {
    SpectrumModifier_2BubblesFrom1(const double distance[3]) {
        std::memcpy(d, distance, sizeof(d));
    }

    virtual void FourierSpaceModifications(
            amrex::MultiFab du_real[6], amrex::MultiFab du_imag[6],
            const double dk, const int dimN) {

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        for (amrex::MFIter mfi(du_real[0], true); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
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
                        int kx = a >= dimN/2 ? a-dimN : a;
                        int ky = b >= dimN/2 ? b-dimN : b;
                        int kz = c >= dimN/2 ? c-dimN : c;
                        double kd = (kx*d[0] + ky*d[1] + kz*d[2]) * dk;
                        double sx = std::sin(kd);
                        double cx = std::cos(kd);

                        for (int i = 0; i < 6; ++i) {
                            double dur = (*du_real_arr)[i](a, b, c);
                            double dui = (*du_imag_arr)[i](a, b, c);
                            (*du_real_arr)[i](a,b,c) = dur + (cx*dur - sx*dui);
                            (*du_imag_arr)[i](a,b,c) = dui + (sx*dur + cx*dui);
                        }
                    }
                }
            }
        }
    };

    double d[3] = {0,0,0};
};

};

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SPECTRUM_MODIFIER_H_
