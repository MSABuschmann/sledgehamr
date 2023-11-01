#ifndef PROJECTS_AXION_STRINGS_COSMOLOGY_H_
#define PROJECTS_AXION_STRINGS_COSMOLOGY_H_

#include <sledgehamr.h>

namespace axion_strings {

ADD_SCALARS(Psi1, Psi2)
ADD_CONJUGATE_MOMENTA(Pi1, Pi2)

/** @brief TODO
 */
AMREX_FORCE_INLINE
double a_prime2(amrex::Array4<amrex::Real const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double r2      = Psi1*Psi1 + Psi2*Psi2;
    double prime_a = (Psi1*Pi2 - Psi2*Pi1)/r2;
    return prime_a*prime_a;
}

AMREX_FORCE_INLINE
double a_prime_screened(amrex::Array4<amrex::Real const> const& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double prime_a = (Psi1*Pi2 - Psi2*Pi1);
    return prime_a;
}

AMREX_FORCE_INLINE
double r_prime2(amrex::Array4<amrex::Real const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double r2      = Psi1*Psi1 + Psi2*Psi2;
    double prime_r = (Psi1*Pi1 + Psi2*Pi2);
    return prime_r*prime_r/r2;
}

/** @brief Checks for zero-crossings between two points in the complex scalar
 *         field.
 * @param   Psi1_1  \Psi_1 of first point.
 * @param   Psi2_1  \Psi_2 of first point.
 * @param   Psi1_2  \Psi_1 of second point.
 * @param   Psi2_2  \Psi_2 of second point.
 * @return Sign of slope of zero-crossing. 0 if no crossing.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ZeroXing(double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2) {
    if (Psi2_1 * Psi2_2 >= 0) return 0;
    if (Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0) return 1;
    return -1;
}

/** @brief Computes the winding factor along a given axis. Will be non-zero if
 *         plaquette is pierced by a string.
 * @return  Winding factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis1(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2),
                    state(i+1,j+1,k  ,Scalar::Psi1),
                    state(i+1,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j+1,k  ,Scalar::Psi1),
                    state(i+1,j+1,k  ,Scalar::Psi2),
                    state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis2(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k+1,Scalar::Psi1),
                    state(i+1,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k+1,Scalar::Psi1),
                    state(i+1,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis3(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2),
                    state(i  ,j+1,k+1,Scalar::Psi1),
                    state(i  ,j+1,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k+1,Scalar::Psi1),
                    state(i  ,j+1,k+1,Scalar::Psi2),
                    state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

/** @brief Function that tags individual cells for refinement.
 * @return  Boolean value as to whether cell should be refined or not.
 */
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool TagCellForRefinement<true>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double* params) {
    // Check all three plaquettes (in positive index direction) for string
    // piercings.
    if (WindingAxis1(state, i, j, k) != 0) return true;
    if (WindingAxis2(state, i, j, k) != 0) return true;
    if (WindingAxis3(state, i, j, k) != 0) return true;

    return false;
};


/** @brief Class to simulate axion strings.
 */
class Cosmology {
  public:
    void Init(sledgehamr::Sledgehamr* owner);

    bool CreateLevelIf(const int lev, const double time) {
        return StringWidth(lev-1, time) <= string_width_threshold;
    }

    double Mr(const double eta) {
        return std::sqrt(2. * lambda) * eta;
    }

    double H(const double eta) {
        return 1./eta;
    }

    double StringWidth(const int lev, const double eta) {
        return 1./(Mr(eta) * sim->GetDx(lev));
    }

    double RefinementTime(const int lev) {
        return sim->GetDimN(lev) / (sqrt(2.*lambda) 
                * string_width_threshold * sim->GetL());
    }

    double Log(const double eta) {
        if (eta <= 0)
            return -DBL_MAX;
        return std::log( Mr(eta) / H(eta) );
    }

    double LogTruncated(const double eta) {
        double log = Log(eta);
        if (log < spectra_log_min) 
            return 0;
        return log;
    }

    double BoxToPhysical(const double L, const double eta, const double T1,
                         const double mpl, const double gStar) {
        return L * eta / Hubble(T1, mpl, gStar);
    }

    double Hubble(const double T, const double mpl, const double gStar) {
        return std::sqrt(4.*std::pow(M_PI, 3) / 45. * gStar * std::pow(T, 4)
                / std::pow(mpl, 2));
    }

    double XiTime(double T, double mpl, double gStar) {
        return 0.3012 / std::sqrt(gStar) * mpl / std::pow(T, 2);
    }

    double XiTemp(double eta, double T1) {
        return T1 / eta;
    }

    double Xi(const int lev, const double eta);

  private:
    void ParseVariables();
    void PrintRefinementTimes();
    void SetProjections();
    void SetSpectra();
    void SetXiMeasurement();

    int GetStringTags(const int lev);
    bool WriteXi(double time, std::string prefix);

    double string_width_threshold;
    double spectra_log_min = 5;
    double interval_xi_log = 0;
    const double lambda = 1;

    sledgehamr::Sledgehamr* sim;
};

}; // namespace axion_strings

#endif // PROJECTS_AXION_STRINGS_COSMOLOGY_H_
