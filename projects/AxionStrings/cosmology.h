#ifndef PROJECTS_AXION_STRINGS_COSMOLOGY_H_
#define PROJECTS_AXION_STRINGS_COSMOLOGY_H_

#include <sledgehamr.h>

#include "setup.h"
#include "kernels_tagging.h"
#include "kernels_energy_densities.h"

namespace AxionStrings {

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

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_COSMOLOGY_H_
