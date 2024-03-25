#ifndef PROJECTS_AXION_STRINGS_COSMOLOGY_H_
#define PROJECTS_AXION_STRINGS_COSMOLOGY_H_

#include <sledgehamr.h>

#include "setup.h"
#include "kernels_tagging.h"
#include "kernels_energy_densities.h"

namespace AxionStrings {

/** @brief Class to simulate axion strings. This is not a sledgehamr class but
 *         will set up the axion string scenario if added to the project class.
 */
class Cosmology {
  public:
    void Init(sledgehamr::Sledgehamr* owner);

    /** @brief We only want to create a new level if string width is below
     *         threshold.
     * @param   lev     Current level.
     * @param   time    Current time.
     * @return  Whether we want to create level lev.
     */
    bool CreateLevelIf(const int lev, const double time) {
        return StringWidth(lev-1, time) <= string_width_threshold;
    }

    /** @brief Current radial mode mass.
     * @param   eta Current time eta.
     * @return m_r(eta)
     */
    double Mr(const double eta) {
        return std::sqrt(2. * lambda) * eta;
    }

    /** @brief Current hubble time.
     * @param   eta Current time eta.
     * @return H(eta).
     */
    double H(const double eta) {
        return 1./eta;
    }

    /** @brief Current string width in units of the grid spacing.
     * @param   lev Current level.
     * @param   eta Current time eta.
     * @return String width.
     */
    double StringWidth(const int lev, const double eta) {
        return 1./(Mr(eta) * sim->GetDx(lev));
    }

    /** @brief Returns time at which a level will be introduced.
     * @param   lev Level.
     * @return Refinement time for level lev.
     */
    double RefinementTime(const int lev) {
        return sim->GetDimN(lev) / (sqrt(2.*lambda)
                * string_width_threshold * sim->GetL());
    }

    /** @brief  Computes scale separation.
     * @param   eta Current time eta.
     * @return  log(m_r/H)
     */
    double Log(const double eta) {
        if (eta <= 0)
            return -DBL_MAX;
        return std::log( Mr(eta) / H(eta) );
    }

    /** @brief  Computes the physical box length.
     * @param   L       Comoving box length.
     * @param   eta     Current time eta.
     * @param   T1      T_1 temperature.
     * @param   mpl     Planck mass.
     * @param   gStar   Degrees of freedom g_*.
     * @return Physical box length.
     */
    double BoxToPhysical(const double L, const double eta, const double T1,
                         const double mpl, const double gStar) {
        return L * eta / Hubble(T1, mpl, gStar);
    }

    /** @brief Computes the physical Hubble parameter.
     * @param   T       Temperature.
     * @param   mpl     Planck mass.
     * @param   gStar   Degrees of freedom g_*.
     * @return Physical Hubble parameter.
     */
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

    long GetStringTags(const int lev);
    bool WriteXi(double time, std::string prefix);

    /** @brief Minimum allowed string width.
     */
    double string_width_threshold;

    /** @brief Minimum log for which we compute spectra.
     */
    double spectra_log_min = 5;

    /** @brief At what interval we want to write string length \xi.
     */
    double interval_xi_log = 0;

    /** @brief lambda parameter in Lagrangian.
     */
    const double lambda = 1;

    /** @brief Pointer to the simulation.
     */
    sledgehamr::Sledgehamr* sim;
};

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_COSMOLOGY_H_
