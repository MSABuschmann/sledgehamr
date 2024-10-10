#ifndef PROJECTS_AXION_ONLY_H_
#define PROJECTS_AXION_ONLY_H_

#include <sledgehamr.h>

#include "kernels_energy_densities.h"
#include "kernels_rhs.h"
#include "kernels_tagging.h"

namespace AxionOnly {

SLEDGEHAMR_FINISH_SETUP

/** @brief Class to simulate axion strings.
 */
class AxionOnly : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(AxionOnly)

    /** @brief Override Init function and let cosmology module handle the setup.
     */
    void Init() override;
    void SetParamsRhs(std::vector<double> &params, const double time,
                      const int lev) override;

  private:
    void ParseVariables();
    void SetProjections();
    void SetSpectrum();

    double eta_c;
    double eta_star;
    double n;
    double N_QCD;
};

}; // namespace AxionOnly

#endif // PROJECTS_AXION_ONLY_H_
