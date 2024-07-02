#ifndef PROJECTS_AXION_STRINGS_PREEVOLUTION_H_
#define PROJECTS_AXION_STRINGS_PREEVOLUTION_H_

#include "../AxionStrings/cosmology.h"
#include "kernels_rhs.h"
#include <sledgehamr.h>

namespace AxionStringsPreevolution {

SLEDGEHAMR_FINISH_SETUP

/** @brief Class to simulate axion strings.
 */
class AxionStringsPreevolution : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(AxionStringsPreevolution)

    void SetParamsRhs(std::vector<double> &params, const double time,
                      const int lev) override;

    void Init() override;

    bool StopRunning(const double time) override;

    void SetRandomState();

  private:
    void ParseConstants();

    double log_0 = 2;
    double eta_0 = 2.3;
    double xi_0 = 0.18;
    double min_eta = 2;
    int random_state = -1;
    Cosmology cosmo;
};

}; // namespace AxionStringsPreevolution

#endif // PROJECTS_AXION_STRINGS_PREEVOLUTION_H_
