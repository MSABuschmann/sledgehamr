#ifndef PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_
#define PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_

#include <sledgehamr.h>
#include "../AxionStrings/cosmology.h"
#include "kernels_rhs.h"

namespace AxionStringsPostevolution {
using AxionStrings::GravitationalWavesRhs;
using AxionStrings::TruncationModifier;
using AxionStrings::TagCellForRefinement;

SLEDGEHAMR_FINISH_SETUP

/** @brief Class to simulate axion strings.
 */
class AxionStringsPostevolution : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(AxionStringsPostevolution)

    void SetParamsRhs(std::vector<double>& params, const double time,
                      const int lev) override;

    void Init() override;

    bool CreateLevelIf(const int lev, const double time) override {
        return cosmo.CreateLevelIf(lev, time);
    };

  private:
    void ParseConstants();
    void GetPreevolutionTime();
    void ReinterpretInitialState();

    double log_0 = 2;
    double eta_0 = 2.3;
    double eta_pre_0 = -1;
    double eta_transition = 2.8;
    double f_transition = 10;

    Cosmology cosmo;
};

}; // namespace AxionStringsPostevolution

#endif // PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_
