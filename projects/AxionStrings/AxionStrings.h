#ifndef PROJECTS_AXION_STRINGS_H_
#define PROJECTS_AXION_STRINGS_H_

#include <sledgehamr.h>

#include "cosmology.h"
#include "kernels_rhs.h"

namespace AxionStrings {

FINISH_SLEDGEHAMR_SETUP

/** @brief Class to simulate axion strings.
 */
class AxionStrings : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(AxionStrings)

    void Init() override {
        cosmo.Init(this);
    };

    bool CreateLevelIf(const int lev, const double time) override {
        return cosmo.CreateLevelIf(lev, time);
    };

  private:
    Cosmology cosmo;
};

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_H_
