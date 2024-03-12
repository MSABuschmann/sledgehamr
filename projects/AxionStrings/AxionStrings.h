#ifndef PROJECTS_AXION_STRINGS_H_
#define PROJECTS_AXION_STRINGS_H_

#include <sledgehamr.h>

#include "cosmology.h"
#include "kernels_rhs.h"

namespace AxionStrings {

SLEDGEHAMR_FINISH_SETUP

/** @brief Class to simulate axion strings.
 */
class AxionStrings : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(AxionStrings)

    /** @brief Override Init function and let cosmology module handle the setup.
     */
    void Init() override {
        cosmo.Init(this);
    };

    /** @brief We want to create a level only when string width is less than
     *         threshold.
     * @param   lev     Level in question.
     * @param   time    Current time.
     * @return  Whether we want to create level or not.
     */
    bool CreateLevelIf(const int lev, const double time) override {
        return cosmo.CreateLevelIf(lev, time);
    };

  private:
    /** @brief Axion cosmology module to set up the scenario.
     */
    Cosmology cosmo;
};

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_H_
