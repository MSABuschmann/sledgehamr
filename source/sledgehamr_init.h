#ifndef SLEDGEHAMR_SLEDGEHAMRINIT_H_
#define SLEDGEHAMR_SLEDGEHAMRINIT_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Determines with project has been requested and creates correct
 *         sledgehamr instance. Also feeds extra derived information to
 *         amrex::AmrCore to make sure it plays nicely.
 */
class SledgehamrInit {
  public:
    SledgehamrInit();
    Sledgehamr* CreateInstance();

  private:
    void DetermineProjectName();
    void FinishAMReXSetup();

    /** @brief Contains 'project.name' from inputs file.
     */
    std::string project_name;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMRINIT_H_
