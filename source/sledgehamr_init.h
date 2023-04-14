#ifndef SLEDGEHAMR_SLEDGEHAMRINIT_H_
#define SLEDGEHAMR_SLEDGEHAMRINIT_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Determines with project has been requested.
 *         Also feeds extra derived information to amrex::AmrCore.
 */
class SledgehamrInit {
  public:
    SledgehamrInit();

    /** @brief Creates instance of requested derived project class.
     * @return Pointer to derived class instance cast to base. Returns NULL if
     *         none found.
     */
    Sledgehamr* CreateInstance();

  private:
    /** @brief Reads project name from inputs file.
     */
    void DetermineProjectName();

    /** @brief Set parameters such as 'amr.n_cell' for amrex::AmrCore.
     */
    void FinishAMReXSetup();

    /** @brief Contains 'project.name' from inputs file.
     */
    std::string project_name;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMRINIT_H_
