#ifndef SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_
#define SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Writes an AMReX plotfile for yt compatability.
 */
class AmrexPlotFile {
  public:
    /** @brief Constructior
     * @param   owner   Pointer to simulation.
     * @param   prefix  Local output folder.
     */
    AmrexPlotFile(Sledgehamr* owner, std::string prefix)
      : sim(owner), folder(prefix) {};

    void Write();

  private:
    /** @brief Local output folder.
     */
    std::string folder;

    /** @brief Pointer to simulation.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_
