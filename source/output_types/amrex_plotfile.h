#ifndef SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_
#define SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_

#include "sledgehamr.h"

namespace sledgehamr {

class AmrexPlotFile {
  public:
    AmrexPlotFile(Sledgehamr* owner, std::string prefix)
      : sim(owner), folder(prefix) {};

    void Write();

  private:
    std::string folder;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_AMREX_PLOTFILE_H_
