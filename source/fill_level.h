#ifndef SLEDGEHAMR_FILL_LEVEL_H_
#define SLEDGEHAMR_FILL_LEVEL_H_

#include "sledgehamr.h"

namespace sledgehamr {

class FillLevel {
  public:
    FillLevel(Sledgehamr* owner, int level) : sim(owner), lev(level) {};

    void FromInitialStateFile();
    void FromCheckpointFile(std::string folder);
    void FromHdf5File(std::string initial_state_file);
    void FromArray(const int comp, double* data, const long long dimN);
    void FromArrayChunks(const int comp, double* data);
    void FromConst(const int comp, const double c);

  private:
    Sledgehamr* sim;
    const int lev;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_FILL_LEVEL_H_
