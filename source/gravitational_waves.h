#ifndef GRAVITATIONAL_WAVES_H_
#define GRAVITATIONAL_WAVES_H_

#include "sledgehamr.h"

namespace sledgehamr {

class GravitationalWaves {
  public:
    GravitationalWaves(Sledgehamr* owner);

    void ComputeSpectrum(hid_t file_id);

  private:
    double IndexToK(int a, int N);
    double GetProjection(int i, int j, int abc[3], int N);
    double GetLambda(int i, int j, int l, int m, int abc[3], int N);

    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // GRAVITATIONAL_WAVES_H_
