#ifndef SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
#define SLEDGEHAMR_GRAVITATIONAL_WAVES_H_

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

    int idx_offset;
    static constexpr int NScalars = 12;

    enum Gw { 
        u_xx = 0, u_yy, u_zz, u_xy, u_xz, u_yz, du_xx, du_yy, du_zz, du_xy,
        du_xz, du_yz, NGwScalars
    };
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
