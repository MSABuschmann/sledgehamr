#ifndef GRAVITATIONAL_WAVES_H_
#define GRAVITATIONAL_WAVES_H_

#include "sledgehamr.h"

namespace sledgehamr {

class GravitationalWaves  {
  public:
    GravitationalWaves(Sledgehamr* owner);

  private:
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // GRAVITATIONAL_WAVES_H_
