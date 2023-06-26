#ifndef SLEDGEHAMR_RKN_H_
#define SLEDGEHAMR_RKN_H_

#include <AMReX_IntegratorBase.H>

#include "integrator.h"

namespace sledgehamr {

class IntegratorRkn : public Integrator {
  public:
    IntegratorRkn(Sledgehamr* owner, const IntegratorType id);

  protected:
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) override;

  private:
    void SetButcherTableau();
    void ReadUserDefinedButcherTableau();

    int number_nodes;
    std::vector<std::vector<double> > tableau;
    std::vector<double> weights_b;
    std::vector<double> weights_bar_b;
    std::vector<double> nodes;

    IntegratorType integrator_type;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_RKN_H_
