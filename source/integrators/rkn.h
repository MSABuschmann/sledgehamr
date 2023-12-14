#ifndef SLEDGEHAMR_RKN_H_
#define SLEDGEHAMR_RKN_H_

#include <AMReX_IntegratorBase.H>

#include "integrator.h"

namespace sledgehamr {

/** @brief Implementation of the Runge-Kutta-Nystroem method of arbitrary order.
 */
class IntegratorRkn : public Integrator {
  public:
    IntegratorRkn(Sledgehamr* owner, const IntegratorType id);

  protected:
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) override;

  private:
    void SetButcherTableau();
    void ReadUserDefinedButcherTableau();

    /** @brief Number of nodes involved.
     */
    int number_nodes;

    /** @brief Butcher Tableau values.
     */
    std::vector<std::vector<double> > tableau;

    /** @brief b-weights.
     */
    std::vector<double> weights_b;

    /** @brief \bar{b} weights.
     */
    std::vector<double> weights_bar_b;

    /** @brief Node values.
     */
    std::vector<double> nodes;

    /** @brief Selected integrator type.
     */
    IntegratorType integrator_type;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_RKN_H_
