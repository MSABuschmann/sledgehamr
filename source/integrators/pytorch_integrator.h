#ifndef SLEDGEHAMR_INTEGRATOR_PYTORCH_H_
#define SLEDGEHAMR_INTEGRATOR_PYTORCH_H_

#include <torch/script.h>

#include "integrator.h"

namespace sledgehamr {

/** @brief Implementation of the low-storage strong stability preserving third
 *         order Runge-Kutta integration scheme (PYTORCH).
 */
class IntegratorPytorch : public Integrator {
  public:
    IntegratorPytorch(Sledgehamr *owner);
    virtual void Integrate(LevelData &mf_old, LevelData &mf_new, const int lev,
                           const double dt, const double dx) override;

  private:
    void ParseParams();
    void LoadPytorchModel();

    std::string pytorch_model_filename;
    torch::jit::script::Module module;
    torch::Dtype dtype0 = torch::kFloat64;
    torch::TensorOptions tensoropt;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_PYTORCH_H_
