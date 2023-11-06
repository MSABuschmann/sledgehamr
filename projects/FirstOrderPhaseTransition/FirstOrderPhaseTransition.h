#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>

#include "setup.h"
#include "bubbles.h"
#include "kernels_rhs.h"
#include "kernels_tagging.h"
#include "kernels_misc.h"

namespace FirstOrderPhaseTransition {

FINISH_SLEDGEHAMR_SETUP

/** @brief Class to simulate a first order phase transition.
 */
class FirstOrderPhaseTransition : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(FirstOrderPhaseTransition)

    void Init() override;
    void SetParamsRhs(std::vector<double>& params, const double time,
                      const int lev) override;
    void BeforeTimestep(const double time) override;

  private:
    void ParseVariables();
    void ParseBubbles();
    void ComputeParameters();
    void SetProjections();

    void InjectBubbles(const double time);
    void InjectBubbleLevels(std::vector<int> ab);
    void FillBubbleLayout(const int lev, std::vector<int> ab);
    void AddBubbleValues(std::vector<int> ab);

    double lambda_bar;
    double quadratic, cubic, quartic;

    std::vector<Bubble> bubbles;
    int next_bubble = 0;
    int idx_perfmon_add_bubbles;
    std::vector<int> bubbles_to_inject;
};

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
