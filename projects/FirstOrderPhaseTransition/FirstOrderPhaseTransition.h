#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>

#include "bubbles.h"
#include "kernels_misc.h"
#include "kernels_rhs.h"
#include "kernels_tagging.h"
#include "setup.h"

namespace FirstOrderPhaseTransition {

SLEDGEHAMR_FINISH_SETUP

/** @brief Class to simulate a first order phase transition.
 */
class FirstOrderPhaseTransition : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(FirstOrderPhaseTransition)

    void Init() override;
    void SetParamsRhs(std::vector<double> &params, const double time,
                      const int lev) override;
    void SetParamsGravitationalWaveRhs(std::vector<double> &params,
                                       const double time,
                                       const int lev) override;
    void SetParamsTruncationModifier(std::vector<double> &params,
                                     const double time, const int lev) override;
    void BeforeTimestep(const double time) override;

  private:
    void ParseVariables();
    void ParseBubbles();
    void ComputeParameters();
    void SetProjections();
    void AddSpectrumModification();

    std::vector<int> FindBubbles(const double time);
    void InjectBubbles(const double time);
    void InjectBubbleLevels(std::vector<int> ab);
    void FillBubbleLayout(const int lev, std::vector<int> ab);
    void AddBubbleValues(std::vector<int> ab);
    void MoveBubblesToCentre();

    bool GwSpectrum_UtimesK(double time, std::string prefix);
    bool GwSpectrum_2BubblesFrom1(double time, std::string prefix);

    int potential_type = 0;
    double lambda_bar = 0;
    double quadratic = 0, cubic = 0, quartic = 0;
    double vbar, vareps, phiesc;

    double tc = -1;
    double t0 = -1;

    std::vector<Bubble> bubbles;
    int next_bubble = 0;
    int idx_perfmon_add_bubbles;
    std::vector<int> bubbles_to_inject;

    std::vector<double> field_maxima;
    amrex::Vector<int> comp_vector;
    double maxima_time = -DBL_MAX;
};

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
