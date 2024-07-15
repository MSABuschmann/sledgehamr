#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SETUP_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SETUP_H_

#include <sledgehamr.h>

namespace FirstOrderPhaseTransition {

SLEDGEHAMR_ADD_SCALARS(Phi)
SLEDGEHAMR_ADD_CONJUGATE_MOMENTA(dPhi)

enum PotentialType { PureLambdaBar = 0, Piecewise = 1 };

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_SETUP_H_
