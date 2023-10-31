/* This file is auto-generated. MODIFICATIONS ARE FUTILE! */

#ifndef SLEDGEHAMR_PROJECTS_H_
#define SLEDGEHAMR_PROJECTS_H_

#include "axion_strings.h"
#include "axion_strings_postevolution.h"
#include "axion_strings_preevolution.h"
#include "first_order_phase_transition.h"

#define SLEDGEHAMR_PROJECT(str) {\
   if (str == "axion_strings") return new axion_strings::axion_strings;\
   if (str == "axion_strings_postevolution") return new axion_strings_postevolution::axion_strings_postevolution;\
   if (str == "axion_strings_preevolution") return new axion_strings_preevolution::axion_strings_preevolution;\
   if (str == "first_order_phase_transition") return new first_order_phase_transition::first_order_phase_transition;\
}

#endif // SLEDGEHAMR_PROJECTS_H_
