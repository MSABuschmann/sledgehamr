/* This file is auto-generated. MODIFICATIONS ARE FUTILE! */

#ifndef SLEDGEHAMR_PROJECTS_H_
#define SLEDGEHAMR_PROJECTS_H_

#include "AxionStrings.h"
#include "AxionStringsPostevolution.h"
#include "AxionStringsPreevolution.h"
#include "FirstOrderPhaseTransition.h"
#include "MinimalExample.h"

#define SLEDGEHAMR_PROJECT(str) {\
   if (str == "AxionStrings") return new AxionStrings::AxionStrings;\
   if (str == "AxionStringsPostevolution") return new AxionStringsPostevolution::AxionStringsPostevolution;\
   if (str == "AxionStringsPreevolution") return new AxionStringsPreevolution::AxionStringsPreevolution;\
   if (str == "FirstOrderPhaseTransition") return new FirstOrderPhaseTransition::FirstOrderPhaseTransition;\
   if (str == "MinimalExample") return new MinimalExample::MinimalExample;\
}

#endif // SLEDGEHAMR_PROJECTS_H_
