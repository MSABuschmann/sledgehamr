#ifndef AxionStrings_H_
#define AxionStrings_H_

#include <SledgeHAMR.H>

/** @brief Class to simulate axion strings.
 */
class AxionStrings
    : public SledgeHAMR
{
	ADD_SCALARS(Psi1, Psi2, Pi1, Pi2);

public:
	AxionStrings ();
};

#endif //AxionStrings_H_
