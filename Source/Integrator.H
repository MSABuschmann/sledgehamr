#ifndef Integrator_H_
#define Integrator_H_

#include <SledgeHAMR.H>

class SledgeHAMR;

class Integrator
{
public:
	Integrator (SledgeHAMR * owner);

	/** @brief Advances a single level.
	 * @param	lev	Level to be advanced.
	 */
	void Advance (int lev);

private:

	/** @brief Pointer to owner on whose data this class operates.
	 */	
	SledgeHAMR * sim;
};

#endif // Integrator_H_
