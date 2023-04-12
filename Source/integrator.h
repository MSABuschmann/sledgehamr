#ifndef Integrator_H_
#define Integrator_H_

#include <sledgehamr.h>

class SledgeHAMR;

/** @brief Pure virtual base class that handles the time
 *	   integration for a single level.
 */
class Integrator
{
public:
	Integrator (SledgeHAMR * owner);

	/** @brief Advances a single level.
	 * @param	lev	Level to be advanced.
	 */
	virtual void Advance (int lev);

private:

	/** @brief Pointer to owner on whose data this class operates.
	 */	
	SledgeHAMR * sim;
};

#endif // Integrator_H_
