#ifndef TimeStepper_H_
#define TimeStepper_H_

#include <SledgeHAMR.H>

class SledgeHAMR;

/** @brief Class that takes care of the sub-cycling in time algorithm
 *	   as well as the scheduling regrids.
 */
class TimeStepper
{
public:

	TimeStepper (SledgeHAMR * owner);	

	/** @brief Recursive function and the core of the
	 *	   sub-cycling time algorithm. Advances a given 
 	 *	   level by dt[level].
	 * @param	lev	Level to be advanced.
	 */
	void Advance (int lev);

	/** @brief Synchronizes the time of all levels.
	 */
	void SynchronizeTimes ();
private:
	
	/** @brief Pointer to owner on whose data this class operates.
	 */
	SledgeHAMR * sim;
};

#endif // TimeStepper_H
