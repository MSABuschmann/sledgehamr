#ifndef TimeStepper_H_
#define TimeStepper_H_

#include <SledgeHAMR.H>
#include <Integrator.H>

class SledgeHAMR;
class Integrator;

/** @brief Class that takes care of the sub-cycling in time algorithm
 *	   as well as the scheduling regrids.
 */
class TimeStepper
{
public:

	TimeStepper (SledgeHAMR * owner);	

	~TimeStepper ();

	/** @brief Recursive function and the core of the
	 *	   sub-cycling in time algorithm. Advances
 	 *	   a given level by dt[level].
	 * @param	lev	Level to be advanced.
	 */
	void Advance (int lev);

private:

	/** @brief Synchronizes two levels by averaging down.
	 *	   Computes truncation errors if regrid is
	 *	   scheduled.
	 * @param	lev	Synchronizes lev with either 
	 *			lev-1 or lev+1 depending on the 
	 *			regrid status.
	 */
	void SynchronizeLevels (int lev);	

	/** @brief Synchronizes the time of all levels.
	 *	   This is needed to avoid de-synchronization
	 *	   due to floating point precission errors
	 * 	   during t -> t + dt operations.
	 */
	void SynchronizeTimes ();

	/** @brief Prints message before and after each level
	 *	   has been advanced.
	 * @param	lev	Level that has/will be advanced.
	 */
	void PreAdvanceMessage (int lev);
	void PostAdvanceMessage (int lev);
	std::string LevelMessage (int lev);

	/** @brief Schedules a regrid ahead of such that we can
	 *	   compute truncations errors first.
	 * @param	lev	Level at which to tag cells.
	 */
	void ScheduleRegrid (int lev);

	/** @brief Performs a regrid if needed if no truncation
	 *	   tags are to be performed.
	 * @param	lev	Level at which to tag cells.
	 */
	void NoShadowRegrid (int lev);

	/** @brief Performs the actual regrid, either local or
	 *	   global as appropriate.
	 * @param	lev	Level at which to tag cells.
	 * @param	time	Current time.
	 */
	void DoRegrid (int lev, double time);

	/** @brief Vector of regridding intervals at each level.
	 */
	std::vector<double> regrid_dt;

	/** @brief Vector of times at which a given level has been
	 *  	   regridded last.
	 */
	std::vector<double> last_regrid_time;

	/** @brief Pointer to module that performs the actual 
	 *	   time integration.
	 */
	Integrator * integrator;

	/** @brief Pointer to owner on whose data this class operates.
	 */
	SledgeHAMR * sim;
};

#endif // TimeStepper_H
