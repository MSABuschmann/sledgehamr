#ifndef LevelSynchronizer_H_
#define LevelSynchronizer_H_

#include <LevelData.H>

/** @brief This class handles all operations between two levels
 *         such as averaging down, interpolation to fine, filling
 *         of ghost cells, etc.
 *	   Class is friend of SledgeHAMR.
 */
class LevelSynchronizer
{
public:

	/** @brief Fills LeveData with information from a coarse level.
         *         This is used e.g. when a new level of refinement is added.
         * @param	lev	New level to be filled with data from lev=1.
	 * @param	time	Time of new level.
	 * @param	levdat	New level data.
	 */
	void FillCoarsePatch (int lev, double time, LevelData& levdat) {};
};

#endif // LevelSynchronizer_H_
