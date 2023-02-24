#ifndef LevelSynchronizer_H_
#define LevelSynchronizer_H_

#include <AMReX_BCUtil.H>
#include <AMReX_Interpolater.H>

#include <SledgeHAMR.H>
#include <LevelData.H>

class SledgeHAMR;

/** @brief This class handles all operations between two levels
 *         such as averaging down, interpolation to fine, filling
 *         of ghost cells, etc.
 *	   Class is friend of SledgeHAMR.
 */
class LevelSynchronizer
{
public:
	LevelSynchronizer (SledgeHAMR * owner);

	/** @brief Fills LeveData with information from a coarse level.
         *         This is used e.g. when a new level of refinement is added.
         * @param	lev	New level to be filled with data from lev=1.
	 * @param	time	Time of new level.
	 * @param	levdat	New level data.
	 */
	void FillCoarsePatch (int lev, double time, LevelData& levdat);

private:	

	/** @brief Fetches level data at a given level and time.
	 * @param	lev	Level at which data is to be fetched.
	 * @param	time	Time at which data is to be fetched. If time 
	 * 			does not align with t_old or t_new both states
	 *			will be returned.
	 * @param	levdat	Vector to be filled with fetched data. 
	 */
	void GetLevelData (int lev, double time, std::vector<LevelData*>& levdat);

	/** @brief Integer array containing the type of 
	 *         boundary condition at each boundary edge
	 */
	amrex::Vector<amrex::BCRec> bcs;

	/** @brief Pointer to AMReX interpolator to be
  	 * 	   used between levels 
	 */
	amrex::Interpolater * mapper = nullptr;

	/** @brief Pointer to owner on whose data this class operates
	 */
	SledgeHAMR * sim;		
};

#endif // LevelSynchronizer_H_
