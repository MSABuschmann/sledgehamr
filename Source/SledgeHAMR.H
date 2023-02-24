#ifndef SledgeHAMR_H_
#define SledgeHAMR_H_

#include <AMReX_AmrCore.H>

#include <LevelSynchronizer.H>

/** @brief Base class for all derived projects. Combines all the ingredients
 *         to make this code work. 
 */
class SledgeHAMR
    : public amrex::AmrCore
{
	friend class LevelSynchronizer;

public:

	SledgeHAMR ();
	virtual ~SledgeHAMR ();

protected:

	/** @brief Make a new level using provided BoxArray and DistributionMapping,
	 *         and fills it with interpolated coarse level data.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be created.
	 * @param	time	Time of new grid.
	 * @param	ba	New amrex::BoxArray.
	 * @param	dm	New amrex::DistributionMapping.
	 */
	virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
						 const amrex::DistributionMapping& dm) override;

	/** @brief Reake a new level using provided BoxArray and DistributionMapping,
	 *         and fills it with interpolated coarse level data.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be remade.
	 * @param	time	Time of new grid.
	 * @param	ba	New amrex::BoxArray.
	 * @param	dm	New amrex::DistributionMapping.
	 */
	virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
				      const amrex::DistributionMapping& dm) override {};
	
	/** @brief Delete level data.
	 * @param	lev	Level to be deleted.
	 */
	virtual void ClearLevel (int lev) override {};

	/** @brief Make a new level from scratch using provided BoxArray and DistributionMapping.
	 *         Only used during initialization.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be created.
	 * @param	time	Time of new grid.
	 * @param	ba	New amrex::BoxArray.
	 * @param	dm	New amrex::DistributionMapping.
	 */
	virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
						  const amrex::DistributionMapping& dm) override {};

	/** @brief Tag cells for refinement.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level on which cells are tagged.
	 * @param	time	Time of said level.
	 * @param	ngrow	Grid growth factor	
	 */
	virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override {};

	/** @brief Instance to perform operations between two levels
	 */
	LevelSynchronizer level_synchronizer;

	/** @brief Holds the actual simulation data for all levels
	 *         at two different states in time.
	 */
	std::vector<LevelData> grid_new;
	std::vector<LevelData> grid_old;
};

#endif //SledgeHAMR_H_l
