#ifndef SledgeHAMR_H_
#define SledgeHAMR_H_

#include <AMReX_AmrCore.H>

#include <LevelData.H>
#include <LevelSynchronizer.H>
#include <TimeStepper.H>
#include <Integrator.H>
#include <IOModule.H>
#include <Scalars.H>

class LevelSynchronizer;
class TimeStepper;
class Integrator;
class IOModule;

/** @brief Base class for all derived projects. Combines all the ingredients
 *         to make this code work. 
 */
class SledgeHAMR
    : public amrex::AmrCore
{
	friend class LevelSynchronizer;
	friend class TimeStepper;
	friend class Integrator;
	friend class IOModule;

public:

	/** @brief Creates instances of submodules and reads input parameters.
	 */
	SledgeHAMR ();

	virtual ~SledgeHAMR ();

	/** @brief Initalizes data from scratch or from checkpoint file. 
	 */
	void Init ();

	/** @brief Starts the evolution
	 */
	void Evolve ();

protected:

	/** @brief Make a new level from scratch using provided BoxArray and DistributionMapping.
	 *         Only used during initialization.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be created.
	 * @param	time	Time of new grid.
	 * @param	ba	New amrex::BoxArray.
	 * @param	dm	New amrex::DistributionMapping.
	 */
	virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
						  const amrex::DistributionMapping& dm) override;

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

	/** @brief Remake a new level using provided BoxArray and DistributionMapping,
	 *         and fills it with interpolated coarse level data.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be remade.
	 * @param	time	Time of new grid.
	 * @param	ba	New amrex::BoxArray.
	 * @param	dm	New amrex::DistributionMapping.
	 */
	virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
				      const amrex::DistributionMapping& dm) override;
	
	/** @brief Delete level data.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level to be deleted.
	 */
	virtual void ClearLevel (int lev) override;

	/** @brief Tag cells for refinement.
	 *         Overrides the pure virtual function in amrex::AmrCore.
	 * @param	lev	Level on which cells are tagged.
	 * @param	time	Time of said level.
	 * @param	ngrow	Grid growth factor	
	 */
	virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override { /* TODO */ };

	/** @brief Instance to perform operations between two levels.
	 */
	LevelSynchronizer * level_synchronizer;

	/** @brief Instance to perform the sub-cycling in time and
	 * 	   that handles regrid calls.
	 */
	TimeStepper * time_stepper;

	/** @brief Module that handles all IO operations
	 * 	   (with the exception of parsing the inputs file
	 */
	IOModule * io_module;

	/** @brief Holds the actual simulation data for all levels
	 *         at two different states in time.
	 */
	std::vector<LevelData> grid_new;
	std::vector<LevelData> grid_old;

	/** @brief Flag whether simulation should use a shadow hierarchy.
	 *         To be set in inputs file.
	 */
	bool shadow_hierarchy = false;

	/** @brief Holds pointers to all simulated scalar fields.
	 */
	std::vector<ScalarField*> scalar_fields;

	/** @brief Number of ghost cells.
	 */
	int nghost = 0;

	/** @brief Start and end times of the simulation.
	 */
	double t_start, t_end;

	/** @brief Time step size and grid spacing at each level.
	 */
	std::vector<double> dt, dx;

	/** @brief CFL criteria.
	 */
	double cfl;

	/** @brief Box length.
	 */
	double L;

	/** @brief Number of cells in each direction for each level.
	 */
	std::vector<unsigned int> dimN;

	/** @brief Number of coarse level cells in each direction. 
	 */
	int coarse_level_grid_size;

private:

	/** @brief Parse various input parameters from inputs file.
	 */
	void ParseInput ();
};

#endif //SledgeHAMR_H_
