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

#define DO_PRAGMA(x) _Pragma(#x)

/** @brief Macro to override the SledgeHAMR::FillRHS function within the project class.
 *	   Will autmatically be added by the ADD_SCALARS macro to the project. This is 
 *	   needed to ensure the RHS function can be inlined by the compiler.
 */
#define FILL_RHS virtual void FillRHS (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const double time,\
					       const amrex::Geometry& geom, int lev) override\
	{\
		DO_PRAGMA(omp parallel)\
		for ( amrex::MFIter mfi(rhs_mf, true); mfi.isValid(); ++mfi ){\
			const amrex::Box& bx = mfi.tilebox();\
			const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi);\
			const amrex::Array4<double const>& state_fab = state_mf.array(mfi);\
			amrex::ParallelFor(bx,\
			[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept\
			{\
				RHS(i, j, k, time, lev, rhs_fab, state_fab);\
			});\
		}\
	};

#define ERROR_EST_WITH_TE virtual void ErrorEstWithTE  (const amrex::Array4<double const>& state_fab,\
							const amrex::Array4<double const>& state_fab_te,\
							const amrex::Array4<char>& tagarr,\
							const amrex::Box& tilebox, double time, int lev,\
							int* ntags_user)\
	{\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_user)++;\
			}\
		});\
	};\

#define ERROR_EST_WITHOUT_TE virtual void ErrorEstWithoutTE  (const amrex::Array4<double const>& state_fab,\
							      const amrex::Array4<double const>& state_fab_te,\
							      const amrex::Array4<char>& tagarr,\
							      const amrex::Box& tilebox, double time, int lev,\
							      int* ntags_user)\
	{\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_user)++;\
			}\
		});\
	};\

/** @brief Base class for all derived projects. Combines all the ingredients
 *         to make this code work. 
 */
class SledgeHAMR
    : public amrex::AmrCore
{
	// Give submodules access to data.
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
	 * @param	lev		Level on which cells are tagged.
	 * @param	time		Time of said level.
	 * @param	ngrow		Grid growth factor.
	 * @param	ntags_user	Counts number of user-defined tags.
	 */
	virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

	virtual void ErrorEstWithTE (const amrex::Array4<double const>& state_fab,
				     const amrex::Array4<double const>& state_fab_te,
				     const amrex::Array4<char>& tagarr,
				     const amrex::Box& tilebox, double time, int lev,
				     int* ntags_user) = 0;
	
	virtual void ErrorEstWithoutTE (const amrex::Array4<double const>& state_fab,
				        const amrex::Array4<double const>& state_fab_te,
				        const amrex::Array4<char>& tagarr,
				        const amrex::Box& tilebox, double time, int lev,
				        int* ntags_user) = 0;


	/** @brief Virtual function that loops over a state to fill the RHS. Has to be defined by
	 *	   derived project class, though this will be hidden from the user. Reason for this
	 *	   is so that the RHS function, declared also by the project class, can be inlined 
	 *	   within this function, which it otherwise couldn't as it is not clear at compile time
	 *	   which project will be run.
	 * @param	rhs_mf		Empty MultiFab to be filled with RHS.
	 * @param	state_mf	State from which the RHS is to be computed.	
	 * @param	time		Current time.
	 * @param	geom		Geometry of the current level.
	 * @param	lev		Currently level.
	 */
	virtual void FillRHS (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const double time, 
					       const amrex::Geometry& geom, int lev) = 0;

	/** @brief Function to check whether a given level is allowed
	 *	   to be created by a given time. Can be overridden by
	 *	   the project class. Will always allow a level to be
	 *	   created by default.
	 * @param 	lev	Level to be created.
	 * @param	time	Time before which the level would be
	 *			created.
	 * @return	If the given level is allowed to be created.
	 */
	virtual bool CanCreateLevel (int lev, double time)
	{ return true; };

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
