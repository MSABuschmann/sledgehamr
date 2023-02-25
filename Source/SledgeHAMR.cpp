#include <SledgeHAMR.H>

SledgeHAMR::SledgeHAMR ()
{
	amrex::Print() << "Starting sledgeHAMR..." << std::endl;
	level_synchronizer = new LevelSynchronizer(this);
}

SledgeHAMR::~SledgeHAMR ()
{
	delete level_synchronizer;
}

void SledgeHAMR::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					     const amrex::DistributionMapping& dm)
{
	const int ncomp = grid_new[lev-1].nComp();
	const int nghost = grid_new[lev-1].nGrow();

	// Define a new level from scratch
	grid_new[lev].define(ba, dm, ncomp, nghost, time);
	grid_old[lev].define(ba, dm, ncomp, nghost);

	// Fill new level with coarse level data
	level_synchronizer->FillCoarsePatch(lev, time, grid_new[lev]);
}

void SledgeHAMR::RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba, 
				  const amrex::DistributionMapping& dm)
{
	const int ncomp = grid_new[lev].nComp();
	const int nghost = grid_new[lev].nGrow();

	// remake new_grid and fill with data
	LevelData new_state(ba, dm, ncomp, nghost, grid_new[lev].t);
	level_synchronizer->FillPatch(lev, time, new_state);
	std::swap(new_state, grid_new[lev]);
	new_state.clear();

	// remake old_grid
	grid_old[lev].clear();
	grid_old[lev].define(ba, dm, ncomp, nghost);
}
