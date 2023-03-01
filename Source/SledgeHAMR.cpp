#include <AMReX_ParmParse.H>

#include <SledgeHAMR.H>

SledgeHAMR::SledgeHAMR ()
{
	amrex::Print() << "Starting sledgeHAMR..." << std::endl;

	level_synchronizer = new LevelSynchronizer(this);

	ParseInput();
	
	grid_new.resize(max_level+1);
	grid_old.resize(max_level+1);
}

SledgeHAMR::~SledgeHAMR ()
{
	delete level_synchronizer;
}

void SledgeHAMR::Init ()
{
	/* TODO: Check for checkpoint file etc. */
	InitFromScratch( t_start );
}

void SledgeHAMR::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					     const amrex::DistributionMapping& dm)
{
	const int ncomp = scalar_fields.size();

	// Define lowest level from scratch
	grid_new[lev].define(ba, dm, ncomp, nghost, time);
	grid_old[lev].define(ba, dm, ncomp, nghost);

	// If shadow hierarchy is used, above level is the shadow level
	// We now need to also make the coarse level.
	if( shadow_hierarchy ){
		++lev;

		// create local copy of ba since ba is const
		amrex::BoxArray rba;
		rba.refine(2);
		grid_new[lev].define(rba, dm, ncomp, nghost, time);
		grid_old[lev].define(rba, dm, ncomp, nghost);

		// already set for shadow level
		SetBoxArray(lev, rba);
		SetDistributionMap(lev, dm);
	}

	// Fill current level lev with initial state data
	/* TODO */

	// fill shadow level with data from coarse level
	if( shadow_hierarchy )
		level_synchronizer->AverageDownTo(0);
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

void SledgeHAMR::ClearLevel (int lev)
{
	grid_new[lev].clear();
	grid_old[lev].clear();
}

void SledgeHAMR::ParseInput ()
{
	{
		amrex::ParmParse pp("amr");
		pp.query("nghost", nghost);
	}
	
	{
		amrex::ParmParse pp("sim");
		pp.get("t_start", t_start);
		pp.get("t_end", t_end);
	}
}
