#include <SledgeHAMR.H>

SledgeHAMR::SledgeHAMR ()
{
	amrex::Print() << "Starting sledgeHAMR..." << std::endl;
}

SledgeHAMR::~SledgeHAMR ()
{}

void SledgeHAMR::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					     const amrex::DistributionMapping& dm)
{
	// Define a new level from scratch
	const int ncomp = grid_new[lev-1].nComp();
	const int nghost = grid_new[lev-1].nGrow();

	grid_new[lev].define(ba, dm, ncomp, nghost);
	grid_old[lev].define(ba, dm, ncomp, nghost);

	grid_new[lev].t = time;
	grid_old[lev].t = time - 1.e200;

	// Fill new level with coarse level data
	level_synchronizer.FillCoarsePatch(lev, time, grid_new[lev]);
}
