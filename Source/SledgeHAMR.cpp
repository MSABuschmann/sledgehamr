#include <AMReX_ParmParse.H>

#include <SledgeHAMR.H>

SledgeHAMR::SledgeHAMR ()
{
	amrex::Print() << "Starting sledgeHAMR..." << std::endl;

	// Initialize modules
	time_stepper = new TimeStepper(this);
	io_module = new IOModule(this);

	ParseInput();

	// Fill various level vectors	
	grid_new.resize(max_level+1);
	grid_old.resize(max_level+1);

	for(int lev=0; lev<=max_level; ++lev){
		dimN.push_back( coarse_level_grid_size * pow(2,lev-shadow_hierarchy) );
		dx.push_back( L/(double)dimN[lev] );
		dt.push_back( dx[lev] * cfl );
	}
}

SledgeHAMR::~SledgeHAMR ()
{
	delete level_synchronizer;
	delete time_stepper;
	delete io_module;
}

void SledgeHAMR::Init ()
{
	// Initialize here and not in the SledgeHAMR constructor such that it
	// knows about the number of scalar fields during construction.
	// Necessary so it can initialize boundary conditions.
	level_synchronizer = new LevelSynchronizer(this);

	/* TODO: Check for checkpoint file etc. */
	InitFromScratch( t_start );
}

void SledgeHAMR::Evolve ()
{
	while( grid_new[0].t < t_end ){
		amrex::Print() << std::endl;

		// Advance all levels starting at lev=0.
		// This performs an entire shadow/coarse level time step.	
		time_stepper->Advance(0);

		// Write any output if requested.
		io_module->Write();
	}
	
	// Force write at the end of simulation.	
	io_module->Write(true);
}

void SledgeHAMR::MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					     const amrex::DistributionMapping& dm)
{
	const int ncomp = scalar_fields.size();

	// Define lowest level from scratch
	amrex::AllPrint() << "define ld" << std::flush<<std::endl;
	grid_new[lev].define(ba, dm, ncomp, nghost, time);
	grid_old[lev].define(ba, dm, ncomp, nghost);
	amrex::AllPrint() << "done def ld" << std::flush <<std::endl;

	SetBoxArray(lev, ba);
	SetDistributionMap(lev, dm);
	
	// If shadow hierarchy is used, above level is the shadow level
	// We now need to also make the coarse level.
	if( shadow_hierarchy ){
		++lev;
		++finest_level;

		// create local copy of ba since ba is const
		amrex::BoxArray rba = ba;
		rba.refine(2);

	amrex::AllPrint() << "define ld" << std::flush<<std::endl;
		grid_new[lev].define(rba, dm, ncomp, nghost, time);
		grid_old[lev].define(rba, dm, ncomp, nghost);
	amrex::AllPrint() << "done def ld" << std::flush <<std::endl;

		// already set for shadow level
		SetBoxArray(lev, rba);
		SetDistributionMap(lev, dm);
	}

	// Fill current level lev with initial state data
	io_module->FillLevelFromFile(lev);

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

void SledgeHAMR::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
	// Current state.
	const amrex::MultiFab& state = grid_new[lev];

	// State containing truncation errors
	// if they have been calculated.
	const amrex::MultiFab& state_te = grid_old[lev];

	// Loop over boxes and cells.
	#pragma omp parallel
	for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
		const amrex::Box& tilebox  = mfi.tilebox();
		const amrex::Array4<double const>& state_fab    = state.array(mfi);
		const amrex::Array4<double const>& state_fab_te = state_te.array(mfi);
	 	const amrex::Array4<char>& tag_arr = tags.array(mfi);

		// Tag with or without truncation errors.
		if ( shadow_hierarchy ) {
			ErrorEstWithTE(state_fab, state_fab_te, tag_arr, tilebox, time, lev);
		}else{
			//ErrorEstWithoutTE(state_fab, tilebox, tags, time);
		}
	}
}

void SledgeHAMR::ParseInput ()
{
	{
		amrex::ParmParse pp("amr");
		pp.query("nghost", nghost);
		pp.query("shadow_hierarchy", shadow_hierarchy);
		pp.query("coarse_level_grid_size", coarse_level_grid_size);
	}
	
	{
		amrex::ParmParse pp("sim");
		pp.get("t_start", t_start);
		pp.get("t_end", t_end);
		pp.get("L", L);
		pp.get("cfl", cfl);
	}
}
