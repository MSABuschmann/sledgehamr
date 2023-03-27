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

	ParseInputScalars();

	/* TODO: Check for checkpoint file etc. */
	InitFromScratch( t_start );
}

void SledgeHAMR::Evolve ()
{
	// Main loop over time.
	while( grid_new[0].t < t_end ){
		// Advance all levels starting at lev=0.
		// This performs an entire shadow/coarse level time step.	
		amrex::Print() << std::endl;
		time_stepper->Advance(0);

		// Write any output if requested.
		amrex::Print() << std::endl;
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
	grid_new[lev].define(ba, dm, ncomp, nghost, time);
	grid_old[lev].define(ba, dm, ncomp, nghost);

	SetBoxArray(lev, ba);
	SetDistributionMap(lev, dm);
	
	// If shadow hierarchy is used, above level is the shadow level
	// We now need to also make the coarse level.
	if( shadow_hierarchy && lev == 0 ){
		++lev;
		++finest_level;

		// create local copy of ba since ba is const
		amrex::BoxArray rba = ba;
		rba.refine(2);

		grid_new[lev].define(rba, dm, ncomp, nghost, time);
		grid_old[lev].define(rba, dm, ncomp, nghost);

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
	// Skip regrid right at the beginning of the sim
	// Allowed to be overwritten if no truncation errors
	// are used (TODO).
	if( time == t_start )
		return;

	// Current state.
	const amrex::MultiFab& state = grid_new[lev];

	// State containing truncation errors
	// if they have been calculated.
	const amrex::MultiFab& state_te = grid_old[lev];

	// Count number of tags.
	int ntags_total = 0;
	int ntags_user = 0;
	std::vector<int> ntags_trunc(scalar_fields.size(), 0);

	// Loop over boxes and cells.
	#pragma omp parallel reduction(+: ntags_total) reduction(+: ntags_user) reduction(vec_int_plus : ntags_trunc) 
	for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
		const amrex::Box& tilebox  = mfi.tilebox();
		const amrex::Array4<double const>& state_fab    = state.array(mfi);
		const amrex::Array4<double const>& state_fab_te = state_te.array(mfi);
	 	const amrex::Array4<char>& tag_arr = tags.array(mfi);

		// Tag with or without truncation errors.
		if ( shadow_hierarchy ) {
			ErrorEstWithTE(state_fab, state_fab_te, tag_arr, tilebox, 
					time, lev, &ntags_total, &ntags_user, &(ntags_trunc[0]));
		}else{
			ErrorEstWithoutTE(state_fab, state_fab_te, tag_arr, tilebox, 
				 	   time, lev, &ntags_total);
		}
	}

	// Collect all tags across MPI ranks.
	amrex::ParallelDescriptor::ReduceIntSum(ntags_total, 0);

	if( shadow_hierarchy ){
		amrex::ParallelDescriptor::ReduceIntSum(ntags_user, 0);
		amrex::ParallelDescriptor::ReduceIntSum(&(ntags_trunc[0]), ntags_trunc.size(), 0);
	}

	// Print statistics.
	long ncells = CountCells(lev);
	double ftotal = (double)ntags_total/(double)ncells;
	double fuser  = (double)ntags_user /(double)ncells;
	amrex::Print()  << "  Tagged cells at level " << lev << ": " << ntags_total << " of " << ncells 
			<< " (" << ftotal*100. << "\%)" << std::endl;
	if( shadow_hierarchy ){
		amrex::Print() << "    User-defined tags: " << ntags_user << std::endl;

		for(int i=0; i<scalar_fields.size(); ++i){
			amrex::Print()  << "    Truncation error tags on " << scalar_fields[i]->name 
					<< ": " << ntags_trunc[i] << std::endl;
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

void SledgeHAMR::ParseInputScalars ()
{
	te_crit.resize( scalar_fields.size() );
	double te_crit_default = 1e99;

	amrex::ParmParse pp("amr");
	pp.query("te_crit", te_crit_default);
	for(int n=0; n<scalar_fields.size(); ++n){
		te_crit[n] = te_crit_default;
		std::string ident = "te_crit_" + scalar_fields[n]->name;
		pp.query(ident.c_str(), te_crit[n]);
	}
}
