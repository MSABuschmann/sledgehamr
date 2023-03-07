#include <Integrator.H>
#include <LevelData.H>

#include <AMReX_TimeIntegrator.H>

Integrator::Integrator (SledgeHAMR * owner)
{
	sim = owner;
}

void Integrator::Advance (int lev)
{
	std::swap(sim->grid_old[lev], sim->grid_new[lev]);
	LevelData& S_new = sim->grid_new[lev];

	// State with ghost cells as the integrator initial condition
	LevelData Sborder(sim->grids[lev], sim->dmap[lev], S_new.nComp(), 
			  	sim->nghost, sim->grid_old[lev].t);

	sim->level_synchronizer->FillPatch(lev, sim->grid_old[lev].t, Sborder);

	/* Integrate from (y,t) = (Sborder, time) by dt_lev to set S_new. */

	// Create integrator with the old state
	amrex::TimeIntegrator<amrex::MultiFab> integrator(Sborder);
	const auto geom_lev = sim->geom[lev];

	// Create a RHS source function we will integrate
	auto source_fun = [&](amrex::MultiFab& rhs, const amrex::MultiFab& state, const double time){
		//FillRHS(rhs, state, time, geom_lev, lev);
		/* TODO */
	};

	// Create a function to call after updating a state
	auto post_update_fun = [&](amrex::MultiFab& S_data, const double time) {
		// Fill ghost cells for S_data from interior & periodic BCs
		// and from interpolating coarser data in space/time at the current stage time.
		
		//sim->level_synchronizer->FillIntermediatePatch(lev, time, S_data);
		/* TODO */
	};

	integrator.set_rhs(source_fun);
	integrator.set_post_update(post_update_fun);
	integrator.advance(Sborder, S_new, sim->grid_old[lev].t, sim->dt[lev]);

	sim->grid_new[lev].t = sim->grid_old[lev].t + sim->dt[lev];
}
