#include <TimeStepper.H>

TimeStepper::TimeStepper (SledgeHAMR * owner)
{
	sim = owner;

	// Initialize the correct integrator	
	amrex::ParmParse pp_inte("integration");

	int inte_type;
	pp_inte.get("type", inte_type);

	if( inte_type == 1 ){
		integrator = new Integrator(sim);
	}else if( inte_type == 5 ){
		// low-storage SSPRK3
		/* TODO */
		amrex::Abort("#error: Integration type not yet implemented");	
	}else{
		amrex::Abort("#error: Unknown integration type: " + std::to_string(inte_type));	
	}

	// Set regridding intervals.
	amrex::ParmParse pp_amr("amr");
	
	double reg_dt = 1e99;
	pp_amr.query("regrid_dt", reg_dt);

	for(int lev=0; lev<=sim->max_level; ++lev){
		regrid_dt.push_back( reg_dt * pow(2, lev - sim->shadow_hierarchy) );
		last_regrid_time.push_back( -1e99 );
	}
}

TimeStepper::~TimeStepper ()
{
	delete integrator;
}

void TimeStepper::Advance (int lev)
{
	// Perform or schedule regrids.
	if ( sim->shadow_hierarchy ) {
		// Schedule regrids ahead of time if needed such
		// that we can compute truncation errors in time.
		// In this case regrids will be performed at the
		// end of two time steps.
		ScheduleRegrid(lev);
	}else{
		// Invoke regridding routine at the beginning of 
		// a time step if we do not use a shadow hierarchy.
		NoShadowRegrid(lev);
	}

	PreAdvanceMessage(lev);

	// Advance this level.
	integrator->Advance(lev);

	PostAdvanceMessage(lev);

	// Advance any finer levels twice.
	if( lev != sim->finest_level ){
		Advance(lev+1);
		Advance(lev+1);
	}

	// Synchronize this level with finer/coarser levels.
	SynchronizeLevels(lev);

	// regrid if needed
	/* TODO */

	// Synchronize times to avoid any floating point
	// precision errors from advancing times on each 
	// level separately.
	if( lev == 0 )	
		SynchronizeTimes();
}

void TimeStepper::SynchronizeLevels (int lev)
{
        if( lev < sim->finest_level )
        {
                if( sim->shadow_hierarchy && false /* Regrid scheduled TODO */ )
                {
                        // a regrid is scheduled so we do not synchronize levels quite yet 
			// in order to compute truncation errors first.
                }else{
			// average lev+1 onto lev.
                        sim->level_synchronizer->AverageDownTo(lev);

			// since we averaged down we do not have truncation errors
			// available at lev+1.
               		sim->grid_old[lev+1].contains_truncation_errors = false; 
                }
        }

        if( lev > 0 && sim->shadow_hierarchy && false /* Regrid scheduled TODO */ )
        {
		// compute truncation errors for level lev and
		// average down between lev and lev-1.
                sim->level_synchronizer->ComputeTruncationErrors(lev);
               	sim->grid_old[lev].contains_truncation_errors = true; 
        }
}

void TimeStepper::SynchronizeTimes ()
{
	for(int lev=1; lev <= sim->finest_level; ++lev)
		sim->grid_new[lev].t = sim->grid_new[0].t;
}

void TimeStepper::PreAdvanceMessage (int lev)
{
	amrex::Print() << "[Level " << lev << "] Advancing ..." << std::endl;
}


void TimeStepper::PostAdvanceMessage (int lev)
{
	amrex::Print() << "[Level " << lev << "] Advanced to t=" << sim->grid_new[lev].t << "." << std::endl;
}

void TimeStepper::ScheduleRegrid (int lev)
{
	/* TODO */
}

void TimeStepper::NoShadowRegrid (int lev)
{
	/* TODO: Implement semi-static case */

        // Regrid changes level "lev+1" so we don't regrid on max_level
        if( lev >= sim->max_level )
		return;

	// Check if enough time since last regrid has passed.
	// We add dt[lev] since we do not want to violate this criteria
	// next time around in case we skip this regrid.
	double time = sim->grid_new[lev].t;
	if( time + sim->dt[lev] <= last_regrid_time[lev] + regrid_dt[lev] )
		return;

	// Check user requirement if we want to invoke a new level.
	// Pass it the level to be created and the time by which 
	// the next regrid could be performed if we were to skip
	// this regrid.
	if( !sim->CanCreateLevel(lev+1, time + sim->dt[lev]) )
		return;

	// Actually do regrid if we made it this far.
	DoRegrid(lev, time);
}

void TimeStepper::DoRegrid (int lev, double time)
{
	/* TODO */
	amrex::Print() << "Perform regrid at level " << lev << " and t=" << time << std::endl;

	for(int l=lev; l <= sim->finest_level; ++l){
		last_regrid_time[l] = time;
	}
}
