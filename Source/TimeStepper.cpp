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
}

TimeStepper::~TimeStepper ()
{
	delete integrator;
}

void TimeStepper::Advance (int lev)
{
	// schedule regrids
	/* TODO */	

	PreAdvanceMessage(lev);

	// Advance this level.
	integrator->Advance(lev);

	PostAdvanceMessage(lev);

	// Advance any finer levels twice.
	if( lev != sim->finest_level ){
		Advance(lev+1);
		Advance(lev+1);
	}

	SynchronizeLevels(lev);

	// regrid if needed
	/* TODO */

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
