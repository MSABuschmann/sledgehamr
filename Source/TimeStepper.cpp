#include <TimeStepper.H>

TimeStepper::TimeStepper (SledgeHAMR * owner)
{
	sim = owner;
}

void TimeStepper::Advance (int lev)
{
	// schedule regrids
	/* TODO */	

	PreAdvanceMessage(lev);

	// integrator->Advance(lev);
	/* TODO */

	PostAdvanceMessage(lev);

	if( lev != sim->finest_level )
		Advance(lev+1);

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
                if( sim->shadow_hierarchy /* && Regrid scheduled TODO */ )
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

        if( lev > 0 && sim->shadow_hierarchy /*&& Regrid scheduled TODO */ )
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
	amrex::Print() << "[Level " << lev << "] Advanced." << std::endl;
}
