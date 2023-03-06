#include <TimeStepper.H>

TimeStepper::TimeStepper (SledgeHAMR * owner)
{
	sim = owner;
}

void TimeStepper::Advance (int lev)
{
	/* TODO */	
}

void TimeStepper::SynchronizeTimes ()
{
	for(int lev=1; lev <= sim->finest_level; ++lev)
		sim->grid_new[lev].t = sim->grid_new[0].t;
}
