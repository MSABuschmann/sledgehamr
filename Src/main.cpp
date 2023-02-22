#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

int main(int argc, char* argv[])
{
	amrex::Initialize(argc,argv);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// wallclock time
	const double strt_total = amrex::ParallelDescriptor::second();



		
	// print runtime
	double end_total = amrex::ParallelDescriptor::second() - strt_total;
	amrex::ParallelDescriptor::ReduceRealMax(end_total, amrex::ParallelDescriptor::IOProcessorNumber());
	amrex::Print() << std::endl << "Total Run Time: " << end_total << "s" << std::endl;

	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
