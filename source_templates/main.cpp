#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include "sledgehamr.h"
#include "sledgehamr_init.h"
#include "projects.h"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // Timer for profiling.
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const double strt_total = amrex::ParallelDescriptor::second();

    // Start sledgeHAMR.
    sledgehamr::SledgehamrInit init;
//    sledgehamr::Sledgehamr *sledge = init.CreateInstance();
//    sledge->Init();
//    sledge->Evolve();
    axion_strings::axion_strings sledge;
    sledge.Init();
    sledge.Evolve();

    // Print runtime.
    double end_total = amrex::ParallelDescriptor::second() - strt_total;
    amrex::ParallelDescriptor::ReduceRealMax(end_total,
            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << std::endl << "Total Run Time: " << end_total << "s"
                   << std::endl;

    // Destroy timer for profiling.
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
