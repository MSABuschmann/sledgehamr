#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "sledgehamr.h"
#include "sledgehamr_init.h"

int main(int argc, char *argv[]) {
    // Fix needed to accomodate change in amrex version 23.08. Will need long
    // term solution for this.
    amrex::ParmParse pp("amrex");
    pp.add("the_arena_is_managed", 1);
    pp.add("use_gpu_aware_mpi", 0);

    amrex::Initialize(argc, argv);

    // Timer for profiling.
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const double strt_total = amrex::ParallelDescriptor::second();

    // Start sledgeHAMR.
    sledgehamr::SledgehamrInit init;
    sledgehamr::Sledgehamr *sledge = init.CreateInstance();
    sledge->InitSledgehamr();
    sledge->Evolve();

    // Print runtime.
    double end_total = amrex::ParallelDescriptor::second() - strt_total;
    amrex::ParallelDescriptor::ReduceRealMax(
        end_total, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << std::endl
                   << "Total Run Time: " << end_total << "s" << std::endl;

    // Destroy timer for profiling.
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
