#ifndef sledgeHAMR_H_
#define sledgeHAMR_H_

#include <AMReX_AmrCore.H>

class SledgeHAMR
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor
    //             
    SledgeHAMR ();
    virtual ~SledgeHAMR ();


protected:

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) override {};

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) override {};

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override {};

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) override {};

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override {};
};

#endif //sledgeHAMR_H_
