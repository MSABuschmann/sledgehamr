#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>

#include <LevelSynchronizer.H>

LevelSynchronizer::LevelSynchronizer(SledgeHAMR * owner)
{
	sim = owner;
}

void LevelSynchronizer::FillCoarsePatch (int lev, double time, LevelData& ld)
{
	BL_ASSERT(lev > 0);
	
	// get lev-1 data
	std::vector<LevelData*> cld = GetLevelData(lev-1, time);
	    
	// boundary conditions
	amrex::CpuBndryFuncFab bndry_func(nullptr);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs, bndry_func);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev  ], bcs, bndry_func);
	
	// interpolate lev from lev-1 
	amrex::InterpFromCoarseLevel(ld, time, *cld[0], 0, 0, ld.nComp(), sim->geom[lev-1], sim->geom[lev],
				     cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1), mapper, bcs, 0);                                         
}

void LevelSynchronizer::FillPatch (int lev, double time, LevelData& ld)
{
	// Get data and boundary conditions for level lev.
	amrex::Vector<LevelData*> fld = GetLevelData(lev, time);
	amrex::Vector<double> ftime = LevelData::getTimes(fld);

	amrex::CpuBndryFuncFab bndry_func(nullptr);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev], bcs, bndry_func);

	if( lev == 0 ){
		// Call FillPatchSingleLevel for the coarse level
		amrex::FillPatchSingleLevel(ld, time, fld, ftime, 0, 0, ld.nComp(), sim->geom[lev], fphysbc, 0);
	}else{
		// Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1) level 
		amrex::Vector<LevelData*> cld = GetLevelData(lev-1, time);
		amrex::Vector<double> ctime = LevelData::getTimes(cld);

		amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs, bndry_func);
	
		amrex::FillPatchTwoLevels(ld, time, cld, ctime, fld, ftime, 0, 0, ld.nComp(),
					  sim->geom[lev-1], sim->geom[lev], cphysbc, 0, fphysbc, 0, 
					  sim->refRatio(lev-1), mapper, bcs, 0);
	}
}

amrex::Vector<LevelData*> LevelSynchronizer::GetLevelData (int lev, double time)
{
	amrex::Vector<LevelData*> ld;	
	LevelData * New = &sim->grid_new[lev]; 
	LevelData * Old = &sim->grid_old[lev]; 

	double teps = (New->t - Old->t)*1.e-3;

	// Add either new, old or both states
	if( time > New->t - teps && time < New->t + teps ){
		ld.push_back( New );
	}else if( time > Old->t - teps && time < Old->t + teps ){
		ld.push_back( Old );
	}else{
		ld.push_back( Old );
		ld.push_back( New );
	}

	return ld;
}
