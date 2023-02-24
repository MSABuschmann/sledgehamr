#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>

#include <LevelSynchronizer.H>

LevelSynchronizer::LevelSynchronizer(SledgeHAMR * owner)
{
	sim = owner;
}

void LevelSynchronizer::FillCoarsePatch (int lev, double time, LevelData& levdat)
{
	BL_ASSERT(lev > 0);
	
	// get lev-1 data
	std::vector<LevelData*> cmf;
	GetLevelData(lev-1, time, cmf);
	    
	// boundary conditions
	amrex::CpuBndryFuncFab bndry_func(nullptr);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs, bndry_func);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev], bcs, bndry_func);
	
	// interpolate lev from lev-1 
	amrex::InterpFromCoarseLevel(levdat, time, *cmf[0], 0, 0, levdat.nComp(), sim->geom[lev-1], sim->geom[lev],
				     cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1), mapper, bcs, 0);                                         
}

void LevelSynchronizer::GetLevelData (int lev, double time, std::vector<LevelData*>& levdat)
{
	levdat.clear();
	LevelData * New = &sim->grid_new[lev]; 
	LevelData * Old = &sim->grid_old[lev]; 

	double teps = (New->t - Old->t)*1.e-3;

	if( time > New->t - teps && time < New->t + teps ){
		levdat.push_back( New );
	}else if( time > Old->t - teps && time < Old->t + teps ){
		levdat.push_back( Old );
	}else{
		levdat.push_back( Old );
		levdat.push_back( New );
	}
}
