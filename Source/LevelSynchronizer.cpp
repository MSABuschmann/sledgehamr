#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

#include <LevelSynchronizer.H>

LevelSynchronizer::LevelSynchronizer(SledgeHAMR * owner)
{
	sim = owner;
	
	const int ncomp = sim->scalar_fields.size();
	
	// Boundary conditions
	bcs.resize(ncomp);
	for(int n = 0; n < ncomp; ++n){
		for(int i = 0; i < AMREX_SPACEDIM; ++i){
			bcs[n].setLo(i, amrex::BCType::int_dir);
			bcs[n].setHi(i, amrex::BCType::int_dir);
		}
	}

	// Set interpolation type between levels
	amrex::ParmParse pp("amr");
	int interpolation_type = InterpType::CellConservativeQuartic;
	pp.query("interpolation_type", interpolation_type);
	
	if( interpolation_type == InterpType::CellConservativeLinear ){
		mapper = &amrex::cell_cons_interp;
	}else if( interpolation_type == InterpType::CellConservativeQuartic ){
		mapper = &amrex::quartic_interp;
	}else if( interpolation_type == InterpType::PCInterp ){
		mapper = &amrex::pc_interp;
	}else{
		amrex::Error("Unsupported interpolation type");
	}
}

void LevelSynchronizer::FillCoarsePatch (int lev, double time, amrex::MultiFab& mf)
{
	BL_ASSERT(lev > 0);
	
	// get lev-1 data
	std::vector<amrex::MultiFab*> cmf = GetLevelData(lev-1, time);
	    
	// boundary conditions
	amrex::CpuBndryFuncFab bndry_func(nullptr);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs, bndry_func);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev  ], bcs, bndry_func);
	
	// interpolate lev from lev-1 
	amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, mf.nComp(), sim->geom[lev-1], sim->geom[lev],
				     cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1), mapper, bcs, 0);                                         
}

void LevelSynchronizer::FillPatch (int lev, double time, amrex::MultiFab& mf)
{
	// Get data and boundary conditions for level lev.
	amrex::Vector<amrex::MultiFab*> fmfs = GetLevelData(lev, time);
	amrex::Vector<double> ftime = LevelData::getTimes(fmfs);

	amrex::CpuBndryFuncFab bndry_func(nullptr);
	amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev], bcs, bndry_func);

	if( lev == 0 ){
		// Call FillPatchSingleLevel for the coarse level
		amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, 0, 0, mf.nComp(), sim->geom[lev], fphysbc, 0);
	}else{
		// Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1) level 
		amrex::Vector<amrex::MultiFab*> cmfs = GetLevelData(lev-1, time);
		amrex::Vector<double> ctime = LevelData::getTimes(cmfs);

		amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs, bndry_func);
	
		amrex::FillPatchTwoLevels(mf, time, cmfs, ctime, fmfs, ftime, 0, 0, mf.nComp(),
					  sim->geom[lev-1], sim->geom[lev], cphysbc, 0, fphysbc, 0, 
					  sim->refRatio(lev-1), mapper, bcs, 0);
	}
}

void LevelSynchronizer::FillIntermediatePatch (int lev, double time, amrex::MultiFab& mf)
{
//	FillPatch(lev, time, mf);
}

void LevelSynchronizer::AverageDownTo (int lev)
{
	amrex::average_down(sim->grid_new[lev+1], sim->grid_new[lev], sim->geom[lev+1], sim->geom[lev], 
			    0, sim->grid_new[lev].nComp(), sim->refRatio(lev));
}

void LevelSynchronizer::ComputeTruncationErrors (int lev)
{
	/* TODO */
}

amrex::Vector<amrex::MultiFab*> LevelSynchronizer::GetLevelData (int lev, double time)
{
	amrex::Vector<amrex::MultiFab*> mfs;	
	LevelData * New = &sim->grid_new[lev]; 
	LevelData * Old = &sim->grid_old[lev]; 

	double teps = fabs(New->t - Old->t)*1.e-3;

	// Add either new, old or both states
	if( time > New->t - teps && time < New->t + teps ){
		mfs.push_back( New );
	}else if( time > Old->t - teps && time < Old->t + teps ){
		mfs.push_back( Old );
	}else{
		mfs.push_back( Old );
		mfs.push_back( New );
	}

	return mfs;
}
