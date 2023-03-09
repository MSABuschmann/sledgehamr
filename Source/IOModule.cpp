#include <IOModule.H>

IOModule::IOModule (SledgeHAMR * owner)
{
	sim = owner;

	// Determine and create output folder
	amrex::ParmParse pp("output");
	std::string output_folder;
	pp.get("output_folder", output_folder);
	amrex::UtilCreateDirectory(output_folder.c_str(), 0755);

	// Add various output formats
	double interval_slices = -1;
	pp.query("interval_slices", interval_slices);

	if( interval_slices >= 0 ){
		OutputModule out1(output_folder + "/slices", OUTPUT_FCT(IOModule::WriteSlices),
					interval_slices); 
		output.push_back( out1 );
	}
}

void IOModule::Write (bool force)
{
	for( OutputModule& out : output )
		out.Write(sim->grid_new[0].t, force);
}

void IOModule::FillLevelFromFile (int lev)
{
	// Figure out how files are organized.
	amrex::ParmParse pp("input");

	int N_chunks = 1;
	pp.query("N_chunks", N_chunks);

	if( N_chunks == 1 ){
		FillLevelFromFile_NoChunks(lev);
	}else if( N_chunks > 1 ){
		FillLevelFromFile_Chunks(lev);
	}
}

void IOModule::FillLevelFromFile_Chunks (int lev)
{
	/* TODO: Include version that supports files split in chunks */
}

void IOModule::FillLevelFromFile_NoChunks (int lev)
{
	amrex::ParmParse pp("input");
	std::string initial_state_file = "";
	pp.query("initial_state", initial_state_file);

	const unsigned int ncomp = sim->scalar_fields.size();
	const unsigned int dimN = sim->dimN[lev];
	const unsigned long long dsetsize = dimN*dimN*dimN;
	double* input_data = new double [dsetsize] ();

	// Iterate over fields but introduce offset such that
	// each node grabs a different file first.
	int lr = amrex::ParallelDescriptor::MyProc();
	int mr = amrex::ParallelDescriptor::NProcs();
	for(unsigned int f=0; f<ncomp; ++f){
		unsigned int f2 = (f+lr) % ncomp;

		// Determine which input file to use if any
		std::string initial_state_file_component = "";
		std::string query = "initial_state_" + sim->scalar_fields[f2]->name;
		pp.query(query.c_str(), initial_state_file_component);

		if( initial_state_file_component == "" )
			initial_state_file_component = initial_state_file;

		// If no file found, fill level with 0's. 
		// Otherwise read file and fill LevelData.
		if( initial_state_file_component == "" ){
			FillLevelFromConst(lev, f2, 0);
		}else{
			ReadFromHDF5(initial_state_file_component, 
					{sim->scalar_fields[f2]->name, "data"},
					input_data);
			FillLevelFromArray(lev, f2, input_data, dimN);
		}
	}

	delete[] input_data;	
}

void IOModule::FillLevelFromArray (int lev, const int comp, double* data, 
				       const unsigned long long dimN)
{
	LevelData& state = sim->grid_new[lev];
	
	#pragma omp parallel
	for ( amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		const auto& state_arr = state.array(mfi);

		amrex::ParallelFor(bx,
		 [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			unsigned long long ind =  (unsigned long long)i * dimN*dimN 
						+ (unsigned long long)j * dimN
						+ (unsigned long long)k;
			state_arr(i,j,k,comp) = data[ind];
		});
	}
}

void IOModule::FillLevelFromConst (int lev, const int comp, const double c)
{
	LevelData& state = sim->grid_new[lev];

	#pragma omp parallel
	for ( amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		const auto& state_arr = state.array(mfi);

		amrex::ParallelFor(bx,
		 [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			state_arr(i,j,k,comp) = c;
		});
	}
}

void IOModule::WriteSlices (double time, std::string prefix)
{
	amrex::Print() << "Write slices: " << prefix << " " << std::endl;
	/* TODO */
}
