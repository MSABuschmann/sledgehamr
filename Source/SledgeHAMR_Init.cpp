#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <SledgeHAMR_Init.H>
#include <Projects.H>

SledgeHAMR_Init::SledgeHAMR_Init ()
{
	determine_project_id ();
	finish_AMReX_setup ();
}

SledgeHAMR* SledgeHAMR_Init::createInstance()
{
	SLEDGEHAMR_PROJECT(project_name);
	amrex::Print() << "Project not found!" << std::endl;
	return NULL;
}

void SledgeHAMR_Init::determine_project_id ()
{
	amrex::ParmParse pp("project");
	pp.get("name",project_name);
}

void SledgeHAMR_Init::finish_AMReX_setup ()
{
	int grid_size = 0;
	bool shadow_hierarchy = false;

	amrex::ParmParse pp("amr");
	pp.get("coarse_level_grid_size",grid_size);
	pp.get("shadow_hierarchy",shadow_hierarchy);

	if( shadow_hierarchy )
		grid_size /= 2;

	std::vector<int> vect(3,grid_size);
	pp.addarr("n_cell",vect);
}
