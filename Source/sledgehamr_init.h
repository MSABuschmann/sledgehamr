#ifndef SledgeHAMR_Init_H_
#define SledgeHAMR_Init_H_

#include<SledgeHAMR.H>

/** @brief Determines with project has been requested.
 *         Also feeds extra derived information to amrex::AmrCore.
 */
class SledgeHAMR_Init
{
public:
	SledgeHAMR_Init ();

	/** @brief Creates instance of requested derived project class.
 	 * @return Pointer to derived class instance cast to base. 
 	 *         Returns NULL if none found.
 	 */
	SledgeHAMR* CreateInstance();
	
private:
	/** @brief Reads project name from inputs file
         */
	void DetermineProjectName ();

	/** @brief Set parameters such as 'amr.n_cell' for amrex::AmrCore
         */
	void FinishAMReXSetup ();

	/** @brief Contains 'project.name' from inputs file
 	 */
	std::string project_name;
};

#endif // SledgeHAMR_Init_H_
