#ifndef LevelData_H_
#define LevelData_H_

#include <AMReX_MultiFab.H>

/** @brief Struct that holds the MultiFab data while 
 *  	   keeping track of time
 */
class LevelData
    : public amrex::MultiFab
{
public:
	LevelData () {};

	LevelData (amrex::BoxArray ba, amrex::DistributionMapping dm, int ncomp, int nghost, double time=0)
		: amrex::MultiFab {ba, dm, ncomp, nghost}
		, t{time} {};

	/** @brief Overload amrex::MultiFab::define function to include time.
	 * @param	ba	Box array.
	 * @param	dm	Distribution mapping.
	 * @param	ncomp	Number of scalar components.
	 * @param	nghost	Number of ghost cells.
	 * @param	time 	t.
	 */
	void define(amrex::BoxArray ba, amrex::DistributionMapping dm, int ncomp, int nghost, double time)
	{
		define(ba, dm, ncomp, nghost);
		t = time;
	};

	using amrex::MultiFab::define;

	/** @brief Time of amrex::MultiFab data.
	 */
	double t = 0.;
};

#endif // LevelData_H_
