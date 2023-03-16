#ifndef LevelData_H_
#define LevelData_H_

#include <AMReX_MultiFab.H>

/** @brief Class that holds the MultiFab data while 
 *  	   keeping track of time
 */
class LevelData
    : public amrex::MultiFab
{
public:
	LevelData () : amrex::MultiFab {} {};

	LevelData (amrex::BoxArray ba, amrex::DistributionMapping dm, int ncomp, int nghost, double time=0)
		: amrex::MultiFab {ba, dm, ncomp, nghost}
		, t{time} {};

	using amrex::MultiFab::MultiFab;

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
	}

	using amrex::MultiFab::define;

	/** @brief Static method that returns vector of times from a given LevelData vector.
	 * @param	mfs	Vector of pointers to MultiFab objects. 
	 *			MultiFab needs to be castable to LevelData. 
	 * @return 		Vector of times. Needs to be amrex::Vector to be
	 *			usable with e.g. amrex::FillPatchTwoLevels.
	 */
	static amrex::Vector<double> getTimes(std::vector<amrex::MultiFab*>& mfs)
	{
		amrex::Vector<double> times;

		for(amrex::MultiFab* mf : mfs){
			times.push_back( static_cast<LevelData*>(mf)->t );
		}

		return times;
	};

	/** @brief Time corresponding to amrex::MultiFab data.
	 */
	double t = 0.;

	/** @brief Flag as to whether the MutliFab has
	 *	   truncation errors encoded into it.
	 */
	bool contains_truncation_errors = false;
};

#endif // LevelData_H_
