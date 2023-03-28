#ifndef AxionStrings_H_
#define AxionStrings_H_

#include <SledgeHAMR.H>
#include <SledgeHAMR_Utils.H>

/** @brief Class to simulate axion strings.
 */
class AxionStrings
    : public SledgeHAMR
{
public:

	ADD_SCALARS(Psi1, Psi2, Pi1, Pi2);

	AxionStrings ();

protected:

	/** @brief Function that calculates the RHS of the EOM at a single cell.
	 * @param	i		i-th cell index.
	 * @param	j		j-th cell index.
	 * @param	k		k-th cell index.
	 * @param	time		Current time.
	 * @param	lev		Currentl level.
	 * @param	rhs_fab		Container to be filled with RHS.
	 * @param	state_fab	Data from which to calculate RHS.
	 */
	void RHS (const int i, const int j, const int k, const double time, const int lev,
	 		    amrex::Array4<double> const& rhs_fab, 
			    amrex::Array4<double const> const& state_fab);

	/** @brief Function that tags individual cells for refinement.
	 * 	   TODO: Make virtual.
	 * @param	i		i-th cell index.
	 * @param	j		j-th cell index.
	 * @param	k		k-th cell index.
	 * @param	time		Current time.
	 * @param	lev		Currentl level.
	 * @param	state_fab	Data.
	 * @return	Boolean value as to whether cell should be refined or not.
	 */
	bool TagCellForRefinement (const int i, const int j, const int k, const double time, const int lev,
				   	     amrex::Array4<double const> const& state_fab);
	
private:

	/** @brief Computes the winding factor along a given axis. Will be non-zero
	 *	   if plaquette is pierced by a string.
	 * @param	i		i-th cell index.
	 * @param	j		j-th cell index.
	 * @param	k		k-th cell index.
	 * @param	state_fab	Data.
	 * @return 	Winding factor.
	 */
	int WindingAxis1 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
	int WindingAxis2 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
	int WindingAxis3 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
	
	/** @brief Checks for zero-crossings between two points in the complex scalar field.
	 * @param	Psi1_1	\Psi_1 of 1. point.
	 * @param	Psi2_1	\Psi_2 of 1. point.
	 * @param	Psi1_2	\Psi_1 of 2. point.
	 * @param	Psi2_2	\Psi_2 of 2. point.
	 * @return Sign of slope of zero-crossing. 0 if no crossing.
	 */
	int ZeroXing (double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2);
};

#endif //AxionStrings_H_
