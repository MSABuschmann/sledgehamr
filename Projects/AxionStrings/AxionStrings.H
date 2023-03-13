#ifndef AxionStrings_H_
#define AxionStrings_H_

#include <SledgeHAMR.H>

/** @brief Class to simulate axion strings.
 */
class AxionStrings
    : public SledgeHAMR
{
public:

	/* Add all four scalar components 	
	 */
	ADD_SCALARS(Psi1, Psi2, Pi1, Pi2);

	AxionStrings ();

protected:

	/** @brief Inline function that calculates the RHS of the EOM at a single cell.
	 *	   Function body should be defined within this class definition.
	 * @param	i		i-th cell.
	 * @param	j		j-th cell.
	 * @param	k		k-th cell.
	 * @param	rhs_fab		Container to be filled with RHS.
	 * @param	state_fab	Data from which to calculate RHS.
	 */
	inline void RHS (const int i, const int j, const int k, const double time,
	 		    amrex::Array4<double> const& rhs_fab, 
			    amrex::Array4<double const> const& state_fab)
	{
		// Fetch field values.
		double Psi1 = state_fab(i, j, k, Scalar::Psi1);
		double Psi2 = state_fab(i, j, k, Scalar::Psi2);
		double Pi1  = state_fab(i, j, k, Scalar::Pi1);
		double Pi2  = state_fab(i, j, k, Scalar::Pi2);

		double eta = time;

		// Compute EOM
		double cross_term = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;
		
		rhs_fab(i, j, k, Scalar::Psi1) =  Pi1;
		rhs_fab(i, j, k, Scalar::Psi2) =  Pi2;
		rhs_fab(i, j, k, Scalar::Pi1)  = -Pi1*2./eta - Psi1 * cross_term;
		rhs_fab(i, j, k, Scalar::Pi2)  = -Pi2*2./eta - Psi2 * cross_term;	
	};
};

#endif //AxionStrings_H_
