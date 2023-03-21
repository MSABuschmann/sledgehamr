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
	 * @param	time		Current time.
	 * @param	lev		Currentl level.
	 * @param	rhs_fab		Container to be filled with RHS.
	 * @param	state_fab	Data from which to calculate RHS.
	 */
	inline void RHS (const int i, const int j, const int k, const double time, const int lev,
	 		    amrex::Array4<double> const& rhs_fab, 
			    amrex::Array4<double const> const& state_fab)
	{
		// Fetch field values.
		double Psi1 = state_fab(i, j, k, Scalar::Psi1);
		double Psi2 = state_fab(i, j, k, Scalar::Psi2);
		double Pi1  = state_fab(i, j, k, Scalar::Pi1);
		double Pi2  = state_fab(i, j, k, Scalar::Pi2);

		double eta = time;

		// Compute Laplacians.
		double dx2 = dx[lev] * dx[lev];
		constexpr unsigned int order = 1;
		double LaplacianPsi1 = SledgeHAMR_Utils::Laplacian<order>(state_fab, i, j, k, Scalar::Psi1, dx2);
		double LaplacianPsi2 = SledgeHAMR_Utils::Laplacian<order>(state_fab, i, j, k, Scalar::Psi2, dx2);

		// Compute EOM
		double cross_term = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;

		rhs_fab(i, j, k, Scalar::Psi1) =  Pi1;
		rhs_fab(i, j, k, Scalar::Psi2) =  Pi2;
		rhs_fab(i, j, k, Scalar::Pi1)  = -Pi1*2./eta + LaplacianPsi1 - Psi1 * cross_term;
		rhs_fab(i, j, k, Scalar::Pi2)  = -Pi2*2./eta + LaplacianPsi2 - Psi2 * cross_term;	

		// DEBUG TODO
		if( i==0 && j==0 && k==0 )
			amrex::Print() << Psi1 << " " << Psi2 << " " << Pi1 << " " << Pi2 << std::endl;
	};

	inline bool TagCellForRefinement(const int i, const int j, const int k, const double time, const int lev,
					  	    amrex::Array4<double const> const& state_fab)
	{
		// Check all three directions for string piercings.
 		if( winding_axis1(i, j, k, state_fab) != 0 )
			return true;

 		if( winding_axis2(i, j, k, state_fab) != 0 )
			return true;

 		if( winding_axis3(i, j, k, state_fab) != 0 )
			return true;

		return false;
	};

	inline int winding_axis1(int i, int j, int k, amrex::Array4<double const> const& state_fab)
	{
		return    winding(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
				  state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2))
			+ winding(state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2),
				  state_fab(i+1,j+1,k  ,Scalar::Psi1),state_fab(i+1,j+1,k  ,Scalar::Psi2))
			+ winding(state_fab(i+1,j+1,k  ,Scalar::Psi1),state_fab(i+1,j+1,k  ,Scalar::Psi2),
				  state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2))
			+ winding(state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2),
				  state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
	};

	inline int winding_axis2(int i, int j, int k, amrex::Array4<double const> const& state_fab)
	{
		return    winding(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
				  state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2))
			+ winding(state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2),
				  state_fab(i+1,j  ,k+1,Scalar::Psi1),state_fab(i+1,j  ,k+1,Scalar::Psi2))
			+ winding(state_fab(i+1,j  ,k+1,Scalar::Psi1),state_fab(i+1,j  ,k+1,Scalar::Psi2),
				  state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2))
			+ winding(state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2),
				  state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
	};

	inline int winding_axis3(int i, int j, int k, amrex::Array4<double const> const& state_fab)
	{
		return    winding(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
				  state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2))
			+ winding(state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2),
				  state_fab(i  ,j+1,k+1,Scalar::Psi1),state_fab(i  ,j+1,k+1,Scalar::Psi2))
			+ winding(state_fab(i  ,j+1,k+1,Scalar::Psi1),state_fab(i  ,j+1,k+1,Scalar::Psi2),
				  state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2))
			+ winding(state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2),
				  state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
	};

	inline int winding(double Psi1_1,double Psi2_1,double Psi1_2,double Psi2_2)
	{
		if( Psi2_1 * Psi2_2 < 0 )
		{
			if( Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0 )
				return 1;
			else
				return -1;
		}
		return 0;
	};
};

#endif //AxionStrings_H_
