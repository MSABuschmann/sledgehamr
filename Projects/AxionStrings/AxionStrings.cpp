#include <AxionStrings.H>

AxionStrings::AxionStrings ()
{
	amrex::Print() << "Starting AxionStrings project..." << std::endl;
	amrex::Print() << std::endl;
}

inline void AxionStrings::RHS (const int i, const int j, const int k, const double time, const int lev,
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
}

inline bool AxionStrings::TagCellForRefinement(const int i, const int j, const int k, const double time, const int lev,
						    amrex::Array4<double const> const& state_fab)
{
	// Check all three plaquettes (in positive index direction) 
	// for string piercings.
	if( WindingAxis1(i, j, k, state_fab) != 0 )
		return true;

	if( WindingAxis2(i, j, k, state_fab) != 0 )
		return true;

	if( WindingAxis3(i, j, k, state_fab) != 0 )
		return true;

	return false;
}

inline int AxionStrings::WindingAxis1(int i, int j, int k, amrex::Array4<double const> const& state_fab)
{
	return    ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
			   state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2))
		+ ZeroXing(state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2),
			   state_fab(i+1,j+1,k  ,Scalar::Psi1),state_fab(i+1,j+1,k  ,Scalar::Psi2))
		+ ZeroXing(state_fab(i+1,j+1,k  ,Scalar::Psi1),state_fab(i+1,j+1,k  ,Scalar::Psi2),
			   state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2))
		+ ZeroXing(state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2),
			   state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

inline int AxionStrings::WindingAxis2(int i, int j, int k, amrex::Array4<double const> const& state_fab)
{
	return    ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
			   state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2))
		+ ZeroXing(state_fab(i+1,j  ,k  ,Scalar::Psi1),state_fab(i+1,j  ,k  ,Scalar::Psi2),
			   state_fab(i+1,j  ,k+1,Scalar::Psi1),state_fab(i+1,j  ,k+1,Scalar::Psi2))
		+ ZeroXing(state_fab(i+1,j  ,k+1,Scalar::Psi1),state_fab(i+1,j  ,k+1,Scalar::Psi2),
			   state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2))
		+ ZeroXing(state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2),
			   state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

inline int AxionStrings::WindingAxis3(int i, int j, int k, amrex::Array4<double const> const& state_fab)
{
	return    ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2),
			   state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2))
		+ ZeroXing(state_fab(i  ,j+1,k  ,Scalar::Psi1),state_fab(i  ,j+1,k  ,Scalar::Psi2),
			   state_fab(i  ,j+1,k+1,Scalar::Psi1),state_fab(i  ,j+1,k+1,Scalar::Psi2))
		+ ZeroXing(state_fab(i  ,j+1,k+1,Scalar::Psi1),state_fab(i  ,j+1,k+1,Scalar::Psi2),
			   state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2))
		+ ZeroXing(state_fab(i  ,j  ,k+1,Scalar::Psi1),state_fab(i  ,j  ,k+1,Scalar::Psi2),
			   state_fab(i  ,j  ,k  ,Scalar::Psi1),state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

inline int AxionStrings::ZeroXing(double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2)
{
	if( Psi2_1 * Psi2_2 < 0 ) {
		if( Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0 ){
			return 1;
		}else{
			return -1;
		}
	}
	return 0;
}
