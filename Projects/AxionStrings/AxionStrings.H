#ifndef AxionStrings_H_
#define AxionStrings_H_

#include <SledgeHAMR.H>
#include <SledgeHAMR_Utils.H>

namespace AxionStrings{

ADD_SCALARS(Psi1, Psi2, Pi1, Pi2);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void RHS (const int i, const int j, const int k, const double time, const int
lev, const double dx, amrex::Array4<double> const& rhs_fab, amrex::Array4<double const> const& state_fab)
{
	// Fetch field values.
	double Psi1 = state_fab(i, j, k, Scalar::Psi1);
	double Psi2 = state_fab(i, j, k, Scalar::Psi2);
	double Pi1  = state_fab(i, j, k, Scalar::Pi1);
	double Pi2  = state_fab(i, j, k, Scalar::Pi2);

	double eta = time;

	// Compute Laplacians.
	double dx2 = dx * dx;
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ZeroXing(double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2)
{
	if( Psi2_1 * Psi2_2 >= 0 ) 
		return 0;

	if( Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0 )
		return 1;
	else
		return -1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis1(int i, int j, int k, amrex::Array4<double const> const& state_fab)
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis2(int i, int j, int k, amrex::Array4<double const> const& state_fab)
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis3(int i, int j, int k, amrex::Array4<double const> const& state_fab)
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool TagCellForRefinement(const int i, const int j, const int k, const double time, const int lev,
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

/** @brief Class to simulate axion strings.
 */
class AxionStrings
    : public SledgeHAMR
{
public:
	AxionStrings (){
		scalar_fields = l_scalar_fields;
		amrex::Print() << "Starting AxionStrings project..." << std::endl;
		amrex::Print() << "Number of field components: " << scalar_fields.size()
<< std::endl;
		amrex::Print() << std::endl;
	};

	virtual void FillRHS (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const double time,
					       const amrex::Geometry& geom, int lev) override
	{
		double l_dx = dx[lev];
		DO_PRAGMA(omp parallel if (amrex::Gpu::notInLaunchRegion()))
		for ( amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ){
			const amrex::Box& bx = mfi.tilebox();
			const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi);
			const amrex::Array4<double const>& state_fab = state_mf.array(mfi);
			amrex::ParallelFor(bx,
			[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			{
				RHS(i, j, k, time, lev, l_dx, rhs_fab, state_fab);
			});
		}
	};

	virtual void ErrorEstWithTECPU  (const amrex::Array4<double const>& state_fab,
							const amrex::Array4<double const>& state_fab_te,
							const amrex::Array4<char>& tagarr,
							const amrex::Box& tilebox, double time, int lev,
							int* ntags_total, int* ntags_user, int* ntags_trunc) override
	{
		const amrex::Dim3 lo = amrex::lbound(tilebox);
		const amrex::Dim3 hi = amrex::ubound(tilebox);

		for (int k = lo.z; k <= hi.z; ++k) {
		for (int j = lo.y; j <= hi.y; ++j) {
		AMREX_PRAGMA_SIMD
		for (int i = lo.x; i <= hi.x; ++i) {
			tagarr(i,j,k) = amrex::TagBox::CLEAR;
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);
			if( res ){
				tagarr(i,j,k) = amrex::TagBox::SET;
				(*ntags_user)++;
				(*ntags_total)++;
			}
			bool te_res = SledgeHAMR_Kernels::TruncationErrorTagCPU(i,
j, k, time, lev, state_fab_te, te_crit, ntags_trunc);
			if( te_res ){
				tagarr(i  ,j  ,k  ) = amrex::TagBox::SET;
				tagarr(i+1,j  ,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j+1,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j  ,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j+1,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j+1,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j  ,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j+1,k+1) = amrex::TagBox::SET;
				(*ntags_total) += 8 - (int)res;
			}
		}}}
	};

	virtual void ErrorEstWithTEGPU  (const amrex::Array4<double const>& state_fab,
							const amrex::Array4<double const>& state_fab_te,
							const amrex::Array4<char>& tagarr,
							const amrex::Box& tilebox, double time, int lev) override
	{
		amrex::Gpu::AsyncArray<double> l_te_crit_arr(&te_crit[0],te_crit.size());
		double* l_te_crit = l_te_crit_arr.data();
		amrex::ParallelFor(tilebox,
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			tagarr(i,j,k) = amrex::TagBox::CLEAR;
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);
			if( res ){
				tagarr(i,j,k) = amrex::TagBox::SET;
			}
			bool te_res = SledgeHAMR_Kernels::TruncationErrorTagGPU(i,
j, k, time, lev, state_fab_te, l_te_crit);
			if( te_res ){
				tagarr(i  ,j  ,k  ) = amrex::TagBox::SET;
				tagarr(i+1,j  ,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j+1,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j  ,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j+1,k  ) = amrex::TagBox::SET;
				tagarr(i  ,j+1,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j  ,k+1) = amrex::TagBox::SET;
				tagarr(i+1,j+1,k+1) = amrex::TagBox::SET;
			}
		});
	};

	virtual void ErrorEstWithoutTECPU  (const amrex::Array4<double const>& state_fab,
							const amrex::Array4<char>& tagarr,
							const amrex::Box& tilebox, double time, int lev,
							int* ntags_total) override
	{
		const amrex::Dim3 lo = amrex::lbound(tilebox);
		const amrex::Dim3 hi = amrex::ubound(tilebox);

		for (int k = lo.z; k <= hi.z; ++k) {
		for (int j = lo.y; j <= hi.y; ++j) {
		AMREX_PRAGMA_SIMD
		for (int i = lo.x; i <= hi.x; ++i) {
			tagarr(i,j,k) = amrex::TagBox::CLEAR;
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);
			if( res ){
				tagarr(i,j,k) = amrex::TagBox::SET;
				(*ntags_total)++;
			}
		}}}
	};

	virtual void ErrorEstWithoutTEGPU  (const amrex::Array4<double const>& state_fab,
					    const amrex::Array4<char>& tagarr,
					    const amrex::Box& tilebox, double time, int lev) override
	{
		amrex::ParallelFor(tilebox,
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			tagarr(i,j,k) = amrex::TagBox::CLEAR;
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);
			if( res ){
				tagarr(i,j,k) = amrex::TagBox::SET;
			}
		});
	};

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
//	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//	void RHS (const int i, const int j, const int k, const double time, const int lev,
//	 		    amrex::Array4<double> const& rhs_fab, 
//			    amrex::Array4<double const> const& state_fab);

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
	//AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
	//bool TagCellForRefinement (const int i, const int j, const int k, const double time, const int lev,
	//			   	     amrex::Array4<double const> const& state_fab);
	
private:

	/** @brief Computes the winding factor along a given axis. Will be non-zero
	 *	   if plaquette is pierced by a string.
	 * @param	i		i-th cell index.
	 * @param	j		j-th cell index.
	 * @param	k		k-th cell index.
	 * @param	state_fab	Data.
	 * @return 	Winding factor.
	 */
//	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//	int WindingAxis1 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
//	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//	int WindingAxis2 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
//	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//	int WindingAxis3 (int i, int j, int k, amrex::Array4<double const> const& state_fab);
	
	/** @brief Checks for zero-crossings between two points in the complex scalar field.
	 * @param	Psi1_1	\Psi_1 of 1. point.
	 * @param	Psi2_1	\Psi_2 of 1. point.
	 * @param	Psi1_2	\Psi_1 of 2. point.
	 * @param	Psi2_2	\Psi_2 of 2. point.
	 * @return Sign of slope of zero-crossing. 0 if no crossing.
	 */
//	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//	int ZeroXing (double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2);
};

}; // namespace AxionStrings

#endif //AxionStrings_H_
