#ifndef Macros_H_
#define Macros_H_

// Force boost to define variadics. 
#define BOOST_PP_VARIADICS 1
#include <boost/preprocessor.hpp>

/** @brief Element-wise vector reduction for OpenMP
 */
#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = omp_orig)

#define DO_PRAGMA(x) _Pragma(#x)

/** @brief Macros to declare and initialize scalar fields within SledgeHAMR class.
 */
#define ADD_SCALAR(field) static ScalarField BOOST_PP_CAT(_s_, field)  = {#field, l_scalar_fields};
#define EXPAND_SCALARS(r, data, field) ADD_SCALAR(field)

/** @brief Macros to create enum of fields for fast and convinient component access
 *	   within the SledgeHAMR class.
 */
#define SCALAR_ENUM_VALUE(r, data, elem) elem,
#define SCALAR_ENUM(name, ...) \
    enum name { \
        BOOST_PP_SEQ_FOR_EACH(SCALAR_ENUM_VALUE, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    };

/** @brief Macro to add multiple scalar fields to project class. All added fields will be
  *	   simulated. Also overrides SledgeHAMR::FillRHS and SledgeHAMR::ErrorEstWithTE.
  *	   Ultimately expands to ADD_SCALARS(a,b,c) -> static ScalarField _s_a = {"a", scalar_fields};
  *						       static ScalarField _s_b = {"b", scalar_fields};
  *						       static ScalarField _s_c = {"c", scalar_fields};
  *						       enum Scalar {a, b, c}; 
  */
#define ADD_SCALARS(...) \
	static std::vector<ScalarField*> l_scalar_fields;\
	BOOST_PP_SEQ_FOR_EACH(EXPAND_SCALARS, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
	SCALAR_ENUM(Scalar, __VA_ARGS__)\

/** @brief TODO
 */
#define PRJ_CONSTRUCTOR(prj) \
	prj (){ \
		scalar_fields = l_scalar_fields;\
		amrex::Print() << "Starting "  << #prj << " project..." << std::endl;\
		amrex::Print() << "Number of field components: " << scalar_fields.size() << std::endl;\
		amrex::Print() << std::endl;\
	};

/** @brief TODO
 */
#define PRJ_FILLRHS virtual void FillRhs \
	(amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const double time,\
					       const amrex::Geometry& geom, int lev) override\
	{\
		double l_dx = dx[lev];\
		DO_PRAGMA(omp parallel if (amrex::Gpu::notInLaunchRegion()))\
		for ( amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ){\
			const amrex::Box& bx = mfi.tilebox();\
			const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi);\
			const amrex::Array4<double const>& state_fab = state_mf.array(mfi);\
			amrex::ParallelFor(bx,\
			[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept\
			{\
				Rhs(i, j, k, time, lev, l_dx, rhs_fab, state_fab);\
			});\
		}\
	};
/** @brief TODO
 */
#define PRJ_ERRORESTWITHTECPU virtual void ErrorEstWithTECPU \
 	(const amrex::Array4<double const>& state_fab,\
	const amrex::Array4<double const>& state_fab_te,\
	const amrex::Array4<char>& tagarr,\
	const amrex::Box& tilebox, double time, int lev,\
	int* ntags_total, int* ntags_user, int* ntags_trunc) override\
	{\
		const amrex::Dim3 lo = amrex::lbound(tilebox);\
		const amrex::Dim3 hi = amrex::ubound(tilebox);\
		for (int k = lo.z; k <= hi.z; ++k) {\
		for (int j = lo.y; j <= hi.y; ++j) {\
		AMREX_PRAGMA_SIMD\
		for (int i = lo.x; i <= hi.x; ++i) {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = false;\
			res = TagCellForRefinement( i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_user)++;\
				(*ntags_total)++;\
			}\
			bool te_res =\
SledgeHAMR_Kernels::TruncationErrorTagCPU(i,\
j, k, time, lev, state_fab_te, te_crit, ntags_trunc);\
			if( te_res ){\
				tagarr(i  ,j  ,k  ) = amrex::TagBox::SET;\
				tagarr(i+1,j  ,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j+1,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j  ,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j+1,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j+1,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j  ,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j+1,k+1) = amrex::TagBox::SET;\
				(*ntags_total) += 8 - (int)res;\
			}\
		}}}\
	};

/** @brief TODO
 */
#define PRJ_ERRORESTWITHTEGPU virtual void ErrorEstWithTEGPU  (const amrex::Array4<double const>& state_fab,\
							const amrex::Array4<double const>& state_fab_te,\
							const amrex::Array4<char>& tagarr,\
							const amrex::Box& tilebox, double time, int lev) override\
	{\
		amrex::Gpu::AsyncArray<double> l_te_crit_arr(&te_crit[0],te_crit.size());\
		double* l_te_crit = l_te_crit_arr.data();\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
			}\
			bool te_res = SledgeHAMR_Kernels::TruncationErrorTagGPU(i,\
j, k, time, lev, state_fab_te, l_te_crit);\
			if( te_res ){\
				tagarr(i  ,j  ,k  ) = amrex::TagBox::SET;\
				tagarr(i+1,j  ,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j+1,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j  ,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j+1,k  ) = amrex::TagBox::SET;\
				tagarr(i  ,j+1,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j  ,k+1) = amrex::TagBox::SET;\
				tagarr(i+1,j+1,k+1) = amrex::TagBox::SET;\
			}\
		});\
	};

/** @brief TODO
 */
#define PRJ_ERRORESTWITHOUTTECPU virtual void ErrorEstWithoutTECPU  (const amrex::Array4<double const>& state_fab,\
							const amrex::Array4<char>& tagarr,\
							const amrex::Box& tilebox, double time, int lev,\
							int* ntags_total) override\
	{\
		const amrex::Dim3 lo = amrex::lbound(tilebox);\
		const amrex::Dim3 hi = amrex::ubound(tilebox);\
		for (int k = lo.z; k <= hi.z; ++k) {\
		for (int j = lo.y; j <= hi.y; ++j) {\
		AMREX_PRAGMA_SIMD\
		for (int i = lo.x; i <= hi.x; ++i) {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_total)++;\
			}\
		}}}\
	};

/** @brief TODO
 */
#define PRJ_ERRORESTWITHOUTTEGPU virtual void ErrorEstWithoutTEGPU  (const amrex::Array4<double const>& state_fab,\
					    const amrex::Array4<char>& tagarr,\
					    const amrex::Box& tilebox, double time, int lev) override\
	{\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
			}\
		});\
	};

/** @brief TODO
 */
#define START_PROJECT(prj) \
	PRJ_CONSTRUCTOR(prj)\
	PRJ_FILLRHS\
	PRJ_ERRORESTWITHTECPU\
	PRJ_ERRORESTWITHTEGPU\
	PRJ_ERRORESTWITHOUTTECPU\
	PRJ_ERRORESTWITHOUTTEGPU

#endif // Macros_H_
