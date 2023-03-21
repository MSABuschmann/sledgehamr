#ifndef Macros_H_
#define Macros_H_

#define DO_PRAGMA(x) _Pragma(#x)

/** @brief Macro to override the SledgeHAMR::FillRHS function within the project class.
 *	   Will autmatically be added by the ADD_SCALARS macro to the project. This is 
 *	   needed to ensure the RHS function can be inlined by the compiler.
 */
#define FILL_RHS virtual void FillRHS (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const double time,\
					       const amrex::Geometry& geom, int lev) override\
	{\
		DO_PRAGMA(omp parallel)\
		for ( amrex::MFIter mfi(rhs_mf, true); mfi.isValid(); ++mfi ){\
			const amrex::Box& bx = mfi.tilebox();\
			const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi);\
			const amrex::Array4<double const>& state_fab = state_mf.array(mfi);\
			amrex::ParallelFor(bx,\
			[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept\
			{\
				RHS(i, j, k, time, lev, rhs_fab, state_fab);\
			});\
		}\
	};

/** @brief Macro to override the SledgeHAMR::ErrorEstWithTE function within the project class.
 *	   Will autmatically be added by the ADD_SCALARS macro to the project. This is needed
 *	   to ensure that e.g. the TagCellForRefinement function can be inlined by the compiler.
 */
#define ERROR_EST_WITH_TE virtual void ErrorEstWithTE  (const amrex::Array4<double const>& state_fab,\
							const amrex::Array4<double const>& state_fab_te,\
							const amrex::Array4<char>& tagarr,\
							const amrex::Box& tilebox, double time, int lev,\
							int* ntags_user)\
	{\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_user)++;\
			}\
		});\
	};\

/** @brief Same as ERROR_EST_WITH_TE macro but for SledgeHAMR::ErrorEstWithoutTE function.
 */
#define ERROR_EST_WITHOUT_TE virtual void ErrorEstWithoutTE  (const amrex::Array4<double const>& state_fab,\
							      const amrex::Array4<double const>& state_fab_te,\
							      const amrex::Array4<char>& tagarr,\
							      const amrex::Box& tilebox, double time, int lev,\
							      int* ntags_user)\
	{\
		amrex::ParallelFor(tilebox,\
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {\
			tagarr(i,j,k) = amrex::TagBox::CLEAR;\
			bool res = TagCellForRefinement(i, j, k, time, lev, state_fab);\
			if( res ){\
				tagarr(i,j,k) = amrex::TagBox::SET;\
				(*ntags_user)++;\
			}\
		});\
	};\

/** @brief Macros to declare and initialize scalar fields within SledgeHAMR class.
 */
#define ADD_SCALAR(field) ScalarField BOOST_PP_CAT(_s_, field)  = {#field, scalar_fields};
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
  *	   Ultimately expands to ADD_SCALARS(a,b,c) -> ScalarField _s_a = {"a", scalar_fields};
  *						       ScalarField _s_b = {"b", scalar_fields};
  *						       ScalarField _s_c = {"c", scalar_fields};
  *						       enum Scalar {a, b, c}; 
  *						       virtual void FillRHS ( ...
  *						       virtual void ErrorEstWithTE ( ...
  */
#define ADD_SCALARS(...) \
	BOOST_PP_SEQ_FOR_EACH(EXPAND_SCALARS, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
	SCALAR_ENUM(Scalar, __VA_ARGS__)\
	FILL_RHS\
	ERROR_EST_WITH_TE\
	ERROR_EST_WITHOUT_TE

#endif // Macros_H_
