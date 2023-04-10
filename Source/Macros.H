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
#define START_PROJECT(prj) \
	PRJ_CONSTRUCTOR(prj)

#endif // Macros_H_
