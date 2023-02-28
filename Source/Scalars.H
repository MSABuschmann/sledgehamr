#ifndef Scalars_H_
#define Scalars_H_

#include <boost/preprocessor.hpp>

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

/** @brief Macro to add multiple scalar fields to SledgeHAMR class. 
  *	   All added fields will be simulated.
  *	   Ultimately expands to ADD_SCALARS(a,b,c) -> ScalarField _s_a = {"a", scalar_fields};
  *						       ScalarField _s_b = {"b", scalar_fields};
  *						       ScalarField _s_c = {"c", scalar_fields};
  *						       enum Scalar {a, b, c}; 
  */
#define ADD_SCALARS(...) \
	BOOST_PP_SEQ_FOR_EACH(EXPAND_SCALARS, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
	SCALAR_ENUM(Scalar, __VA_ARGS__)

/** @brief Class to keep track of scalar fields.
 */
class ScalarField{
public:
	/** @brief Upon construction this class will automatically add itself to a 
	 *	   scalar field vector sfv in order to be simulated by the code.
	 *	   It will store its amrex::MultiFab component number as id.
	 * @param	str	Name of scalar field.
	 * @param	sfv	Vector of scalar fields to which this field 
	 *			will be added.
	 */
	ScalarField(std::string str, std::vector<ScalarField*>& sfv)
		: name(str) 
		, id(sfv.size())
	{
		sfv.push_back(this);	
	};

	/** @brief Convenient operator to return internal id
	 */
	operator int() const
	{
		return id;
	};

	/** @brief Name of field by which it will be referred 
	 *	   to in any input and output.
	 */
	std::string name;

	/** @brief Internal id. Corresponds to component in 
	 *	   amrex::Multifab.
	 */
	const int id;	
};

#endif // Scalars_H_
