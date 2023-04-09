#ifndef Scalars_H_
#define Scalars_H_

#include <SledgeHAMR.H>
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
