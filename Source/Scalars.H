#ifndef Scalars_H_
#define Scalars_H_

#include <boost/preprocessor.hpp>

#define ADD_SCALAR(field) ScalarField BOOST_PP_CAT(_s_, field)  = {#field, scalar_fields};
#define EXPAND_SCALARS(r, data, field) ADD_SCALAR(field)

#define SCALAR_ENUM_VALUE(r, data, elem) elem,
#define SCALAR_ENUM(name, ...) \
    enum name { \
        BOOST_PP_SEQ_FOR_EACH(SCALAR_ENUM_VALUE, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    };

#define ADD_SCALARS(...) \
	BOOST_PP_SEQ_FOR_EACH(EXPAND_SCALARS, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
	SCALAR_ENUM(Scalar, __VA_ARGS__)

class ScalarField{
public:
	ScalarField(std::string str, std::vector<ScalarField*>& sfv)
		: name(str) 
		, id(sfv.size())
	{
		sfv.push_back(this);	
	};

	operator int() const
	{
		return id;
	};

	std::string name;
	const int id;	
};

#endif // Scalars_H_
