#ifndef SledgeHAMR_Utils_H_
#define SledgeHAMR_Utils_H_

namespace SledgeHAMR_Utils{

	/** @brief Calculates the Laplacian at a given order. Only 0th, 1st and 2nd order 
	 *	   aren currently implemented.
	 * @param	state_fab	Data from which Laplacian is to be calculated.
	 * @param	i		i-th center cell.
	 * @param	j		j-th center cell.
	 * @param	k		k-th center cell.
	 * @param	c		Scalar component.
	 * @param	dx2		Squared grid spacing.
	 */
	template<unsigned int>
	inline double Laplacian(amrex::Array4<amrex::Real const> const& state_fab,
					const int i, const int j, const int k, 
					const int c, const double dx2);

	template <>
	inline double Laplacian<0>(amrex::Array4<amrex::Real const> const& state_fab,
					const int i, const int j, const int k, 
					const int c, const double dx2)
	{
		return 0;
	};

	template <>
	inline double Laplacian<1>(amrex::Array4<amrex::Real const> const& state_fab,
					const int i, const int j, const int k, 
					const int c, const double dx2)
	{
		return (state_fab(i + 1, j, k, c) + state_fab(i - 1, j, k, c) +
			state_fab(i, j + 1, k, c) + state_fab(i, j - 1, k, c) +
			state_fab(i, j, k + 1, c) + state_fab(i, j, k - 1, c) -
			6.*state_fab(i,j,k,c)) / dx2;
	};

	template <>
	inline double Laplacian<2>(amrex::Array4<amrex::Real const> const& state_fab,
					const int i, const int j, const int k, 
					const int c, const double dx2)
	{
		return ( - (state_fab(i + 2, j, k, c) + state_fab(i - 2, j, k, c) +
			    state_fab(i, j + 2, k, c) + state_fab(i, j - 2, k, c) +
			    state_fab(i, j, k + 2, c) + state_fab(i, j, k - 2, c) )
		      + 16.*(state_fab(i + 1, j, k, c) + state_fab(i - 1, j, k, c) +
			     state_fab(i, j + 1, k, c) + state_fab(i, j - 1, k, c) +
			     state_fab(i, j, k + 1, c) + state_fab(i, j, k - 1, c) )
		      - 90.*state_fab(i,j,k,c) ) / (12.*dx2);
	};

	/** @brief Returns string of ordinal number suffix for small positive
	 *	   numbers.
	 * @param	num	Number.
	 */
	static std::string OrdinalNumberSuffix(int num)
	{
		if( num == 1 ) return "st";
		if( num == 2 ) return "nd";
		if( num == 3 ) return "rd";
		return "th";
	}

	/** @brief Returns string containing the name of the level.
	 * @param	lev		Level.
	 * @param	shadow_hierachy Whether a shadow hierarchy 
	 *				has been used.
	 */
	static std::string LevelName(int lev, bool shadow_hierarchy)
	{
		if( lev==0 && shadow_hierarchy )
			return "shadow level";

		if( (lev==0 && !shadow_hierarchy) ||
		    (lev==1 &&  shadow_hierarchy) )
			return "coarse level";

		int refnum = lev - shadow_hierarchy;
		return std::to_string(refnum) + OrdinalNumberSuffix(refnum)
			+ " refinement";
	};
};

#endif // SledgeHAMR_Utils_H_
