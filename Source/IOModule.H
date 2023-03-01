#ifndef IOModule_H_
#define IOModule_H_

#include <SledgeHAMR.H>

class SledgeHAMR;

/** @brief Class that handles all I/O operations besides
 *	   parsing the inputs file.
 */
class IOModule
{
public:
	IOModule (SledgeHAMR * owner);

	/** @brief Fills a given level with data from hdf5 file(s).
	 * @param	lev	Level to be filled with data.
	 */
	void FillLevelFromFile (int lev);

private:
	/** @brief Pointer to owner on whose data this class operates.
	 */
	SledgeHAMR * sim;
};

#endif // IOModule_H_
