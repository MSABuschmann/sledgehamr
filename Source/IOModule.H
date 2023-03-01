#ifndef IOModule_H_
#define IOModule_H_

#include "hdf5.h"

#include <AMReX_ParmParse.H>

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
	/** @brief Fills a given level with chunked data 
	 *	   distributed over many files.
	 * @param	lev	Level to be filled with data.
	 */
	void FillLevelFromFile_Chunks (int lev);

	/** @brief Fills a given level with data from a single hdf5 file.
	 * @param	lev	Level to be filled with data.
	 */
	void FillLevelFromFile_NoChunks (int lev);

	/** @brief Copies data from array into LevelData.
	 * @param	lev	Level to be filled with data.
	 * @param	d	Data array.
	 */
	void FillLevelFromArray (int lev, double* d);

	/** @brief Fills LevelData with a constant value.
	 * @param	lev	Level to be filled with data.
	 * @param	c	Constant.
	 */
	void FillLevelFromConst (int lev, const double c);

	/** @brief Reads dataset from HDF5 file. 
	 * @param	filename	Name of HDF5 file.
	 * @param	dnames		Vector of datasets to be tried.
	 *				First dataset to be found will be read.
	 * @param	data		Data pointer to be filled with data.
	 */
	template <typename T> 
	void ReadFromHDF5(std::string filename, std::vector<std::string> dnames,
				      T* data);

	/** @brief Pointer to owner on whose data this class operates.
	 */
	SledgeHAMR * sim;
};

#include "IOModule.tcc"

#endif // IOModule_H_
