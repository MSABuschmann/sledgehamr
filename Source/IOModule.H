#ifndef IOModule_H_
#define IOModule_H_

#include "hdf5.h"

#include <AMReX_ParmParse.H>

#include <SledgeHAMR.H>
#include <OutputModule.H>

class SledgeHAMR;

/** @brief Class that handles all I/O operations besides
 *	   parsing the inputs file.
 */
class IOModule
{
public:
	IOModule (SledgeHAMR * owner);

	/** @brief Writes output if requested.
	 */
	void Write (bool force=false);

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
	 * @param	comp	Component to be filled.
	 * @param	data	Data array.
	 * @param	dimN	Number of cells in each direction
	 *			of level lev.
	 */
	void FillLevelFromArray (int lev, const int comp, double* data, 
				     const unsigned long long dimN);

	/** @brief Fills LevelData with a constant value.
	 * @param	lev	Level to be filled with data.
	 * @param	comp	Component to be filled.
	 * @param	c	Constant.
	 */
	void FillLevelFromConst (int lev, const int comp, const double c);

	/** @brief OUTPUT_FCT. Writes slices along all three directions
	 *	   and all scalar fields.	
	 * @param	time	Current time.
	 * @param	prefix	Output path.
	 */
	void WriteSlices (double time, std::string prefix);

	/** @brief Writes an individual slice along one direction to file
	 *	   for all scalar fields.
	 * @param	time	Current time.
	 * @param	state	Pointer to full grid from which the slice shall be taken.
	 * @param	file_id	ID of HDF5 file.
	 * @param	ident	String identification for slice, e.g. 'x'.
	 * @param	d1	Axis along which the slice shall be taken.
	 * @param	d2	First orthogonal direction.
	 * @param	d3	Second. orthogonal direction.
	 */
	void WriteSingleSlice (double time, const LevelData* state, int lev, hid_t file_id, std::string ident, 
					unsigned int d1, unsigned int d2, unsigned int d3);

	/** @brief Reads dataset from HDF5 file. 
	 * @param	filename	Name of HDF5 file.
	 * @param	dnames		Vector of datasets to be tried.
	 *				First dataset to be found will be read.
	 * @param	data		Data pointer to be filled with data.
	 *				Can be double, float or int.
	 */
	template <typename T> 
	void ReadFromHDF5 (std::string filename, std::vector<std::string> dnames,
				       T* data);

	/** @brief Write dataset to HDF5 file. 
	 * @param	file_id		ID of HDF5 file.
	 * @param	dset		Name of dataset.
	 * @param	data		Pointer to data.
	 *				Can be double, float or int.
	 * @param	size		Size of data.
	 */
	template <typename T> 
	void WriteToHDF5(hid_t file_id, std::string dset, T* data, 
			       unsigned long long size);

	/** @brief Pointer to owner on whose data this class operates.
	 */
	SledgeHAMR * sim;

	/** @brief Vector of output modules
	 */
	std::vector<OutputModule> output;
};

#include "IOModule.tcc"

#endif // IOModule_H_
