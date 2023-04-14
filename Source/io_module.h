#ifndef SLEDGEHAMR_IO_MODULE_H_
#define SLEDGEHAMR_IO_MODULE_H_

#include <hdf5.h>

#include <AMReX_ParmParse.H>

#include "sledgehamr.h"
#include "output_module.h"

namespace sledgehamr {

class Sledgehamr;

/** @brief Class that handles all I/O operations besides parsing the inputs
 *         file.
 */
class IOModule {
  public:
    IOModule (Sledgehamr* owner);

    /** @brief Writes output if requested.
     */
    void Write(bool force=false);

    /** @brief Fills a given level with data from hdf5 file(s).
     * @param   lev Level to be filled with data.
     */
    void FillLevelFromFile(int lev);

    /** @brief Fills LevelData with a constant value.
     * @param   lev     Level to be filled with data.
     * @param   comp    Component to be filled.
     * @param   c       Constant.
     */
    void FillLevelFromConst(int lev, const int comp, const double c);

  private:

    /** @brief Fills a given level with chunked data distributed over many
     *         files.
     * @param   lev Level to be filled with data.
     */
    void FillLevelFromFile_Chunks(int lev);

    /** @brief Fills a given level with data from a single hdf5 file.
     * @param   lev Level to be filled with data.
     */
    void FillLevelFromFile_NoChunks(int lev);

    /** @brief Copies data from array into LevelData.
     * @param   lev Level to be filled with data.
     * @param   comp    Component to be filled.
     * @param   data    Data array.
     * @param   dimN    Number of cells in each direction of level lev.
     */
    void FillLevelFromArray(int lev, const int comp, double* data,
                            const long long dimN);

    /** @brief OUTPUT_FCT. Writes slices along all three directions and all
     *         scalar fields.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteSlices(double time, std::string prefix);

    /** @brief Writes an individual slice along one direction to file for all
     *         scalar fields.
     * @param   time    Current time.
     * @param   state   Pointer to full grid from which the slice shall be
     *                  taken.
     * @param   file_id ID of HDF5 file.
     * @param   ident   String identification for slice, e.g. 'x'.
     * @param   d1      Axis along which the slice shall be taken.
     * @param   d2      First orthogonal direction.
     * @param   d3      Second. orthogonal direction.
     */
    void WriteSingleSlice(double time, const LevelData* state, int lev,
                          hid_t file_id, std::string ident, int d1, int d2,
                          int d3);

    /** @brief Reads dataset from HDF5 file.
     * @param   filename    Name of HDF5 file.
     * @param   dnames      Vector of datasets to be tried. First dataset to be
     *                      found will be read.
     * @param   data        Data pointer to be filled with data. Can be double,
     *                      float or int. TODO: Add static_assert.
     */
    template <typename T>
    void ReadFromHDF5(std::string filename, std::vector<std::string> dnames,
                      T* data);

    /** @brief Write dataset to HDF5 file.
     * @param   file_id ID of HDF5 file.
     * @param   dset    Name of dataset.
     * @param   data    Pointer to data. Can be double, float or int. TODO: Add
     *                  static_assert
     * @param   size    Size of data.
     */
    template <typename T>
    void WriteToHDF5(hid_t file_id, std::string dset, T* data, long long size);

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;

    /** @brief Vector of output modules
     */
    std::vector<OutputModule> output;
};

template <typename T>
void IOModule::ReadFromHDF5(std::string filename,
                            std::vector<std::string> dnames, T* data) {
    // Identify datatype.
    hid_t mem_type_id;
    if (std::is_same<T, float>::value) {
        mem_type_id = H5T_NATIVE_FLOAT;
    } else if (std::is_same<T, double>::value) {
        mem_type_id = H5T_NATIVE_DOUBLE;
    } else if (std::is_same<T, int>::value) {
        mem_type_id = H5T_NATIVE_INT;
    }

    // Try and open HDF5 file.
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == H5I_INVALID_HID) {
        amrex::Abort("#error: Could not open initial state file: " + filename);
    }

    // Try and find dataset. Iterate over vector, use first to be found.
    std::string dname_conc = " |";
    std::string dname_found = "";
    for (std::string dname : dnames) {
        htri_t exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);

        if (exists > 0) {
            dname_found = dname;
            break;
        }

        dname_conc += " " + dname;
    }

    if (dname_found  == "") {
        H5Fclose(file_id);
        amrex::Abort(
                "#error: Could not find correct dataset in initial state file: "
                + filename + dname_conc);
    }

    // Read dataset.
    hid_t dataset_id = H5Dopen2(file_id, dname_found.c_str(), H5P_DEFAULT);
    H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
}

template <typename T>
void IOModule::WriteToHDF5(hid_t file_id, std::string dset, T* data,
                           long long size) {
    // Identify datatype.
    hid_t mem_type_id, dset_type_id;
    if (std::is_same<T, float>::value) {
        mem_type_id = H5T_NATIVE_FLOAT;
        dset_type_id = H5T_IEEE_F32LE;
    } else if (std::is_same<T, double>::value) {
        mem_type_id = H5T_NATIVE_DOUBLE;
        dset_type_id = H5T_IEEE_F64LE;
    } else if (std::is_same<T, int>::value) {
        mem_type_id = H5T_NATIVE_INT;
        dset_type_id = H5T_IEEE_F64LE;
    } else {
        amrex::Abort("#error: Writing of dataset " + dset
                     + " failed due to unknown datatype.");
    }

    // Create and write dataset.
    hsize_t dims[1] = {size};
    hid_t space = H5Screate_simple(1,dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, dset.c_str(), dset_type_id, space,
                                 H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset_id);
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_IO_MODULE_H_
