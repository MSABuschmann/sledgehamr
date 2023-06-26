#ifndef SLEDGEHAMR_IO_MODULE_H_
#define SLEDGEHAMR_IO_MODULE_H_

#include <hdf5.h>

#include <AMReX_ParmParse.H>

#include "sledgehamr.h"
#include "output_module.h"
#include "projection.h"

namespace sledgehamr {

class Sledgehamr;
class Projection;
class Spectrum;

/** @brief Class that handles all I/O operations besides parsing the inputs
 *         file.
 */
class IOModule {
  public:
    IOModule (Sledgehamr* owner);

    /** @brief Writes output if requested.
     */
    void Write(bool force=false);

    /** @brief Reads dataset from HDF5 file.
     * @param   filename    Name of HDF5 file.
     * @param   dnames      Vector of datasets to be tried. First dataset to be
     *                      found will be read.
     * @param   data        Data pointer to be filled with data. Can be double,
     *                      float or int. TODO: Add static_assert.
     */
    template <typename T>
    static void ReadFromHDF5(std::string filename,
                             std::vector<std::string> dnames, T* data);

    /** @brief Write dataset to HDF5 file.
     * @param   file_id ID of HDF5 file.
     * @param   dset    Name of dataset.
     * @param   data    Pointer to data. Can be double, float or int. TODO: Add
     *                  static_assert
     * @param   size    Size of data.
     */
    template <typename T>
    static void WriteToHDF5(hid_t file_id, std::string dset, T* data,
                     unsigned long long size);

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

    /** @brief Vectors containing instructions for projections and spectra.
     */
    std::vector<Projection> projections;
    std::vector<Spectrum> spectra;

  private:
    /** @brief Copies data from array into LevelData.
     * @param   lev Level to be filled with data.
     * @param   comp    Component to be filled.
     * @param   data    Data array.
     * @param   dimN    Number of cells in each direction of level lev.
     */
    void FillLevelFromArray(int lev, const int comp, double* data,
                            const long long dimN);

    /** @brief Fills a given level with chunked data distributed over many
     *         files.
     * @param   lev Level to be filled with data.
     */
    void FillLevelFromFile_Chunks(int lev);

    /** @brief Fills a given level with data from a single hdf5 file.
     * @param   lev Level to be filled with data.
     */
    void FillLevelFromFile_NoChunks(int lev);

    /** @brief OUTPUT_FCT. Wrapper to write slices along all three directions
     *         and all scalar fields.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteSlices(double time, std::string prefix);

    /** @brief OUTPUT_FCT. Wrapper to write slices of truncation errors and
     *         corresponding fields.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteSlicesTruncationError(double time, std::string prefix);

    /** @brief Writes slices along all three directions and all scalar fields
               with or without truncation errors.
     * @param   time                    Current time.
     * @param   prefix                  Output path.
     * @param   with_truncation_errors  Whether truncation errors shall be
     *                                  written.
     */
    void DoWriteSlices(double time, std::string prefix,
                       bool with_truncation_errors);

    /** @brief Writes an individual slice along one direction to file for all
     *         scalar fields.
     * @param   time                    Current time.
     * @param   state                   Pointer to full grid from which the
     *                                  slice shall be taken.
     * @param   file_id                 ID of HDF5 file.
     * @param   ident                   String identifier, e.g. 'x'.
     * @param   d1                      Axis along which the slice shall be
     *                                  taken.
     * @param   d2                      First orthogonal direction.
     * @param   d3                      Second. orthogonal direction.
     * @param   is_truncation_errors    Whether the data contains truncation
     *                                  errors.
     */
    void WriteSingleSlice(double time, const LevelData* state, int lev,
                          hid_t file_id, std::string ident, int d1, int d2,
                          int d3, bool is_truncation_errors);

    /** @brief OUTPUT_FCT. Wrapper to write the coarse level.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteCoarseBox(double time, std::string prefix);

    /** @brief OUTPUT_FCT. Wrapper to write the coarse level with truncation
     *         errors.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteCoarseBoxTruncationError(double time, std::string prefix);

    /** @brief Writes the coarse level possibly with truncation errors.
     * @param   time                    Current time.
     * @param   prefix                  Output path.
     * @param   downsample_factor       Downsampling factor.
     * @param   with_truncation_errors  Whether truncation errors shall be
     *                                  written.
     */
    void DoWriteCoarseBox(double time, std::string prefix,
                          int downsample_factor, bool with_truncation_errors);

    /** @brief Writes an entire level to file.
     * @param   time                    Current time.
     * @param   state                   Pointer to full grid from which the
     *                                  slice shall be taken.
     * @param   ident                   String identifier, e.g. "data" or "te".
     * @param   file_id                 ID of HDF5 file.
     * @param   downsample_factor       Downsample level by this factor.
     * @param   is_truncation_errors    Whether the data contains truncation
     *                                  errors.
     */
    void WriteLevel(double time, const LevelData* state, int lev,
                    hid_t file_id, std::string ident, int downsample_factor,
                    bool is_truncation_errors);

    /** @brief OUTPUT_FCT. Wrapper to write all levels.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteFullBox(double time, std::string prefix);

    /** @brief OUTPUT_FCT. Wrapper to write truncation errors on all levels.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteFullBoxTruncationError(double time, std::string prefix);

    /** @brief Writes all levels possibly with truncation errors.
     * @param   time                    Current time.
     * @param   prefix                  Output path.
     * @param   downsample_factor       Downsampling factor.
     * @param   with_truncation_errors  Whether truncation errors shall be
     *                                  written.
     */
    void DoWriteFullBox(double time, std::string prefix, int downsample_factor,
                        bool with_truncation_errors);

    /** @brief OUTPUT_FCT. Write projections.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteProjections(double time, std::string prefix);

    /** @brief OUTPUT_FCT. Write spectra.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteSpectra(double time, std::string prefix);

    /** @brief OUTPUT_FCT. Write gravitational wave spectrum.
     * @param   time   Current time.
     * @param   prefix Output path.
     */
    void WriteGravitationalWaveSpectrum(double time, std::string prefix);


    /** Downsampling factors for coarse/full level output.
     */
    int coarse_box_downsample_factor = 1;
    int coarse_box_truncation_error_downsample_factor = 1;
    int full_box_downsample_factor = 1;
    int full_box_truncation_error_downsample_factor = 1;

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
        amrex::Abort("#error: Could not open file: " + filename);
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
                "#error: Could not find correct dataset in file: "
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
                           unsigned long long size) {
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
