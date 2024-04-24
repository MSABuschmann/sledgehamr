#ifndef SLEDGEHAMR_HDF5_UTILS_H_
#define SLEDGEHAMR_HDF5_UTILS_H_

#include <hdf5.h>

#include <AMReX_AmrCore.H>

namespace sledgehamr {
namespace utils {
namespace hdf5 {

/** @brief Template function that reads dataset of type T from an HDF5 file.
 * @param   filename    HDF5 filename.
 * @param   dnames      List of data set names to try.
 * @param   data        Array pointer to data.
 * @return Whether the read was successfull.
 */
template <typename T>
static bool Read(std::string filename, std::vector<std::string> dnames,
                 T *data) {
    // Identify datatype.
    hid_t mem_type_id;
    if (std::is_same<T, float>::value) {
        mem_type_id = H5T_NATIVE_FLOAT;
    } else if (std::is_same<T, double>::value) {
        mem_type_id = H5T_NATIVE_DOUBLE;
    } else if (std::is_same<T, int>::value) {
        mem_type_id = H5T_NATIVE_INT;
    }

    if (!amrex::FileExists(filename))
        return false;

    // Try and open HDF5 file.
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == H5I_INVALID_HID) {
        amrex::Abort("#error: Could not open file: " + filename);
    }

    // Try and find dataset. Iterate over vector, use first to be found.
    const bool verbose = false;
    std::string dname_conc = " |";
    std::string dname_found = "";
    for (std::string dname : dnames) {
        htri_t exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);
        if (verbose) {
            amrex::Print() << "Attempt reading dataset " << dname << ":"
                           << exists << std::endl;
        }

        if (exists > 0) {
            dname_found = dname;
            break;
        }

        dname_conc += " " + dname;
    }

    if (dname_found == "") {
        H5Fclose(file_id);

        return false;
    }

    // Read dataset.
    hid_t dataset_id = H5Dopen2(file_id, dname_found.c_str(), H5P_DEFAULT);
    if (H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) <
        0) {
        return false;
    }

    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return true;
}

/** @brief Template function to write a dataset of type T to an HDF5 file.
 * @param   file_id HDF5 file id.
 * @param   dset    Dataset name.
 * @param   data    Array pointer to data.
 * @param   size    Length of data.
 */
template <typename T>
static void Write(hid_t file_id, std::string dset, T *data,
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
        amrex::Abort("#error: Writing of dataset " + dset +
                     " failed due to unknown datatype.");
    }

    // Create and write dataset.
    hsize_t dims[1] = {size};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, dset.c_str(), dset_type_id, space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset_id);
}

/** @brief Checks from a list of datasets whether one of them exists in a given
 *         HDF5 file.
 * @param   filename    HDF5 filename.
 * @param   dnames      List of datasets.
 * @return The existing dataset or an empty string of none exist.
 */
static std::string FindDataset(std::string filename,
                               std::vector<std::string> dnames) {
    // Try and open HDF5 file.
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == H5I_INVALID_HID) {
        amrex::Abort("#error: Could not open file: " + filename);
    }

    // Try and find dataset. Iterate over vector, use first to be found.
    std::string dname_found = "";
    for (std::string dname : dnames) {
        htri_t exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);

        if (exists > 0)
            return dname;
    }

    return "";
}

}; // namespace hdf5
}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_HDF5_UTILS_H_
