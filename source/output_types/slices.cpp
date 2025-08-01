#include "slices.h"
#include "hdf5_utils.h"

namespace sledgehamr {

void Slices::ParseParams() {
    amrex::ParmParse pp_prj("output.slices");
    pp_prj.queryarr("location", slice_location, 0, 3);
}

/** @brief Writes slices along all the directions through all scalar fields and
 *         levels.
 */
void Slices::Write() {
    for (int lev = 0; lev <= sim->GetFinestLevel(); ++lev) {
        const LevelData *state = &sim->GetLevelData(lev);
        if (with_truncation_errors && !state->contains_truncation_errors)
            continue;

        // Create folder and file.
        std::string subfolder = folder + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(subfolder.c_str(), 0755);

        std::string filename =
            subfolder + "/" +
            std::to_string(amrex::ParallelDescriptor::MyProc()) + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        // Write field data.
        WriteSingleSlice(state, lev, file_id, "x", 0, 1, 2, false);
        WriteSingleSlice(state, lev, file_id, "y", 1, 0, 2, false);
        WriteSingleSlice(state, lev, file_id, "z", 2, 0, 1, false);

        if (with_truncation_errors) {
            // Write truncation errors.
            const LevelData *state_old = &sim->GetOldLevelData(lev);
            WriteSingleSlice(state_old, lev, file_id, "te_x", 0, 1, 2, true);
            WriteSingleSlice(state_old, lev, file_id, "te_y", 1, 0, 2, true);
            WriteSingleSlice(state_old, lev, file_id, "te_z", 2, 0, 1, true);
        }

        H5Fclose(file_id);
    }
}

/** @brief Writes a single slices.
 * @param   state   State data.
 * @param   lev     Current level.
 * @param   file_id HDF5 file id.
 * @param   ident   Unique identifier string.
 * @param   d1      Orientation 1.
 * @param   d2      Orientation 2.
 * @param   d3      Orientation 3.
 * @param   is_truncation_error Whether state data contains truncation error
 *                              estimates.
 */
void Slices::WriteSingleSlice(const LevelData *state, int lev, hid_t file_id,
                              std::string ident, int d1, int d2, int d3,
                              bool is_truncation_error) {
    std::vector<int> le1, he1, le2, he2;
    const int ndist = is_truncation_error ? 2 : 1;

    const double nDim = static_cast<double>(sim->GetDimN(lev));
    int slice_ind = static_cast<int>(nDim * slice_location[d1]);
    // In case we are writing truncation errors we need to make sure we pick a
    // slice that has them.
    if (slice_ind % 2 == 1) {
        slice_ind++;
    }

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox();

        // Find box that cuts through selected plane.
        if (bx.smallEnd(d1) <= slice_ind && bx.bigEnd(d1) >= slice_ind) {
            const auto &state_arr = state->array(mfi);

            // Get box dimensions
            int l1 = bx.smallEnd(d2);
            int l2 = bx.smallEnd(d3);
            int h1 = bx.bigEnd(d2) + 1;
            int h2 = bx.bigEnd(d3) + 1;
            le1.push_back(l1);
            le2.push_back(l2);
            he1.push_back(h1);
            he2.push_back(h2);

            int dim1 = (h1 - l1) / ndist;
            int dim2 = (h2 - l2) / ndist;

            // Copy data into flattened array for each scalar field.
            long len = dim1 * dim2;
            // TODO Adjust output type.
            //            std::unique_ptr<float[]> output_arr(new float[len]);
            //            std::fill_n(output_arr.get(), len, 0.0f);
            std::vector<float> output_arr(len);

            for (int f = 0; f < state->nComp(); ++f) {
                for (int j = l2; j < h2; ++j) {
                    for (int i = l1; i < h1; ++i) {
                        if (is_truncation_error && (i % 2 != 0 || j % 2 != 0))
                            continue;

                        int ind = (i - l1) / ndist * dim2 + (j - l2) / ndist;

                        if (d1 == 0) {
                            output_arr[ind] = state_arr(slice_ind, i, j, f);
                        } else if (d1 == 1) {
                            output_arr[ind] = state_arr(i, slice_ind, j, f);
                        } else if (d1 == 2) {
                            output_arr[ind] = state_arr(i, j, slice_ind, f);
                        }
                    }
                }

                std::string dset_name = sim->GetScalarFieldName(f) + "_" +
                                        ident + "_" +
                                        std::to_string(le1.size());
                utils::hdf5::Write(file_id, dset_name, &(output_arr[0]), len);
                // utils::hdf5::Write(file_id, dset_name, output_arr.get(),
                // len);
            }
            //            output_arr.reset();
        }
    }

    // Write header information for this slice.
    const int nparams = 5;
    double header_data[nparams] = {
        state->t, (double)amrex::ParallelDescriptor::NProcs(),
        (double)(sim->GetFinestLevel()), (double)sim->GetDimN(lev),
        (double)le1.size()};
    utils::hdf5::Write(file_id, "Header_" + ident, header_data, nparams);

    // Write box dimensions so we can reassemble slice.
    if (le1.size() == 0)
        return;

    utils::hdf5::Write(file_id, "le1_" + ident, (int *)&(le1[0]), le1.size());
    utils::hdf5::Write(file_id, "le2_" + ident, (int *)&(le2[0]), le2.size());
    utils::hdf5::Write(file_id, "he1_" + ident, (int *)&(he1[0]), he1.size());
    utils::hdf5::Write(file_id, "he2_" + ident, (int *)&(he2[0]), he2.size());
}
}; // namespace sledgehamr
