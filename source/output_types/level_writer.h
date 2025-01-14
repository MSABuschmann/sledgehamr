#ifndef SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
#define SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_

#include "hdf5_utils.h"
#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Writes level data to disk.
 */
class LevelWriter {
  public:
    LevelWriter(Sledgehamr *owner, std::string prefix, int output_type);
    void Write();

  private:
    void DetermineSetup();
    void ParseParams();
    void CheckDownsampleFactor();

    template <typename T>
    void WriteSingleLevel(const LevelData *state, int lev, hid_t file_id,
                          std::string ident, bool is_truncation_error);

    /** @brief Pointer to simulation.
     */
    Sledgehamr *sim;

    /** @brief Local output folder.
     */
    std::string folder;

    /** @brief Mininum level to write (inclusive).
     */
    int level_min;

    /** @brief Maximum level to write (inclusive).
     */
    int level_max;

    /** @brief Local unique output id. Must be either
     *          - io_module->idx_coarse_box
     *          - io_module->idx_coarse_box_truncation_error
     *          - io_module->idx_full_box
     *          - io_module->idx_full_box_truncation_error
     */
    const int output_id;

    /** @brief Display name of output type.
     */
    std::string name;

    /** @brief Information to print to screen for this output type.
     */
    std::string info;

    /** @brief Whether we are including truncation errors estimates.
     */
    bool with_truncation_errors;

    /** @brief Downsample data by this factor before writing.
     */
    int downsample_factor = 1;

    /** @brief floating point precision.
     */
    int precision = 32;
};

/** @brief Writes a single level to disk.
 * @param   state               State.
 * @param   lev                 Level.
 * @param   file_id             HDF5 file to write to.
 * @param   ident               String identifier for dataset.
 * @param   is_truncation_error Whether state contains truncation errors.
 */
template <typename T>
void LevelWriter::WriteSingleLevel(const LevelData *state, int lev,
                                   hid_t file_id, std::string ident,
                                   bool is_truncation_error) {
    const int ndist = is_truncation_error ? 2 : 1;
    const int grid_density = downsample_factor * ndist;

    std::vector<int> lex, hex, ley, hey, lez, hez;

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox();
        const auto &state_arr = state->array(mfi);

        // Get box dimensions
        int lx = bx.smallEnd(0);
        int ly = bx.smallEnd(1);
        int lz = bx.smallEnd(2);
        int hx = bx.bigEnd(0) + 1;
        int hy = bx.bigEnd(1) + 1;
        int hz = bx.bigEnd(2) + 1;
        lex.push_back(lx);
        ley.push_back(ly);
        lez.push_back(lz);
        hex.push_back(hx);
        hey.push_back(hy);
        hez.push_back(hz);

        int dimx = (hx - lx) / grid_density;
        int dimy = (hy - ly) / grid_density;
        int dimz = (hz - lz) / grid_density;
        double volfac = 1. / std::pow(downsample_factor, 3);

        // Copy data into flattened array for each scalar field.
        long len = dimx * dimy * dimz;

        // TODO Adjust output type.
        for (int f = 0; f < state->nComp(); ++f) {
            // std::unique_ptr<float[]> output_arr(new float[len]);
            // std::fill_n(output_arr.get(), len, 0.0f);
            std::vector<T> output_arr(len);

            for (int k = lz; k < hz; ++k) {
                for (int j = ly; j < hy; ++j) {
                    for (int i = lx; i < hx; ++i) {
                        if (is_truncation_error &&
                            (i % 2 != 0 || j % 2 != 0 || k % 2 != 0))
                            continue;

                        int ind = (i - lx) / grid_density * dimy * dimz +
                                  (j - ly) / grid_density * dimz +
                                  (k - lz) / grid_density;

                        if (is_truncation_error) {
                            output_arr[ind] =
                                std::max(static_cast<T>(state_arr(i, j, k, f)),
                                         output_arr[ind]);
                        } else {
                            output_arr[ind] += state_arr(i, j, k, f) * volfac;
                        }
                    }
                }
            }

            std::string dset_name = sim->GetScalarFieldName(f) + "_" + ident +
                                    "_" + std::to_string(lex.size());
            //    utils::hdf5::Write(file_id, dset_name, output_arr.get(), len);
            utils::hdf5::Write(file_id, dset_name, &(output_arr[0]), len);
            //    output_arr.reset();
        }
    }

    // Write header information for this slice.
    const int nparams = 6;
    double header_data[nparams] = {state->t,
                                   (double)amrex::ParallelDescriptor::NProcs(),
                                   (double)sim->GetFinestLevel(),
                                   (double)sim->GetDimN(lev),
                                   (double)downsample_factor,
                                   (double)lex.size()};
    utils::hdf5::Write(file_id, "Header_" + ident, header_data, nparams);

    // Write box dimensions so we can reassemble slice.
    if (lex.size() == 0)
        return;

    utils::hdf5::Write(file_id, "lex_" + ident, (int *)&(lex[0]), lex.size());
    utils::hdf5::Write(file_id, "ley_" + ident, (int *)&(ley[0]), ley.size());
    utils::hdf5::Write(file_id, "lez_" + ident, (int *)&(lez[0]), lez.size());
    utils::hdf5::Write(file_id, "hex_" + ident, (int *)&(hex[0]), hex.size());
    utils::hdf5::Write(file_id, "hey_" + ident, (int *)&(hey[0]), hey.size());
    utils::hdf5::Write(file_id, "hez_" + ident, (int *)&(hez[0]), hez.size());
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
