#include "level_writer.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

/** @brief Gathers metadata.
 * @param   owner       Pointer to simulation.
 * @param   prefix      Local ouptut folder.
 * @param   output_type Id of output type, see output_id for details.
 */
LevelWriter::LevelWriter(Sledgehamr *owner, std::string prefix, int output_type)
    : sim(owner), folder(prefix), output_id(output_type) {
    DetermineSetup();
    ParseParams();
}

/** @brief Sets metadata according to output type.
 */
void LevelWriter::DetermineSetup() {
    if (output_id == sim->io_module->idx_coarse_box) {
        name = "coarse_box";
        info = "coarse level";
        level_min = 0;
        level_max = 0;
        with_truncation_errors = false;
    } else if (output_id == sim->io_module->idx_coarse_box_truncation_error) {
        name = "coarse_box_truncation_error";
        info = "coarse level truncation error estimates";
        level_min = 0;
        level_max = 0;
        with_truncation_errors = true;
    } else if (output_id == sim->io_module->idx_full_box) {
        name = "full_box";
        info = "full box (all levels)";
        level_min = 0;
        level_max = sim->GetFinestLevel();
        with_truncation_errors = false;
    } else if (output_id == sim->io_module->idx_full_box_truncation_error) {
        name = "full_box_truncation_error";
        info = "full box (all levels) truncation error estimates";
        level_min = 0;
        level_max = sim->GetFinestLevel();
        with_truncation_errors = true;
    } else {
        std::string msg = "LevelWriter:DetermineSetup: Unknown setup!";
        amrex::Abort(msg);
    }
}

/** @brief Reads relevant parameters.
 */
void LevelWriter::ParseParams() {
    std::string pre = "output." + name;
    amrex::ParmParse pp(pre);
    pp.query("downsample_factor", downsample_factor);
    pp.query("precision", precision);

    if (precision != 32 && precision != 64) {
        amrex::Print() << "Warning: Unknown precision requested for " << name
                       << ": " << precision << "\n Defaulting to 32-bit."
                       << std::endl;
    }

    CheckDownsampleFactor();
}

/** @brief Checks if downsameple factor is valid.
 */
void LevelWriter::CheckDownsampleFactor() {
    if (!amrex::ParallelDescriptor::IOProcessor())
        return;

    if (!utils::IsPowerOfTwo(downsample_factor)) {
        std::string msg = "LevelWriter::CheckDownsampleFactor: Downsample "
                          "factor output." +
                          name + " is not a power of 2!";
        amrex::Abort(msg);
    }

    for (int lev = 0; lev <= level_max; ++lev) {
        if (downsample_factor > sim->GetBlockingFactor(lev)) {
            std::string msg = "LevelWriter::CheckDownsampleFactor: Downsample "
                              "factor output." +
                              name +
                              " exceeds blocking "
                              "factor!";
            amrex::Abort(msg);
        }
    }
}

/** @brief Writes levels according to metdata constraints.
 */
void LevelWriter::Write() {
    for (int lev = level_min; lev <= level_max; ++lev) {
        // Create folder and file.
        std::string subfolder = folder + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(subfolder.c_str(), 0755);

        std::string filename =
            subfolder + "/" +
            std::to_string(amrex::ParallelDescriptor::MyProc()) + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        const LevelData *state = &sim->GetLevelData(lev);
        if (precision == 32) {
            WriteSingleLevel<float>(state, lev, file_id, "data", false);
        } else if (precision == 64) {
            WriteSingleLevel<double>(state, lev, file_id, "data", false);
        }

        if (with_truncation_errors) {
            const LevelData *state = &sim->GetOldLevelData(lev);

            if (precision == 32) {
                WriteSingleLevel<float>(state, lev, file_id, "te", true);
            } else if (precision == 64) {
                WriteSingleLevel<double>(state, lev, file_id, "te", true);
            }
        }

        H5Fclose(file_id);
    }
}

}; // namespace sledgehamr
