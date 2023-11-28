#include "level_writer.h"
#include "sledgehamr_utils.h"
#include "hdf5_utils.h"

namespace sledgehamr {

LevelWriter::LevelWriter(Sledgehamr* owner, std::string prefix, int output_type)
    : sim(owner), folder(prefix), output_id(output_type) {
    DetermineSetup();
    ParseParams();
}

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
    } else if (output_id ==  sim->io_module->idx_full_box) {
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
    } else{
        std::string msg = "LevelWriter:DetermineSetup: Unknown setup!";
        amrex::Abort(msg);
    }
}

void LevelWriter::ParseParams() {
    std::string pre = "output." + name;
    amrex::ParmParse pp(pre);
    pp.query("downsample_factor", downsample_factor);
    CheckDownsampleFactor();
}

void LevelWriter::CheckDownsampleFactor() {
    if (!amrex::ParallelDescriptor::IOProcessor())
        return;

    if (!utils::IsPowerOfTwo(downsample_factor)) {
        std::string msg = "LevelWriter::CheckDownsampleFactor: Downsample "
                          "factor output." + name + " is not a power of 2!";
        amrex::Abort(msg);
    }

    for (int lev = 0; lev <= level_max; ++lev) {
        if (downsample_factor > sim->GetBlockingFactor(lev)) {
            std::string msg = "LevelWriter::CheckDownsampleFactor: Downsample "
                              "factor output." + name + " exceeds blocking "
                              "factor!";
            amrex::Abort(msg);
        }
    }
}

void LevelWriter::Write() {
    amrex::Print() << "Write " << info << ": " << folder << std::endl;

    for (int lev = level_min; lev <= level_max; ++lev) {
        // Create folder and file.
        std::string subfolder = folder + "/Level_" + std::to_string(lev);
        amrex::UtilCreateDirectory(subfolder.c_str(), 0755);

        std::string filename = subfolder + "/"
                        + std::to_string(amrex::ParallelDescriptor::MyProc())
                        + ".hdf5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        const LevelData* state = &sim->GetLevelData(lev);
        WriteSingleLevel(state, lev, file_id, "data", false);

        if (with_truncation_errors) {
            const LevelData* state = &sim->GetOldLevelData(lev);
            WriteSingleLevel(state, lev, file_id, "te", true);
        }

        H5Fclose(file_id);
    }
}

void LevelWriter::WriteSingleLevel(
        const LevelData* state, int lev, hid_t file_id, std::string ident,
        bool is_truncation_error) {
    const int ndist = is_truncation_error ? 2 : 1;
    const int grid_density = downsample_factor * ndist;

    std::vector<int> lex, hex, ley, hey, lez, hez;

    // Not performance critical. Do not use OpenMP because not thread-safe.
    for (amrex::MFIter mfi(*state, false); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();
        const auto& state_arr = state->array(mfi);

        // Get box dimensions
        int lx = bx.smallEnd(0);
        int ly = bx.smallEnd(1);
        int lz = bx.smallEnd(2);
        int hx = bx.bigEnd(0)+1;
        int hy = bx.bigEnd(1)+1;
        int hz = bx.bigEnd(2)+1;
        lex.push_back(lx);
        ley.push_back(ly);
        lez.push_back(lz);
        hex.push_back(hx);
        hey.push_back(hy);
        hez.push_back(hz);


        int dimx = (hx-lx) / grid_density;
        int dimy = (hy-ly) / grid_density;
        int dimz = (hz-lz) / grid_density;
        double volfac = 1./std::pow(downsample_factor, 3);

        // Copy data into flattened array for each scalar field.
        long len = dimx * dimy * dimz;

        // TODO Adjust output type.
        for (int f=0; f<state->nComp(); ++f) {
            std::unique_ptr<float[]> output_arr(new float[len]);
            std::fill_n(output_arr.get(), len, 0.0f);

            for (int k=lz; k<hz; ++k) {
                for (int j=ly; j<hy; ++j) {
                    for (int i=lx; i<hx; ++i) {
                        if (is_truncation_error
                            && (i%2 != 0 || j%2 != 0 || k%2 != 0))
                            continue;

                        int ind = (i-lx)/grid_density*dimy*dimz
                                + (j-ly)/grid_density*dimz
                                + (k-lz)/grid_density;

                        if (is_truncation_error) {
                            output_arr[ind] = std::max(
                                    static_cast<float>(state_arr(i,j,k,f)),
                                    output_arr[ind]);
                        } else {
                            output_arr[ind] += state_arr(i,j,k,f) * volfac;
                        }
                    }
                }
            }

            std::string dset_name = sim->GetScalarFieldName(f) + "_" + ident
                                  + "_" + std::to_string(lex.size());
            utils::hdf5::Write(file_id, dset_name, output_arr.get(), len);
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

    utils::hdf5::Write(file_id, "lex_"+ident, (int*)&(lex[0]), lex.size());
    utils::hdf5::Write(file_id, "ley_"+ident, (int*)&(ley[0]), ley.size());
    utils::hdf5::Write(file_id, "lez_"+ident, (int*)&(lez[0]), lez.size());
    utils::hdf5::Write(file_id, "hex_"+ident, (int*)&(hex[0]), hex.size());
    utils::hdf5::Write(file_id, "hey_"+ident, (int*)&(hey[0]), hey.size());
    utils::hdf5::Write(file_id, "hez_"+ident, (int*)&(hez[0]), hez.size());
}

}; // namespace sledgehamr
