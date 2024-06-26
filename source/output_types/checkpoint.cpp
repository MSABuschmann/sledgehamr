#include <filesystem>

#include <AMReX_FileSystem.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include "checkpoint.h"
#include "hdf5_utils.h"

namespace sledgehamr {

/** @brief Writes an entire checkpoint file.
 */
void Checkpoint::Write() {
    const int nlevels = sim->finest_level + 1;
    const int noutput = sim->io_module->output.size();

    for (int lev = 0; lev < nlevels; ++lev) {
        std::string level_dir = GetLevelDirName(lev);
        amrex::UtilCreateCleanDirectory(level_dir, true);
    }

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName = GetBoxArrayName();
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
        if (!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        for (int lev = 0; lev < nlevels; ++lev) {
            sim->boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }

        std::string filename = GetHeaderName();
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);

        const int nparams = 8;
        double header_data[nparams] = {
            sim->grid_new[0].t,
            (double)amrex::ParallelDescriptor::NProcs(),
            (double)sim->finest_level,
            (double)sim->dimN[0],
            (double)sim->nghost,
            (double)sim->scalar_fields.size(),
            (double)noutput,
            (double)sim->io_module->idx_checkpoints};
        utils::hdf5::Write(file_id, "Header", header_data, nparams);

        // levels: blocking_factor, istep
        std::vector<int> bf(nlevels), istep(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            bf[lev] = sim->blocking_factor[lev][0];
            istep[lev] = sim->grid_new[lev].istep;
        }
        utils::hdf5::Write(file_id, "blocking_factors", &(bf[0]), nlevels);
        utils::hdf5::Write(file_id, "isteps", &(istep[0]), nlevels);

        // outputs: last id, last time written
        std::vector<int> next_id(noutput);
        std::vector<double> last_time_written(noutput);
        for (int i = 0; i < noutput; ++i) {
            next_id[i] = sim->io_module->output[i].GetNextId();
            last_time_written[i] =
                sim->io_module->output[i].GetLastTimeWritten();
        }

        // correct checkpoint itself since it hasn't been updated yet.
        next_id[sim->io_module->idx_checkpoints]++;
        last_time_written[sim->io_module->idx_checkpoints] = sim->grid_new[0].t;

        utils::hdf5::Write(file_id, "next_id", &(next_id[0]), noutput);
        utils::hdf5::Write(file_id, "last_time_written",
                           &(last_time_written[0]), noutput);

        H5Fclose(file_id);
    }

    // Write the MultiFab data.
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::VisMF::Write(
            sim->grid_new[lev],
            amrex::MultiFabFileFullPrefix(lev, folder, "Level_", "Cell"));
    }

    amrex::ParallelDescriptor::Barrier();
}

/** @brief Reads a checkpoint header.
 * @return Whether read was successfull.
 */
bool Checkpoint::ReadHeader() {
    const int nparams = 8;
    double header[nparams];
    std::string filename = GetHeaderName();
    if (!utils::hdf5::Read(filename, {"Header"}, header)) {
        return false;
    }

    time = header[0];
    MPIranks = static_cast<int>(header[1]);
    finest_level = std::min(static_cast<int>(header[2]), sim->max_level);
    dim0 = static_cast<int>(header[3]);
    nghost = static_cast<int>(header[4]);
    nscalars = static_cast<int>(header[5]);
    noutput = static_cast<int>(header[6]);
    npredefoutput = static_cast<int>(header[7]);

    return true;
}

/** @brief Reads a checkpoint and sets data accordingly.
 */
void Checkpoint::Read() {
    if (sim->restart_sim)
        amrex::Print() << "Restarting from checkpoint: " << folder << std::endl;

    if (!ReadHeader()) {
        const char *msg = "Sledgehamr::Checkpoint::Read: "
                          "Could not find checkpoint header!";
        amrex::Abort(msg);
    }

    if (nscalars != sim->scalar_fields.size()) {
        const char *msg = "Sledgehamr::Checkpoint::Read: "
                          "Number of scalar fields has changed!";
        amrex::Abort(msg);
    }

    std::string File(GetBoxArrayName());
    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());
    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    sim->finest_level = finest_level;
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};
        sim->SetBoxArray(lev, ba);
        sim->SetDistributionMap(lev, dm);

        // In case nghost changed, we can create grid_old already with the new
        // value. grid_new has to be set to the old value first since ghost
        // cells are saved in the checkpoint. We change nghost for this MultiFab
        // later below.
        sim->grid_old[lev].define(ba, dm, nscalars, sim->nghost);
        sim->grid_new[lev].define(ba, dm, nscalars, nghost, time);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            sim->grid_new[lev],
            amrex::MultiFabFileFullPrefix(lev, folder, "Level_", "Cell"));
    }

    if (nghost != sim->nghost) {
        amrex::Print() << "#warning: Number of ghost cells has changed!\n"
                       << "checkpoint: " << nghost
                       << " vs input file: " << sim->nghost << std::endl;
        sim->level_synchronizer->ChangeNGhost(nghost);
    }

    if (MPIranks != amrex::ParallelDescriptor::NProcs()) {
        amrex::Print() << "#warning: Number of MPI ranks has changed. Will "
                       << "regrid coarse level to satisfy new constraint."
                       << std::endl;
        sim->level_synchronizer->RegridCoarse();
    }

    UpdateLevels();
}

/** @brief Moves to next line in file stream.
 * @param   is  filestream.
 */
void Checkpoint::GotoNextLine(std::istream &is) {
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

/** @brief Updates the meta data of all output modules with the meta data in
 *         checkpoint.
 */
void Checkpoint::UpdateOutputModules() {
    if (!sim->restart_sim)
        return;

    if (!ReadHeader()) {
        const char *msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find checkpoint header!";
        amrex::Abort(msg);
    }

    if (noutput != sim->io_module->output.size() ||
        npredefoutput != sim->io_module->idx_checkpoints) {
        const char *msg = "Sledgehamr::Checkpoint::Read: "
                          "Number of output types changed!";
        amrex::Abort(msg);
    }

    std::string filename = GetHeaderName();
    std::vector<int> next_id(noutput);
    std::vector<double> last_time_written(noutput);
    if (!utils::hdf5::Read(filename, {"next_id"}, &(next_id[0]))) {
        const char *msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find next_id!";
        amrex::Abort(msg);
    }

    if (!utils::hdf5::Read(filename, {"last_time_written"},
                           &(last_time_written[0]))) {
        const char *msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find last_time_written!";
        amrex::Abort(msg);
    }

    for (int i = 0; i < noutput; ++i) {
        sim->io_module->output[i].SetNextId(next_id[i]);
        sim->io_module->output[i].SetLastTimeWritten(last_time_written[i]);
    }
}

/** @brief Updates meta data of levels with that of the checkpoint.
 */
void Checkpoint::UpdateLevels() {
    std::string filename = GetHeaderName();
    std::vector<int> blocking_factor(sim->finest_level + 1);
    std::vector<int> istep(sim->finest_level + 1);
    if (!utils::hdf5::Read(filename, {"isteps"}, &(istep[0]))) {
        const char *msg = "Sledgehamr::Checkpoint::UpdateLevels: "
                          "Could not find isteps!";
        amrex::Abort(msg);
    }

    if (!utils::hdf5::Read(filename, {"blocking_factors"},
                           &(blocking_factor[0]))) {
        const char *msg = "Sledgehamr::Checkpoint::UpdateLevels: "
                          "Could not find blocking_factors!";
        amrex::Abort(msg);
    }

    // Check if blocking factor changed and react accordingly.
    for (int lev = 0; lev <= sim->finest_level; ++lev) {
        if (blocking_factor[lev] != sim->blocking_factor[lev][0]) {
            amrex::Print() << "#warning: Blocking factor on level " << lev
                           << "changed from " << blocking_factor[lev] << " to "
                           << sim->blocking_factor[lev][0] << std::endl;
        }

        // If blocking factor is now larger we cannot do a local regrid. Need
        // to do a global regrid first.
        if (blocking_factor[lev] < sim->blocking_factor[lev][0]) {
            for (int l = 0; l <= lev; ++l) {
                sim->time_stepper->local_regrid->do_global_regrid[l] = true;
            }
        }

        sim->grid_new[lev].istep = istep[lev];
    }
}

/** @brief Deletes the checkpoint.
 */
void Checkpoint::Delete() {
    amrex::Print() << "Deleting checkpoint " << folder << " ..." << std::endl;

    if (!amrex::ParallelDescriptor::IOProcessor())
        return;

    // First verify the given folder is actually a checkpoint file ...
    if (!ReadHeader()) {
        amrex::Print() << "Not a valid checkpoint! How did this happen ??"
                       << std::endl;
    }

    // Delete files/folders individually to make sure we aren't accidentally
    // deleting user data.

    std::string header = GetHeaderName();
    amrex::FileSystem::Remove(header);

    std::string ba_file = GetBoxArrayName();
    amrex::FileSystem::Remove(ba_file);

    for (int lev = 0; lev <= finest_level; ++lev) {
        std::string lev_dir = GetLevelDirName(lev);
        amrex::FileSystem::RemoveAll(lev_dir);
    }

    // Specifically use std::filesystem here and not amrex::FileSystem.
    // The amrex version does not delete emtpy folders but std::filesystem
    // does.
    std::filesystem::remove(folder);
}

}; // namespace sledgehamr
