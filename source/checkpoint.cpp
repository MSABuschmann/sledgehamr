#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>

#include "checkpoint.h"

namespace sledgehamr {

void Checkpoint::Write(std::string prefix) {
    amrex::Print() << "Writing checkpoint " << prefix << std::endl;

    const int nlevels = sim->finest_level + 1;
    const int noutput = sim->io_module->output.size();
    amrex::PreBuildDirectorHierarchy(prefix, "Level_", nlevels, true);

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(prefix + "/BoxArrays");
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                           std::ofstream::trunc |
                                           std::ofstream::binary);
        if( !HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        for (int lev = 0; lev < nlevels; ++lev) {
            sim->boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }

        std::string filename = prefix + "/Meta.hdf5";
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
                (double)sim->io_module->idx_checkpoints
        };
        IOModule::WriteToHDF5(file_id, "Header", header_data, nparams);

        // levels: blocking_factor, istep
        std::vector<int> bf(nlevels), istep(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            bf[lev] = sim->blocking_factor[lev][0];
            istep[lev] = sim->grid_new[lev].istep;
        }
        IOModule::WriteToHDF5(file_id, "blocking_factors", &(bf[0]), nlevels);
        IOModule::WriteToHDF5(file_id, "isteps", &(istep[0]), nlevels);

        // outputs: last id, last time written
        std::vector<int> next_id(noutput);
        std::vector<double> last_time_written(noutput);
        for (int i = 0; i < noutput; ++i) {
            next_id[i] = sim->io_module->output[i].GetNextId();
            last_time_written[i] =
                    sim->io_module->output[i].GetLastTimeWritten();
        }
        IOModule::WriteToHDF5(file_id, "next_id", &(next_id[0]), noutput);
        IOModule::WriteToHDF5(file_id, "last_time_written",
                              &(last_time_written[0]), noutput);

        H5Fclose(file_id);
    }

    // Write the MultiFab data.
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::VisMF::Write(
                sim->grid_new[lev],
                amrex::MultiFabFileFullPrefix(
                        lev, prefix, "Level_", "Cell"));
    }
}

void Checkpoint::Read(std::string prefix, int id) {
    std::string folder = prefix + "/checkpoints/" + std::to_string(id);
    amrex::Print() << "Restarting from checkpoint: " << folder << std::endl;

    const int nparams = 8;
    double header[nparams];
    std::string filename = folder + "/Meta.hdf5";
    IOModule::ReadFromHDF5(filename, {"Header"}, header);

    double time = header[0];
    int MPIranks = static_cast<int>(header[1]);
    int finest_level = static_cast<int>(header[2]);
    int dim0 = static_cast<int>(header[3]);
    int nghost = static_cast<int>(header[4]);
    int nscalars = static_cast<int>(header[5]);
    int noutput = static_cast<int>(header[6]);
    int npredefoutput = static_cast<int>(header[7]);

    if (nscalars != sim->scalar_fields.size()) {
        const char* msg = "Sledgehamr::IOModule::ReadCheckpoint: "
                          "Number of scalar fields has changed!";
        amrex::Abort(msg);
    }

    if (noutput != sim->io_module->output.size() ||
        npredefoutput != sim->io_module->idx_checkpoints) {
        const char* msg = "Sledgehamr::IOModule::ReadCheckpoint: "
                          "Number of output types changed!";
        amrex::Abort(msg);
    }

    std::string File(folder + "/BoxArrays");
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

        amrex::DistributionMapping dm {ba, amrex::ParallelDescriptor::NProcs()};
        sim->SetBoxArray(lev, ba);
        sim->SetDistributionMap(lev, dm);

        // In case nghost changed, we can create grid_old already with the new
        // value. grid_new has to be set to the old value first since ghost
        // cells are saved in the checkpoint. We change nghost for this MultiFab
        // later below.
        sim->grid_old[lev].define(ba, dm, nscalars, sim->nghost, time);
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
        // TODO
    }

    if (MPIranks != amrex::ParallelDescriptor::NProcs()) {
        amrex::Print() << "#warning: Number of MPI ranks has changed. Will "
                       << "regrid coarse level to satisfy new constraint."
                       << std::endl;
        // TODO
    }
}

void Checkpoint::GotoNextLine(std::istream& is) {
    constexpr std::streamsize bl_ignore_max {100000};
    is.ignore(bl_ignore_max, '\n');
}

}; // namespace sledgehamr
