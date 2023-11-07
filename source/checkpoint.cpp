#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

#include "checkpoint.h"

namespace sledgehamr {

void Checkpoint::Write(std::string prefix) {
    amrex::Print() << "Writing checkpoint " << prefix << std::endl;

    const int nlevels = sim->finest_level + 1;
    const int noutput = sim->io_module->output.size();

    for(int lev = 0; lev < nlevels; ++lev) {
        amrex::Print() << "Create : " << prefix + "/Level_" + std::to_string(lev) << std::endl;
        amrex::UtilCreateCleanDirectory(
                prefix + "/Level_" + std::to_string(lev), true);
    }

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

        std::string filename = GetHeaderName(prefix);
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

        // correct checkpoint itself since it hasn't been updated yet.
        next_id[sim->io_module->idx_checkpoints]++;
        last_time_written[sim->io_module->idx_checkpoints] = sim->grid_new[0].t;

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

bool Checkpoint::ReadHeader(std::string folder) {
    const int nparams = 8;
    double header[nparams];
    std::string filename = GetHeaderName(folder);
    if (!IOModule::ReadFromHDF5(filename, {"Header"}, header)) {
        return false;
    }

    double time = header[0];
    int MPIranks = static_cast<int>(header[1]);
    int finest_level = static_cast<int>(header[2]);
    int dim0 = static_cast<int>(header[3]);
    int nghost = static_cast<int>(header[4]);
    int nscalars = static_cast<int>(header[5]);
    int noutput = static_cast<int>(header[6]);
    int npredefoutput = static_cast<int>(header[7]);
    return true;
}

void Checkpoint::Read(std::string prefix, int id) {
    std::string folder = prefix + "/checkpoints/" + std::to_string(id);
    Read(folder);
}

void Checkpoint::Read(std::string folder) {
    if (sim->restart_sim)
        amrex::Print() << "Restarting from checkpoint: " << folder << std::endl;

    if (!ReadHeader(folder)) {
        const char* msg = "Sledgehamr::Checkpoint::Read: "
                          "Could not find checkpoint header!";
        amrex::Abort(msg);
    }

    if (nscalars != sim->scalar_fields.size()) {
        const char* msg = "Sledgehamr::Checkpoint::Read: "
                          "Number of scalar fields has changed!";
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
        ChangeNGhost(nghost);
    }

    if (MPIranks != amrex::ParallelDescriptor::NProcs()) {
        amrex::Print() << "#warning: Number of MPI ranks has changed. Will "
                       << "regrid coarse level to satisfy new constraint."
                       << std::endl;
        RegridCoarse();
    }

    UpdateLevels(folder);
}

void Checkpoint::GotoNextLine(std::istream& is) {
    constexpr std::streamsize bl_ignore_max {100000};
    is.ignore(bl_ignore_max, '\n');
}

void Checkpoint::ChangeNGhost(int new_nghost) {
    // Create new LevelData.
    const int lev                        = 0;
    LevelData& ld_old                    = sim->grid_new[lev];
    const amrex::BoxArray& ba            = ld_old.boxArray();
    const amrex::DistributionMapping& dm = ld_old.DistributionMap();
    const int ncomp                      = ld_old.nComp();
    const double time                    = ld_old.t;
    const amrex::Geometry& geom          = sim->geom[lev];

    // Allocate and fill.
    LevelData ld_new(ba, dm, ncomp, new_nghost, time);

    amrex::CpuBndryFuncFab bndry_func(nullptr); 
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(
            geom, sim->level_synchronizer->bcs, bndry_func);

    amrex::Vector<amrex::MultiFab*> smf{static_cast<amrex::MultiFab*>(&ld_old)};
    amrex::Vector<double> stime{time};
    amrex::MultiFab& mf = ld_new;

    amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp, geom,
                                physbc, 0);

    // Swap and update.
    std::swap(sim->grid_new[lev], ld_new);
    sim->nghost = new_nghost;
}

void Checkpoint::RegridCoarse() {
    const int lev               = 0;
    LevelData& ld_old           = sim->grid_new[lev];
    const int ncomp             = ld_old.nComp();
    const double time           = ld_old.t;
    const amrex::Geometry& geom = sim->geom[lev];

    sim->grid_old[lev].clear();

    // New layout.
    amrex::BoxArray ba = ld_old.boxArray();
    amrex::Box bx = ba.minimalBox();
    ba = amrex::BoxArray(bx);
    sim->ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs()); 
    amrex::DistributionMapping dm(ba, amrex::ParallelDescriptor::NProcs());

    // Allocate and fill.
    LevelData ld_new(ba, dm, ncomp, sim->nghost, time);

    amrex::CpuBndryFuncFab bndry_func(nullptr); 
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(
            geom, sim->level_synchronizer->bcs, bndry_func);

    amrex::Vector<amrex::MultiFab*> smf{static_cast<amrex::MultiFab*>(&ld_old)};
    amrex::Vector<double> stime{time};
    amrex::MultiFab& mf = ld_new;

    amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp, geom,
                                physbc, 0);

    std::swap(sim->grid_new[lev], ld_new);
    sim->grid_old[lev] = LevelData(ba, dm, ncomp, sim->nghost, time);
    sim->SetBoxArray(lev, ba);
    sim->SetDistributionMap(lev, dm);
}

void Checkpoint::UpdateOutputModules(std::string prefix, int id) {
    std::string folder = prefix + "/checkpoints/" + std::to_string(id);
    UpdateOutputModules(folder);
}

void Checkpoint::UpdateOutputModules(std::string folder) {
    if (!sim->restart_sim)
        return;

    if (!ReadHeader(folder)) {
        const char* msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find checkpoint header!";
        amrex::Abort(msg);
    }

    if (noutput       != sim->io_module->output.size() ||
        npredefoutput != sim->io_module->idx_checkpoints) {
        const char* msg = "Sledgehamr::Checkpoint::Read: "
                          "Number of output types changed!";
        amrex::Abort(msg);
    }

    std::string filename = GetHeaderName(folder);
    std::vector<int> next_id(noutput);
    std::vector<double> last_time_written(noutput);
    if (!IOModule::ReadFromHDF5(filename, {"next_id"}, &(next_id[0]))) {
        const char* msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find next_id!";
        amrex::Abort(msg);
    }

    if (!IOModule::ReadFromHDF5(filename, {"last_time_written"},
                                &(last_time_written[0]))) {
        const char* msg = "Sledgehamr::Checkpoint::UpdateOutputModules: "
                          "Could not find last_time_written!";
        amrex::Abort(msg);
    }

    for (int i = 0; i < noutput; ++i) {
        sim->io_module->output[i].SetNextId(next_id[i]);
        sim->io_module->output[i].SetLastTimeWritten(last_time_written[i]);
    }
}

void Checkpoint::UpdateLevels(std::string folder) {
    std::string filename = GetHeaderName(folder);
    std::vector<int> blocking_factor(sim->finest_level+1);
    std::vector<int> istep(sim->finest_level+1);
    if (!IOModule::ReadFromHDF5(filename, {"isteps"}, &(istep[0]))) {
        const char* msg = "Sledgehamr::Checkpoint::UpdateLevels: "
                          "Could not find isteps!";
        amrex::Abort(msg);
    }

    if (!IOModule::ReadFromHDF5(filename, {"blocking_factors"},
                            &(blocking_factor[0]))) {
        const char* msg = "Sledgehamr::Checkpoint::UpdateLevels: "
                          "Could not find blocking_factors!";
        amrex::Abort(msg);
    }

    // Check if blocking factor changed and react accordingly.
    for(int lev = 0; lev <= sim->finest_level; ++lev) {
        if (blocking_factor[lev] != sim->blocking_factor[lev][0]) {
            amrex::Print() << "#warning: Blocking factor on level " << lev
                           << "changed from " << blocking_factor[lev]
                           << " to " << sim->blocking_factor[lev][0]
                           << std::endl;
        }

        // If blocking factor is now larger we cannot do a local regrid. Need
        // to do a global regrid first.
        if (blocking_factor[lev] < sim->blocking_factor[lev][0]) {
            for(int l = 0; l <= lev; ++l) {
                sim->time_stepper->local_regrid->do_global_regrid[l] = true;
            }
        }

        sim->grid_new[lev].istep = istep[lev];
    }
}

}; // namespace sledgehamr
